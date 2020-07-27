'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import sys
import multiprocessing as mp
import os
import pysam
import time
import subprocess
import statistics
import itertools

# Internal
import log
import unix
import databases
import formats
import bamtools
import structures
import events
import clusters
import output
import annotation
import bkp
import filters
import alignment
import gRanges

#from modules.callers import SV_caller
from callers import SV_caller

## FUNCTIONS ##

## CLASSES ##
class SV_caller_long(SV_caller):
    '''
    Structural variation (SV) caller for long read sequencing data
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

    def call(self):
        '''
        Search for structural variants (SV) genome wide or in a set of target genomic regions
        '''

        ### 1. Create SV clusters 
        msg = '1. Create SV clusters'
        log.header(msg)
        allMetaclusters = self.make_clusters()

        ### 2. Annotate SV clusters intervals  
        msg = '2. Annotate SV clusters intervals'
        log.header(msg)

        # Load annotations
        self.load_annotations()

        # Empty list harcoded, because it is not needed por viruses and it crashes otherwise.
        self.annotations['TRANSDUCTIONS'] = []

        # Create output directory
        annotDir = self.outDir + '/ANNOT/'
        unix.mkdir(annotDir)
        
        # Define annotation steps
        steps = ['REPEAT']

        if self.confDict['annovarDir'] is not None:
            steps.append('GENE')

        # For each cluster type
        for SV_type in allMetaclusters:
            
            metaclusters = allMetaclusters[SV_type]
            annotation.annotate(metaclusters, steps, self.annotations, self.confDict['annovarDir'], annotDir)

        # Remove annotation directory
        unix.rm([annotDir])

        ### 3. Determine what type of sequence has been inserted for INS metaclusters
        msg = '3. Determine what type of sequence has been inserted for INS metaclusters'
        log.header(msg)

        # Create output directory
        outDir = self.outDir + '/INS_TYPE/'
        unix.mkdir(outDir)

        if 'INS' in allMetaclusters:

            ## Infer insertion type
            clusters.INS_type_metaclusters(allMetaclusters['INS'], self.reference, self.annotations, 1, self.confDict['viralDb'], outDir)

        # Remove output directory
        unix.rm([outDir])
            
        ### 4. Resolve structure for solo, partnered and orphan transductions
        msg = '4. Resolve structure for solo, partnered and orphan transductions'
        log.header(msg)
        
        if 'INS' in allMetaclusters:
            consensus = self.refDir + '/consensusDb.fa'
            transduced = self.refDir + '/transducedDb.fa.masked'

            # Create output directory
            outDir = self.outDir + '/STRUCTURE/'
            unix.mkdir(outDir)

            # Structure inference
            allMetaclusters['INS'] = clusters.structure_inference_parallel(allMetaclusters['INS'], consensus, transduced, self.confDict['transductionSearch'], self.confDict['processes'], outDir)

            # Remove output directory
            unix.rm([outDir])
        
        ### 5. Identify BND junctions
        msg = '5. Identify BND junctions'
        log.header(msg)
        allJunctions = []

        if 'BND' in allMetaclusters:
        
            ### Search for repeat, transduction or viral bridges
            # Create output directory
            outDir = self.outDir + '/BND_JUNCTIONS/'
            unix.mkdir(outDir)   
            allMetaclusters['BND'] = clusters.search4bridges_metaclusters_parallel(allMetaclusters['BND'], 10000, 80, self.confDict['minReads'], 25, self.annotations, self.refDir, self.confDict['viralDb'], self.confDict['processes'], outDir)

            ### Search for BND junctions
            allJunctions = clusters.search4junctions_metaclusters(allMetaclusters['BND'], self.refLengths, self.confDict['processes'], self.confDict['minReads'], 25, self.reference, self.refDir, self.confDict['viralDb'], outDir)

            # NOTE 2020: New June 2020. For keeping those BNDs without pair
            
            for metaclusterBND in allMetaclusters['BND']:
                if metaclusterBND not in allJunctions:
                    if 'solo-BND' not in allMetaclusters:
                        allMetaclusters['solo-BND'] = []
                    allMetaclusters['solo-BND'].append(metaclusterBND)
            
            if allMetaclusters['solo-BND']:
                clusters.soloBND_type_metaclusters(allMetaclusters['solo-BND'], self.confDict, self.reference, self.refLengths, self.refDir, self.confDict['transductionSearch'], 1, self.confDict['viralDb'], outDir)
            


            # Remove output directory
            unix.rm([outDir])

        ### 6. Assess MEI novelty
        msg = '6. Assess MEI novelty'
        log.header(msg)

        if self.confDict['germlineMEI'] is not None:
            steps = ['GERMLINE-MEI']
            annotation.annotate(allMetaclusters['INS'], steps, self.annotations, self.confDict['annovarDir'], annotDir)

        ### 7. Apply second round of filtering 
        msg = '7. Apply second round of filtering'
        log.header(msg)
        filters2Apply = ['PERC-RESOLVED']
        metaclustersPass, metaclustersFailed = filters.filter_metaclusters(allMetaclusters, filters2Apply, self.confDict)
                
        ### 8. Report SV calls into output files
        msg = '8. Report SV calls into output files'
        log.header(msg)
        
        ##  8.1 Report INS
        if 'INS' in metaclustersPass:
            outFileName = 'INS_MEIGA.PASS'
            output.INS2VCF(metaclustersPass['INS'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        if 'INS' in metaclustersFailed:
            outFileName = 'INS_MEIGA.FAILED.2'
            output.INS2VCF(metaclustersFailed['INS'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        ## 8.2 Report BND junctions
        if allJunctions:
            outFileName = 'BND_MEIGA.tsv'
            output.write_junctions(allJunctions, outFileName, self.outDir)

        ## 8.2 Report solo-BND junctions
        if allMetaclusters['solo-BND']:
            outFileName = 'soloBND_MEIGA.tsv'
            output.INS2VCF_junction(allMetaclusters['solo-BND'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)
            # TODO 2020: Hacer un write especifico.
            for metasolo in allMetaclusters['solo-BND']:
                print ('solo-BND ' + str(metasolo.beg) + ' '+ str(metasolo.SV_features))
        
        
    def make_clusters(self):
        '''
        Search for structural variant (SV) clusters 

        Output:
            1. metaclustersPass: dictionary containing one key per SV type and the list of metaclusters identified of each given SV type
        '''
        ### 1. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])
        
        ### 2. Search for SV clusters in each bin ##
        # Create output directory
        unix.mkdir(self.outDir + '/CLUSTER/')

        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        metaclustersPassList, metaclustersFailedList = zip(*pool.starmap(self.make_clusters_bin, bins))
        # NOTE 2020: New 2020
        #metaclustersPass = pool.starmap(self.make_clusters_bin, bins)
        pool.close()
        pool.join()

        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        ### 3. Collapse metaclusters in a single dict and report metaclusters that failed filtering
        # NOTE 2020: In 2020:
        #metaclustersPass = structures.merge_dictionaries(metaclustersPass)

        metaclustersPass = structures.merge_dictionaries(metaclustersPassList)
        metaclustersFailed = structures.merge_dictionaries(metaclustersFailedList)
        
        if 'INS' in metaclustersFailed:
            outFileName = 'INS_MEIGA.FAILED.1.tsv'
            output.write_INS(metaclustersFailed['INS'], outFileName, self.outDir)
   
                
        return metaclustersPass

    def make_clusters_bin(self, ref, beg, end):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''
        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)
        start = time.time()

        binDir = self.outDir + '/CLUSTER/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:
            eventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        SV_types = sorted(eventsDict.keys())
        counts = [str(len(eventsDict[SV_type])) for SV_type in SV_types]
       
        step = 'COLLECT'
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)
        
        ## 2. Organize all the SV events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior metaclustering'
        log.step(step, msg)

        ## Define bin database sizes 
        minBinSize = min([self.confDict['maxInsDist'], self.confDict['maxBkpDist']])
        binSizes = [minBinSize, 1000, 10000, 100000, 1000000]

        ## Create bins
        eventsBinDb = structures.create_bin_database_interval(ref, beg, end, eventsDict, binSizes)

        ## 3. Group events into clusters ##
        step = 'CLUSTERING'
        msg = 'Group events into clusters' 
        log.step(step, msg)
        clustersBinDb = clusters.create_clusters(eventsBinDb, self.confDict)
        
        msg = 'Number of created clusters: ' + str(clustersBinDb.nbEvents()[0])
        log.step(step, msg)

        ## 4. Polish clusters ##
        step = 'POLISH'
        msg = 'Polish SV clusters' 
        log.step(step, msg)
        clusters.polish_clusters(clustersBinDb, self.confDict['minClusterSize'])

        ## 5. Group events into metaclusters ##
        step = 'META-CLUSTERING'
        msg = 'Group events into metaclusters' 
        log.step(step, msg)

        metaclusters = clusters.create_metaclusters(clustersBinDb)        
        msg = 'Number of created metaclusters: ' + str(len(metaclusters)) 
        log.step(step, msg)

        ## 6. Infer structural variant type ##
        step = 'SV-TYPE'
        msg = 'Infer structural variant type' 
        log.step(step, msg)

        outDir = binDir + '/SV_TYPE/' 
        metaclustersSVType = clusters.SV_type_metaclusters(metaclusters, self.confDict['minINDELlen'], self.confDict['technology'], outDir)
        
        # Do cleanup
        unix.rm([outDir])

        ## 7. Filter metaclusters ##
        step = 'FILTER'
        msg = 'Filter out metaclusters' 
        log.step(step, msg)
        filters2Apply = ['MIN-NBREADS', 'MAX-NBREADS', 'CV', 'SV-TYPE']
        # NOTE 2020: New 2020
        metaclustersSVType, metaclustersSVTypeFailed = filters.filter_metaclusters(metaclustersSVType, filters2Apply, self.confDict)
        #metaclustersSVType = filters.filter_metaclusters(metaclustersSVType, filters2Apply, self.confDict)[0]

        ## 8. Generate consensus event for SV metaclusters ##
        step = 'CONSENSUS'
        msg = 'Generate consensus event for SV metaclusters' 
        log.step(step, msg)

        targetSV = ['INS']
        outDir = binDir + '/CONSENSUS/' 
        clusters.create_consensus(metaclustersSVType, self.confDict, self.reference, targetSV, outDir)       

        ## 9. Lighten up metaclusters  ##
        ## Metaclusters passing all the filters
        clusters.lighten_up_metaclusters(metaclustersSVType)
        
        ## Filtered metaclusters
        clusters.lighten_up_metaclusters(metaclustersSVTypeFailed)

        # Do cleanup
        unix.rm([outDir, binDir])

        ## Print time taken to process bin
        end = time.time()
        time_taken = end - start
        msg = 'SV calling in bin: ' + binId + ' finished in ' + str(time_taken)
        log.info(msg)

        return metaclustersSVType, metaclustersSVTypeFailed
        #return metaclustersSVType