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
import retrotransposons

## FUNCTIONS ##

## CLASSES ##
class SV_caller():
    '''
    Structural variation (SV) caller 
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        self.mode = mode
        self.bam = bam
        self.normalBam = normalBam
        self.reference = reference
        self.refDir = refDir
        self.confDict = confDict
        self.outDir = outDir
        self.repeatsBinDb = None

        ## Compute reference lengths
        self.refLengths = bamtools.get_ref_lengths(self.bam)

    def minimap2_index(self):
        '''
        Return path to minimap2 index file
        '''
        index = os.path.splitext(self.reference)[0] + '.mmi' 

        return index

    def load_annotations(self):
        '''
        Load set of annotations into bin databases. Set 'annotation' attribute with one key per annotation type
        and bin databases as values
        '''
        # NOTE MERGE SR2020: For SV_caller_short this load is managed inside the caller
        # SONIA: why??
        if self.confDict['technology'] != 'ILLUMINA':
            annotDir = self.outDir + '/LOAD_ANNOT/'
            unix.mkdir(annotDir)
            annotations2load = ['REPEATS']
    
            if self.confDict['transductionSearch']:    
                annotations2load.append('TRANSDUCTIONS')
    
            if True: # at one point include flag for pseudogene search
                annotations2load.append('EXONS')
    
            if self.confDict['germlineMEI'] is not None:
                annotations2load.append('GERMLINE-MEI')
    
            self.annotations = annotation.load_annotations(annotations2load, self.refLengths, self.refDir, self.confDict['germlineMEI'], self.confDict['processes'], annotDir)        
            unix.rm([annotDir])


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
            '''
            for metaclusterBND in allMetaclusters['BND']:
                if metaclusterBND not in allJunctions:
                    if 'solo-BND' not in allMetaclusters:
                        allMetaclusters['solo-BND'] = []
                    allMetaclusters['solo-BND'].append(metaclusterBND)
            
            if allMetaclusters['solo-BND']:
                clusters.soloBND_type_metaclusters(allMetaclusters['solo-BND'], self.confDict, self.reference, self.refLengths, self.refDir, self.confDict['transductionSearch'], 1, self.confDict['viralDb'], outDir)
            '''


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
            outFileName = 'INS_MEIGA.FAILED'
            output.INS2VCF(metaclustersFailed['INS'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        ## 8.2 Report BND junctions
        if allJunctions:
            outFileName = 'BND_MEIGA.tsv'
            output.write_junctions(allJunctions, outFileName, self.outDir)

        ## 8.2 Report solo-BND junctions
        '''
        if allMetaclusters['solo-BND']:
            outFileName = 'soloBND_MEIGA.tsv'
            output.INS2VCF_junction(allMetaclusters['solo-BND'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)
            # TODO 2020: Hacer un write especifico.
            for metasolo in allMetaclusters['solo-BND']:
                print ('solo-BND ' + str(metasolo.beg) + ' '+ str(metasolo.SV_features))
        '''
        
        
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
        metaclustersPass = pool.starmap(self.make_clusters_bin, bins)
        pool.close()
        pool.join()

        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        ### 3. Collapse metaclusters in a single dict
        metaclustersPass = structures.merge_dictionaries(metaclustersPass)       
                
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
        metaclustersSVType = filters.filter_metaclusters(metaclustersSVType, filters2Apply, self.confDict)[0]

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
        
        # Do cleanup
        unix.rm([outDir, binDir])

        ## Print time taken to process bin
        end = time.time()
        time_taken = end - start
        msg = 'SV calling in bin: ' + binId + ' finished in ' + str(time_taken)
        log.info(msg)

        return metaclustersSVType

# Multiprocessing lock as global variable. Useful for safely writting from different processes to the same output file.
def init(l):
    global lock
    lock = l

class SV_caller_short(SV_caller):
    '''
    Structural variation (SV) caller for Illumina short read sequencing data
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

        self.repeatsBinDb = None
        self.viralSeqs = {}
        
    def call(self):
        '''
        Search for integrations genome wide or in a set of target genomic regions
        '''
        # NOTE SR: if 'ME' in self.confDict['targetINT2Search'] annotRepeats is not and option
        if 'ME' in self.confDict['targetINT2Search']:
            annotDir = self.outDir + '/ANNOT/'
            unix.mkdir(annotDir)
            annotations2load = ['REPEATS']
    
            if self.confDict['transductionSearch']:    
                annotations2load.append('TRANSDUCTIONS')
    
            # if True: # at one point include flag for pseudogene search
            #     annotations2load.append('EXONS')
    
            if self.confDict['germlineMEI'] is not None:
                annotations2load.append('GERMLINE-MEI')
    
            self.annotations = annotation.load_annotations(annotations2load, self.refLengths, self.refDir, self.confDict['germlineMEI'], self.confDict['processes'], annotDir)        

            
        else:
            self.annotations = {}
            self.annotations['REPEATS'], self.annotations['TRANSDUCTIONS'] = None, None

        ### 1. Create integration clusters 
        msg = '1. Create integration clusters. PID: ' + str(os.getpid())
        log.header(msg)      
        '''
        EXPLANATION
        Tuple of lists of lists -> metaclustersListofLists ([], [], [], [], [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>',
        '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}]], [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}]], [], [], [])
        '''
        metaclustersListofLists = self.make_clusters()

        # Flat metaclustersList
        '''
        EXPLANATION
        Lists of lists -> metaclustersList [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, 56715108, '.', '<INS>', '.', 'PASS', {INFO_Dict}]]
        '''
        metaclustersList = list(itertools.chain(*metaclustersListofLists))

        # NOTE SR: Perform repeats_annotation for both, MEIs and VIRUSES.
        #if 'ME' in self.confDict['targetINT2Search']:

        # If ME analysis is done, this is loaded at the beggining
        # DESILENCE
        if 'ME' not in self.confDict['targetINT2Search'] and self.confDict['annotRepeats']:
            annotDir = self.outDir + '/ANNOT/'
            self.annotations = annotation.load_annotations(['REPEATS', 'TRANSDUCTIONS'], self.refLengths, self.refDir, self.confDict['germlineMEI'], self.confDict['processes'], annotDir)

        # NOTE SR: if 'ME' in self.confDict['targetINT2Search'] annotRepeats is not and option
        if 'ME' in self.confDict['targetINT2Search'] or self.confDict['annotRepeats']:
            ### 2. Annotate SV clusters intervals
            msg = '2. Annotate SV clusters intervals'
            log.header(msg)

            ## 5. Check if annotated retrotransposon on the reference genome at cluster intervals ##
            # COMMENT: This is temporary and will be incorporated into the filtering function at one point
            step = 'ANNOTATE-REPEATS'
            msg = 'Check if annotated retrotransposon on the reference genome at cluster intervals'
            log.step(step, msg)

            ## Annotate
            buffer = 100
            annotation.repeats_annotation_lighter(metaclustersList, self.annotations['REPEATS'], buffer)
            
            ## 6. Perform gene-based annotation with ANNOVAR of discordant read pair clusters ##
            # Do gene-based annotation step if enabled
            # Remove annotation directory
            unix.rm([annotDir])
            
        if self.confDict['annovarDir'] is not None:
            annotDir = self.outDir + '/ANNOT/'
            step = 'ANNOTATE'
            msg = 'Perform gene-based annotation with ANNOVAR of discordant read-pair clusters'
            log.step(step, msg)

            ## Annotate
            annotation.gene_annotation_lighter(metaclustersList, self.confDict['annovarDir'], annotDir)

            # Remove annotation directory
            unix.rm([annotDir])

        # TODO SR: Think if is worth it to make viral db inside MEIGA (headers, etc)
        # TODO SR: Index viral db
        ## 1.2 Create and index viral database
        #self.viralDb, self.viralDbIndex = databases.buildVirusDb(self.refDir, dbDir)
        
  
        # Report integrations calls into output VCF files
        outFileName = 'metaclusters.PASS'
        output.INS2VCF_SR(metaclustersList, self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], self.confDict['VCFInfoFields'], self.confDict['VCFREF'], outFileName, self.outDir)


    def make_clusters(self):

        ### If viruses option is selected, collect read name and sequence of discordant low quality reads from all bam refs ##
        if 'VIRUS' in self.confDict['targetINT2Search']:
            # TEMP SR: DESILENCE
            # Make genomic bins
            bins = bamtools.makeGenomicBins(self.bam, self.confDict['binSize'], None)
            
            l = mp.Lock()

            # Collect read name and sequence of discordant low quality reads from all bam refs
            collectVirusDir = self.outDir + '/COLLECT_VIRUS'
            unix.mkdir(collectVirusDir)

            pool = mp.Pool(processes=self.confDict['processes'], initializer=init, initargs=(l,))
            pool.starmap(self.callCollectSeq, bins)
            pool.close()
            pool.join()

            # If normal bam is present, collect also its reads
            if self.mode == "PAIRED":
                bins = bamtools.makeGenomicBins(self.normalBam, self.confDict['binSize'], None)

                # TODO SR: Pass more arguments
                pool = mp.Pool(processes=self.confDict['processes'], initializer=init, initargs=(l,))
                pool.starmap(self.callCollectSeqNormal, bins)
                pool.close()
                pool.join()
        
            # Filter viral discodants and store them in a dictionary. Algo make a fasta file with those viral sequences from viralDB found in bam.
            self.viralSeqs, self.identDbPath = virus.find_virus_discordants(self.bam, self.normalBam, self.confDict['viralDb'], self.confDict['komplexityThreshold'], self.confDict['minTotalMatchVirus'], self.confDict['minParcialMatchVirus'], self.confDict['maxMatchCheckMAPQVirus'], self.confDict['minMAPQVirus'], self.confDict['maxBasePercVirus'], self.confDict['minLccVirus'], self.confDict['processes'], self.outDir)
            # TEMP SR: Remove allfastas
            #unix.rm([collectVirusDir])
            
        ### 1. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 2. Search for SV clusters in each bin ##
        # Create output directory
        unix.mkdir(self.outDir + '/CLUSTER/')

        # Genomic bins will be distributed into X processes
        # TODO: mirar que pasa cuando tienes 2 dictionarios
        pool = mp.Pool(processes=self.confDict['processes'])
        metaclustersPassListofLists, metaclustersFailedListofLists = zip(*pool.starmap(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        ### 3. Clean up
        if 'VIRUS' in self.confDict['targetINT2Search']:
            del self.viralSeqs
            
            if not self.confDict['keepIdentDb']:
                unix.rm([self.identDbPath])
            
        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        '''
        EXPLANATION
        Tuple of lists of lists -> metaclustersFailedListofLists ([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], {}, {}, [], [], [], [],
        [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>',
        '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}]], [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}]], [], [], [])
        '''
        if metaclustersFailedListofLists:
            # Flat metaclustersFailedList
            metaclustersFailedList = list(itertools.chain(*metaclustersFailedListofLists))
            if metaclustersFailedList:
                outFileName = 'metaclusters.FAILED'
                '''
                EXPLANATION
                Lists of lists -> metaclustersFailedList [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, 56715108, '.', '<INS>', '.', 'PASS', {INFO_Dict}]]
                '''

                # Write VCF output
                output.INS2VCF_SR(metaclustersFailedList, self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], self.confDict['VCFInfoFields'], self.confDict['VCFREF'], outFileName, self.outDir)

        return metaclustersPassListofLists
    
    def callCollectSeq(self, ref, binBeg, binEnd):
        '''
        Wrapper to call function for collecting read name and sequence of discordant low quality reads from all TEST bam refs
        '''
        bamtools.collectDiscodantsLowMAPQSeq(ref, binBeg, binEnd, self.bam, self.confDict['discordantMatesMaxMAPQ'], self.confDict['discordantMatesCheckUnmapped'], self.confDict['discordantMatesSupplementary'], self.confDict['discordantMatesMaxBasePerc'], self.confDict['discordantMatesMinLcc'], self.outDir)

    def callCollectSeqNormal(self, ref, binBeg, binEnd):
        '''
        Wrapper to call function for collecting read name and sequence of discordant low quality reads from all NORMAL bam refs
        '''
        bamtools.collectDiscodantsLowMAPQSeq(ref, binBeg, binEnd, self.normalBam, self.confDict['discordantMatesMaxMAPQ'], self.confDict['discordantMatesCheckUnmapped'], self.confDict['discordantMatesSupplementary'], self.confDict['discordantMatesMaxBasePerc'], self.confDict['discordantMatesMinLcc'], self.outDir)

    def make_clusters_bin(self, ref, beg, end):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''
        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'INSERTION calling in bin: ' + binId + ' PID: ' + str(os.getpid())
        log.subHeader(msg)
        start = time.time()

        binDir = self.outDir + '/CLUSTER/' + binId
        unix.mkdir(binDir)

        ## 1. Search for integration candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            discordants = bamtools.collectOnlyDISCORDANT(ref, beg, end, self.bam, self.confDict, None, False)

        # b) Paired sample mode (tumour & matched normal)
        else:
            discordants = bamtools.collectOnlyDISCORDANT_paired(ref, beg, end, self.bam, self.normalBam, self.confDict, False)

        counts = str(len(discordants))
        step = 'COLLECT'
        msg = 'Number of DISCORDANT events in bin ' + binId + ': ' + counts + '. PID: ' + str(os.getpid())
        log.step(step, msg)
        
        ## 2. Discordant read pair identity ##
        ## Determine identity
        discordantsIdentity = events.determine_discordant_identity(discordants, self.annotations['REPEATS'], self.annotations['TRANSDUCTIONS'], self.bam, None, binDir, self.confDict['targetINT2Search'], self.viralSeqs)

        del discordants

        step = 'IDENTITY'
        SV_types = sorted(discordantsIdentity.keys())
        counts = [str(len(discordantsIdentity[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events per identity in bin ' + binId + '(' + ','.join(SV_types) + '): ' + ','.join(counts) + '. PID: ' + str(os.getpid())
        log.step(step, msg)

        ## 3. Group discordant read pairs into metaclusters based on their mate identity ##
        # NOTE SR: discordantsIdentity list is already sorted by beg (because is the way that pysam reads the bam file).
        metaclusters = clusters.create_discordantClusters(discordantsIdentity, self.confDict['minClusterSize'], self.confDict['equalOrientBuffer'], self.confDict['oppositeOrientBuffer'], self.confDict['libraryReadLength'])

        del discordantsIdentity

        step = 'META-CLUSTERING'
        counts = len(metaclusters)
        msg = 'Number of created metaclusters clusters in bin ' + binId + ': ' + str(counts) + ' . PID: ' + str(os.getpid())
        log.step(step, msg)

        # Return if no DISCODANT metaclusters found.
        if counts == []:
            unix.rm([binDir])
            metaclustersSVType, metaclustersSVTypeFailed = {}, {}
            return metaclustersSVType, metaclustersSVTypeFailed

        # TODO SR: ANNOTATE-REPEATS step: Desilence and put in the right place (now it is repeated in two different places) in case we want to analyse RT. If not, decide if it is neccessary or not.
        '''
        ## 5. Check if annotated retrotransposon on the reference genome at cluster intervals ##
        # COMMENT: This is temporary and will be incorporated into the filtering function at one point
        step = 'ANNOTATE-REPEATS'
        msg = 'Check if annotated retrotransposon on the reference genome at cluster intervals'
        log.step(step, msg)

        ## Create a list containing all discordant read pair events:
        allDiscordantClusters = []

        for eventType in discordantClustersDict.keys():
            allDiscordantClusters.extend(discordantClustersDict[eventType])

        ## Annotate
        buffer = 100
        annotation.repeats_annotation(allDiscordantClusters, self.repeatsBinDb, buffer)
        
        ## 6. Perform gene-based annotation with ANNOVAR of discordant read pair clusters ##
        # Do gene-based annotation step if enabled
        if self.confDict['annovarDir'] is not None:

            step = 'ANNOTATE'
            msg = 'Perform gene-based annotation with ANNOVAR of discordant read-pair clusters'
            log.step(step, msg)

            ## Annotate
            annotDir = binDir + '/ANNOT/' 
            annotation.gene_annotation(allDiscordantClusters, self.confDict['annovarDir'], annotDir)
        '''

        ## 4. Filter metaclusters before adding clipping reads ##

        filters2Apply = {}
        # Filters to apply to non-identified clusters
        filters2Apply['GENERIC'] = self.confDict['filtersBfClip']
        # Filters to apply to VIRUS clusters
        filters2Apply['VIRUS'] = self.confDict['filtersBfClip']
        # TODO SR: Add here proper filters for MEs clusters
        filters2Apply['ME']  = self.confDict['filtersBfClip']
        metaclustersPass1, metaclustersFailed1 = filters.filter_metaclusters_SR(metaclusters, filters2Apply, self.confDict, self.bam)

        # Remove those metaclusters that didnt pass the filteres if the option for not printting them is selected
        if not self.confDict['printFiltered']:
            metaclustersFailed1 = []

        step = 'FILTER METACLUSTERS BEFORE ADDING CLIPPING READS'
        msg = 'Filter out metaclusters. Nb PASS metaclusters: '+ str(len(metaclustersPass1)) +'. Nb NOT PASS: ' + str(len(metaclustersFailed1)) + '. PID: ' + str(os.getpid())
        log.step(step, msg)

        ## 5. Analyse metaclusters features and add supporting clipping reads ##
        # TODO SR: Now adding clipping step is better than before. Even though maybe it is not neccessary, it would be great to choose in a better way which clipping we should add to the metacluster.

        # Make viral identities DB
        if 'VIRUS' in self.confDict['targetINT2Search'] and 'ME' in self.confDict['targetINT2Search']:

            filenames = [self.confDict['MEDB'], self.identDbPath]
            mergedDbPath = self.outDir + '/ME_VIRUS_DB.fa'
            with open(mergedDbPath, 'w') as outfile:
                for fname in filenames:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            
            outfile.close()

            bkp.analyzeMetaclusters(metaclustersPass1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, mergedDbPath)
            if self.confDict['analyseFiltered']:
                bkp.analyzeMetaclusters(metaclustersFailed1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, mergedDbPath)

        # If only VIRUS are being analysed, pick viral identities db which is created internally.
        elif 'VIRUS' in self.confDict['targetINT2Search'] and not 'ME' in self.confDict['targetINT2Search']:
            bkp.analyzeMetaclusters(metaclustersPass1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, self.identDbPath)
            if self.confDict['analyseFiltered']:
                bkp.analyzeMetaclusters(metaclustersFailed1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, self.identDbPath)

        # If only ME are analysed, pick a ME db that is given to the pipeline as an argument
        elif not 'VIRUS' in self.confDict['targetINT2Search'] and 'ME' in self.confDict['targetINT2Search']:
            bkp.analyzeMetaclusters(metaclustersPass1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, self.confDict['MEDB'])
            if self.confDict['analyseFiltered']:
                bkp.analyzeMetaclusters(metaclustersFailed1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, self.confDict['MEDB'])

        step = 'BKP-ANALYSIS'
        msg = 'Analysed BKP and added clipping supporting reads. PID: ' + str(os.getpid())
        log.step(step, msg)

        # NOTE SR: BY THIS WAY WE ARE ALSO ASSEGING MUT.ORIGIN TO THE METACLUSTER.
        # NOTE SR: All failed metaclusters are going to only one file.

        ## 6. Filter metaclusters before adding clipping reads ##

        filters2Apply = {}
        # Filters to apply to non-identified clusters
        filters2Apply['GENERIC'] = self.confDict['filtersAfClip']
        # Filters to apply to VIRUS clusters
        filters2Apply['VIRUS'] = self.confDict['filtersAfClip']
        # TODO SR: Add here proper filters for MEs clusters
        filters2Apply['ME']  = self.confDict['filtersAfClip']

        metaclustersList, metaclustersFailed2 = filters.filter_metaclusters_SR(metaclustersPass1, filters2Apply, self.confDict, self.bam)

        # Remove those metaclusters that didnt pass the filteres if the option for not printting them is selected
        if not self.confDict['printFiltered']:
            del metaclustersFailed2
            metaclustersFailedList = []
        else:
            metaclustersFailedList = metaclustersFailed1 + metaclustersFailed2

        step = 'FILTER METACLUSTERS AFTER ADDING CLIPPING READS'
        msg = 'Filter out metaclusters. Nb PASS metaclusters: '+ str(len(metaclustersList)) +'. Nb NOT PASS: ' + str(len(metaclustersFailedList)) + '. PID: ' + str(os.getpid())
        log.step(step, msg)

        ## 7. Extract VCF fileds ##

        # Get from metaclusters those fields that should be printed in VCF output
        metaclustersFields = output.VCFMetaclustersFields(metaclustersList)

        # Get from metaclusters those fields that should be printed in VCF output
        if not self.confDict['printFiltered']:
            metaclustersFailedFields = []
        else:
            metaclustersFailedFields = output.VCFMetaclustersFields(metaclustersFailedList)

        ### 8. Do cleanup
        unix.rm([binDir])

        # TODO SR: Not sure if this dels are needed, because function is finishing and I think they get deleted automatically.
        '''
        del metaclustersList
        del dictMetaclustersSecondFilterFailedTempList
        del metaclustersSVTypeFailedList
        del metaclustersFailedList
        del dictMetaclustersSecondFilter
        del dictMetaclustersSecondFilterFailedTemp
        '''

        return metaclustersFields, metaclustersFailedFields

class SV_caller_sureselect(SV_caller):
    '''
    Structural variation (SV) caller for Illumina sureselect data targetering source element´s downstream regions
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

    def call(self):
        '''
        Search for structural variants (SV) genome wide or in a set of target genomic regions
        '''

        ### 1. Create bed file containing source element´s transduced intervals 
        tdDir = self.outDir + '/TRANSDUCED/'
        unix.mkdir(tdDir)
        sourceBed = self.refDir + '/srcElements.bed'
        transducedPath = databases.create_transduced_bed(sourceBed, 10000, tdDir)
        
        ### 2. Infer read size
        self.confDict['readSize'] = self.infer_readSize()
                
        ### 3. Define genomic bins to search for SV (will correspond to transduced areas)
        bins = bamtools.binning(transducedPath, None, None, None)

        ## Organize bins into a dictionary
        self.confDict['rangesDict'] = gRanges.rangeList2dict(bins)

        ### 4. Associate to each bin the src identifier
        BED = formats.BED()
        BED.read(transducedPath, 'List', None)   

        for index, coords in enumerate(bins):
            coords.append(BED.lines[index].optional['cytobandId'])
        
        unix.rm([tdDir])
        
        ### 5. Search for SV clusters in each bin 
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        clusterPerSrc = pool.starmap(self.make_clusters_bin, bins)
        pool.close()
        pool.join()
        
        # Convert into dictionary
        clusterPerSrcDict = {srcId:clusters for srcId,clusters in clusterPerSrc}

        ### 6. Write calls to file
        ## 6.1 Transduction counts per source element
        output.write_tdCounts_sureselect(clusterPerSrcDict, self.outDir)

        ## 6.2 Transduction calls
        output.write_tdCalls_sureselect(clusterPerSrcDict, self.outDir)
    
    def infer_readSize(self):
        '''
        Infer read size from bam file
        '''
        
        # take the first 500 reads of the bam file
        command = 'samtools view ' + self.bam + '| awk \'{print length($10)}\' | head -500 | tr \'\n\' \' \''
        result = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        
        # if command fails, exit
        if result.returncode != 0:
            step = 'infer_readSize'
            msg = 'readSize inference failed' 
            log.step(step, msg)
            sys.exit(1)
        
        # save the result in a list of integers
        readSizes_str = result.stdout.decode('utf-8').split(" ")
        readSizes_str.remove("")
        readSizes_int = [int(i) for i in readSizes_str]
        
        # calculate the mode
        readSize = statistics.mode(readSizes_int)
        
        return(readSize)
        

    def make_clusters_bin(self, ref, beg, end, srcId):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''
        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        ## 1. Search for discordant and clipped read events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:
            eventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        step = 'COLLECT'
        SV_types = sorted(eventsDict.keys())
        counts = [str(len(eventsDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)

        ## If search for supplementary alignments selected:
        if self.confDict['blatClip']:
            
            ## 2. Search for supplementary alignments by realigning the clipped sequences
            step = 'SEARCH4SUPPL'
            msg = 'Search for supplementary alignments by realigning the clipped sequences'
            log.step(step, msg)

            ## Create output directory 
            supplDir = self.outDir + '/SUPPLEMENTARY/' + srcId
            unix.mkdir(supplDir)
            
            ## Left-clippings
            events.search4supplementary(eventsDict['LEFT-CLIPPING'], self.reference, srcId, supplDir)
                
            ## Rigth-clippings
            events.search4supplementary(eventsDict['RIGHT-CLIPPING'], self.reference, srcId, supplDir)

            ## Remove output directory
            unix.rm([supplDir])

        ## 3. Discordant and clipping clustering ##
        ## 3.1 Organize discordant and clipping events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize discordant and clipping events into genomic bins prior clustering'
        log.step(step, msg)
        
        ## Create bin database with discordants
        discordantsDict = {}
        discordantsDict['DISCORDANT'] = eventsDict['DISCORDANT']

        binSizes = [500, 1000, 10000, 100000, 1000000]  
        discordantsBinDb = structures.create_bin_database_interval(ref, beg, end, discordantsDict, binSizes)

        ## Create bin database with clippings 
        clippingsDict = {}
        clippingsDict['LEFT-CLIPPING'] = eventsDict['LEFT-CLIPPING']
        clippingsDict['RIGHT-CLIPPING'] = eventsDict['RIGHT-CLIPPING']

        binSizes = [self.confDict['maxBkpDist'], 100, 500, 1000, 10000, 100000, 1000000]
        clippingsBinDb = structures.create_bin_database_interval(ref, beg, end, clippingsDict, binSizes)

        ## 3.2 Group discordant and clipping events into clusters ##
        step = 'CLUSTERING'
        msg = 'Group discordant and clipping events into clusters'
        log.step(step, msg)
 
        ## Discordant clustering
        discordantClustersBinDb = clusters.create_clusters(discordantsBinDb, self.confDict)

        ## Clipping clustering
        clippingClustersBinDb = clusters.create_clusters(clippingsBinDb, self.confDict)

        ## 3.3 Group discordant read pairs based on mate position ##
        step = 'GROUP-BY-MATE'
        msg = 'Group discordant read pairs based on mate position'
        log.step(step, msg)

        ## Make groups
        discordants = clusters.cluster_by_matePos(discordantClustersBinDb.collect(['DISCORDANT']), self.refLengths, self.confDict['minClusterSize'])

        ## 3.4 Group clipping events based on suppl alignment position ##
        step = 'GROUP-BY-SUPPL'
        msg = 'Group clipping events based on suppl alignment position'
        log.step(step, msg)

        ## Left clipping
        leftClippingClusters = clusters.cluster_by_supplPos(clippingClustersBinDb.collect(['LEFT-CLIPPING']), self.refLengths, self.confDict['minClusterSize'], 'LEFT-CLIPPING')

        ## Right clipping
        rightClippingClusters = clusters.cluster_by_supplPos(clippingClustersBinDb.collect(['RIGHT-CLIPPING']), self.refLengths, self.confDict['minClusterSize'], 'RIGHT-CLIPPING')

        ## 4. Cluster filtering ##
        ## 4.1 Discordant cluster filtering ##
        step = 'FILTER-DISCORDANT'
        msg = 'Discordant cluster filtering'
        log.step(step, msg)

        filters2Apply = ['MIN-NBREADS', 'MATE-REF', 'MATE-SRC', 'MATE-MAPQ', 'GERMLINE', 'UNESPECIFIC', 'READ-DUP', 'CLUSTER-RANGE']
        filteredDiscordants = filters.filter_discordants(discordants, filters2Apply, self.bam, self.normalBam, self.confDict)

        ## 4.2 Clipping cluster filtering ##
        step = 'FILTER-CLIPPING'
        msg = 'Clipping cluster filtering'
        log.step(step, msg)
        
        filters2Apply = ['MIN-NBREADS', 'SUPPL-REF', 'SUPPL-SRC', 'SUPPL-MAPQ', 'GERMLINE', 'READ-DUP', 'CLUSTER-RANGE']
        filteredLeftClippings = filters.filter_clippings(leftClippingClusters, filters2Apply, self.confDict)
        filteredRightClippings = filters.filter_clippings(rightClippingClusters, filters2Apply, self.confDict)
            
        ## 5. Create metaclusters ##
        step = 'META-CLUSTERING'
        msg = 'Group discordant mates and suplementary clusters into metaclusters'
        log.step(step, msg)
        metaclusters = clusters.metacluster_mate_suppl(filteredDiscordants, filteredLeftClippings, filteredRightClippings, self.confDict['minReads'], self.refLengths)
        
        ## 6. Determine metaclusters precise coordinates ##
        bkp.bkp_retroTest(metaclusters, self.bam, self.confDict['readSize'])
        
        ## 7. Determine metaclusters identity ##           
        retrotransposons.identity_metaclusters_retrotest(metaclusters, self.bam, self.outDir)
        
        return [srcId, metaclusters]
