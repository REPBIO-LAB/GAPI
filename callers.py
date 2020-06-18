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
            clusters.INS_type_metaclusters(allMetaclusters['INS'], self.reference, self.annotations, 1, outDir)

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
            allMetaclusters['BND'] = clusters.search4bridges_metaclusters_parallel(allMetaclusters['BND'], 10000, 80, self.confDict['minReads'], 25, self.annotations, self.refDir, self.confDict['processes'], outDir)

            ### Search for BND junctions
            allJunctions = clusters.search4junctions_metaclusters(allMetaclusters['BND'], self.refLengths, self.confDict['processes'], self.confDict['minReads'], 25)
            
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
        clusters.lighten_up_metaclusters(metaclustersSVType)
        
        # Do cleanup
        unix.rm([outDir, binDir])

        ## Print time taken to process bin
        end = time.time()
        time_taken = end - start
        msg = 'SV calling in bin: ' + binId + ' finished in ' + str(time_taken)
        log.info(msg)

        return metaclustersSVType


class SV_caller_short(SV_caller):
    '''
    Structural variation (SV) caller for Illumina short read sequencing data
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

    def call(self):
        '''
        Search for structural variants (SV) genome wide or in a set of target genomic regions
        '''
        ### 1. Create, index and load reference databases prior SV calling ##
        dbDir = self.outDir + '/DATABASES/'
        unix.mkdir(dbDir)

        ## 1.1 Load annotated retrotransposons into a bin database
        ## Read bed
        rtAnnotBed = self.refDir + '/repeats_repeatMasker.bed'
        #rtAnnotBed = self.refDir + '/repeats_repeatMasker.L1.bed'

        rtAnnot = formats.BED()
        rtAnnot.read(rtAnnotBed, 'nestedDict', None) 

        ## Create bin database
        self.repeatsBinDb = structures.create_bin_database_parallel(self.refLengths, rtAnnot.lines, 1)
        
        ## 1.2 Create and index viral database
        #self.viralDb, self.viralDbIndex = databases.buildVirusDb(self.refDir, dbDir)
        
        ## 1.3 Create transduced regions database
        # a) Create database if transduction search enabled
        if self.confDict['transductionSearch']:

            ## Create bed file
            sourceBed = self.refDir + '/srcElements.bed'
            transducedPath = databases.create_transduced_bed(sourceBed, 15000, dbDir)

            ## Read bed
            transducedBed = formats.BED()
            transducedBed.read(transducedPath, 'nestedDict', None)

            ## Create bin database
            self.transducedBinDb = structures.create_bin_database_parallel(self.refLengths, transducedBed.lines, 1)

        # b) Skip database creation
        else:
            self.transducedBinDb = None

        ### 2. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 3. Search for SV clusters in each bin ##
        # Genomic bins will be distributed into X processes
        # TODO: mirar que pasa cuando tienes 2 dictionarios
        pool = mp.Pool(processes=self.confDict['processes'])
        discordantClusters = pool.starmap(self.make_clusters_bin, bins)
        #self.make_clusters_bin('22',16050000,16150000)   #TODO DEBUG ONLY
        #self.make_clusters_bin('9', 18720898, 18722291)   #TODO DEBUG ONLY
    
        pool.close()
        pool.join()
    
        # Report SV calls into output files
        output.write_DISCORDANT(discordantClusters, self.outDir)

        ### 5. Do cleanup
        unix.rm([dbDir])

    def make_clusters_bin(self, ref, beg, end):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''
        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        binDir = self.outDir + '/CLUSTER/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            discordantDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:
            discordantDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        step = 'COLLECT'
        SV_types = sorted(discordantDict.keys())
        counts = [str(len(discordantDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)

        if counts == []:
            unix.rm([binDir])
            return None
                
        ## 2. Discordant read pair identity ##
        ## Determine identity
        discordantsIdentity = events.determine_discordant_identity(discordantDict['DISCORDANT'], self.repeatsBinDb, self.transducedBinDb)

        step = 'IDENTITY'
        SV_types = sorted(discordantsIdentity.keys())
        counts = [str(len(discordantsIdentity[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events per identity in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)

        if counts == []:
            unix.rm([binDir])
            return None
                
        ## 3. Organize discordant read pairs into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize discordant read pairs into genomic bins prior clustering'
        log.step(step, msg)

        ## Define bin database sizes 
        ## Note: bigger window sizes are needed for SR (see comments, ask Eva where are the comments?)
        binSizes = [1000, 10000, 100000, 1000000]

        ## Create bins
        discordantsBinDb = structures.create_bin_database_interval(ref, beg, end, discordantsIdentity, binSizes)
        
        ## 4. Group discordant read pairs into clusters based on their mate identity ##
        buffer = 100
        discordantClustersDict = clusters.create_discordantClusters(discordantsBinDb, self.confDict['minClusterSize'], buffer)
    
        step = 'DISCORDANT-CLUSTERING'
        SV_types = sorted(discordantClustersDict.keys())
        counts = [str(len(discordantClustersDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of created discordant clusters in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)

        if counts == []:
            unix.rm([binDir])
            return None

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

        ### Do cleanup
        unix.rm([binDir])

        ##
        return discordantClustersDict
        
        '''
        Eva will further polish next steps!

        ## 6. Filter discordant metaclusters ##

        filters.filterClusters(discordantClustersBinDb, ['DISCORDANT'], self.confDict, self.bam)

        ## Remove those clusters that fail in one or more filters
        newDiscordantClustersDict = filters.applyFilters(discordantClustersBinDb)

        step = 'DISCORDANT-FILTERING'
        SV_types = sorted(discordantClustersDict.keys())

        counts = [str(len(discordantClustersDict[SV_type])) for SV_type in SV_types]
        msg = '[DISCORDANT-FILTERING] Number of created discordant clusters after filtering in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.subHeader(msg)

        if counts == []:
            unix.rm([binDir])
            return None

        # vuelvo a hacer la bindb que contiene ya solo los clusters que pasaron los filtros
        discordantClustersBinDb = structures.create_bin_database_interval(ref, beg, end, newDiscordantClustersDict, binSizes)

        ## 7. Making reciprocal clusters ##
        # TODO: AJUSTAR ESTOS PARAMETROS!!! (PASARLOS SI ESO COMO OPCION EN LOS ARGUMENTOS)
        reciprocalEventsDict = clustering.reciprocal(discordantClustersBinDb, 1, 1, 300)

        ## 8. Put reciprocal and independent clusters into bins ##

        reciprocalEventsBinDb = structures.create_bin_database_interval(ref, beg, end, reciprocalEventsDict, binSizes)
        '''
      
        '''
        # HASTA AQUI YO CREO QUE ESTA TODO BIEN!!!! A PARTIR DE AQUI HAY QUE REPASAR!!

        ## 8. Create metaclusters from reciprocal and independent clusters ##

        metaclustersBinDb = clusters.create_metaclusters(reciprocalEventsBinDb, self.confDict)

        step = 'META-CLUSTERING'
        msg = '[META-CLUSTERING] Number of created metaclusters: ' + str(metaclustersBinDb.nbEvents()[0])
        log.subHeader(msg)
        '''
        '''
        LO ULTIMO QUE HICE FUE RETORNAR EVENTS DE LA RECIPROCAL EN VEZ DE CLUSTERS, Y FUNCIONA, PERO HAY EN ALGUN MOMENTO QUE SE MEZCLAN LOS DE DISTINTO TIPO AL HACER LA RECIPROCAL, ASI QUE TENGO QUE REPASARLO!
        
        {'DISCORDANT-Hepatitis': [<events.DISCORDANT object at 0x7f73aeb5c208>, <events.DISCORDANT object at 0x7f73aeb5c780>, <events.DISCORDANT object at 0x7f73aeb55780>, <events.DISCORDANT object at 0x7f73aeb55a58>], 'DISCORDANT-UNVERIFIED:': [<events.DISCORDANT object at 0x7f73aeb5c208>, <events.DISCORDANT object at 0x7f73aeb5c780>, <events.DISCORDANT object at 0x7f73aeb55780>, <events.DISCORDANT object at 0x7f73aeb55a58>], 'DISCORDANT-HBV': [<events.DISCORDANT object at 0x7f73aeb5c2b0>, <events.DISCORDANT object at 0x7f73aeb5c2b0>, <events.DISCORDANT object at 0x7f73aeb55780>, <events.DISCORDANT object at 0x7f73aeb55a58>]}
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6390> 2 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6390> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6390> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> 5 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:7:2108:9935:21524/1 2 105457573 DISCORDANT UNVERIFIED: MINUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6400> 2 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6400> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6400> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> 5 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:7:2108:9935:21524/1 2 105457573 DISCORDANT UNVERIFIED: MINUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6470> 2 2 105457243 105457720
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6470> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6470> HWI-ST672:120:D0CF5ACXX:8:2104:17314:11147/1 2 105457619 DISCORDANT HBV MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> 5 2 105457243 105457720
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:120:D0CF5ACXX:8:2104:17314:11147/1 2 105457619 DISCORDANT HBV MINUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:120:D0CF5ACXX:8:2104:17314:11147/1 2 105457619 DISCORDANT HBV MINUS

        '''

        '''
        dictMetaclustersLEFT = bkp.analizeBkp(metaclustersBinDb, self.viralDb, self.reference, 'LEFT', binDir)
        dictMetaclustersRIGHT = bkp.analizeBkp(metaclustersBinDb, self.viralDb, self.reference, 'RIGHT', binDir)

        
        print ('LEFT')
        print (dictMetaclustersLEFT)

        print ('RIGHT')
        print (dictMetaclustersRIGHT)
        '''
        '''
        dictMetaclusters = bkp.analyzeMetaclusters(metaclustersBinDb, self.confDict, self.bam, self.normalBam, self.mode, self.viralDb, self.viralDbIndex, binDir)

        ### Do cleanup
        unix.rm([binDir])

        return dictMetaclusters
        '''

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

        ### Do cleanup
        supplDir = self.outDir + '/SUPPLEMENTARY/'
        unix.rm([supplDir])    
    
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

        return [srcId, metaclusters]
        