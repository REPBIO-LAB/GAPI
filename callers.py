'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import sys
import multiprocessing as mp
import subprocess
import os
import pysam

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
        ### 1. Define genomic bins to search for SV ##
        msg = '1. Define genomic bins to search for SV'
        log.header(msg)

        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])
        
        ### 2. Search for SV clusters in each bin ##
        # Genomic bins will be distributed into X processes
        msg = '2. Search for SV clusters in each bin'
        log.header(msg)

        pool = mp.Pool(processes=self.confDict['processes'])
        metaclustersPassList, metaclustersFailedList = zip(*pool.starmap(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        ### 3. Collapse metaclusters in a single dict and report metaclusters that failed first round of filtering
        msg = '3. Collapse metaclusters in a single dict and report metaclusters that failed first round of filtering'
        log.header(msg)

        metaclustersPass_round1 = structures.merge_dictionaries(metaclustersPassList)
        metaclustersFailed_round1 = structures.merge_dictionaries(metaclustersFailedList)

        if 'INS' in metaclustersFailed_round1:
            outFileName = 'INS_MEIGA.FAILED_1.tsv'
            output.write_INS(metaclustersFailed_round1['INS'], outFileName, self.outDir)

        ## Remove metaclusters that failed filtering to release memory 
        del metaclustersPassList
        del metaclustersFailedList
        del metaclustersFailed_round1
        
        ## 4. Perform repeat-based annotation of INS target region
        msg = '4. Perform repeat-based annotation of INS target region'
        log.header(msg)

        if ('INS' in metaclustersPass_round1):

            ## 4.1. Load repeats database ##
            step = 'REPEAT-ANNOT'
            msg = '4.1. Load repeats database'
            log.step(step, msg)

            annotDir = self.outDir + '/ANNOT/'
            annotations2load = ['REPEATS']   

            refLengths = bamtools.get_ref_lengths(self.bam)
            annotations = annotation.load_annotations(annotations2load, refLengths, self.refDir, self.confDict['processes'], annotDir)
            repeatsAnnot = annotations['REPEATS']

            ## 4.2. Perform repeats annotation ##
            msg = '4.2. Perform repeats annotation'
            log.step(step, msg)

            annotation.repeats_annotation(metaclustersPass_round1['INS'], repeatsAnnot, 200)

        ### 5. Determine what type of sequence has been inserted for INS metaclusters
        msg = '5. Determine what type of sequence has been inserted for INS metaclusters'
        log.header(msg)

        if 'INS' in metaclustersPass_round1:

            ## 5.1. Load transduced regions and exons database ##
            step = 'INS-TYPE'
            msg = '5.1. Load transduced regions and exons database if appropiate'
            log.step(step, msg)

            annotDir = self.outDir + '/ANNOT/'
            annotations2load = []

            if self.confDict['transductionSearch']:    
                annotations2load.append('TRANSDUCTIONS')

            #if True: # at one point include flag for pseudogene search
                #annotations2load.append('EXONS')

            annotations = annotation.load_annotations(annotations2load, refLengths, self.refDir, self.confDict['processes'], annotDir)
            annotations['REPEATS'] = repeatsAnnot

            ## 5.2. Insertion type inference
            msg = '5.2. Insertion type inference'
            log.step(step, msg)

            outDir = self.outDir + '/INS_TYPE/'
            unix.mkdir(outDir)

            clusters.INS_type_metaclusters(metaclustersPass_round1['INS'], self.reference, annotations['REPEATS'], annotations['TRANSDUCTIONS'], annotations['EXONS'], self.confDict, outDir)
            
            # Cleanup
            unix.rm([outDir])

            # Remove annotation databases to release memory
            del annotations 
            del repeatsAnnot
            
        ### 6. Resolve structure for solo, partnered and orphan transductions
        msg = '6. Resolve structure for solo, partnered and orphan transductions'
        log.header(msg)
        
        if 'INS' in metaclustersPass_round1:
            consensus = self.refDir + '/consensusDb.fa'
            transduced = self.refDir + '/transducedDb.fa.masked'

            outDir = self.outDir + '/STRUCTURE/'
            unix.mkdir(outDir)

            ## Reduce the number of processes to the half to decrease RAM usage
            processes = int(self.confDict['processes']/2)

            metaclustersPass_round1['INS'] = clusters.structure_inference_parallel(metaclustersPass_round1['INS'], consensus, transduced, self.confDict['transductionSearch'], processes, outDir)

            # Cleanup
            unix.rm([outDir])

        ### 7. Apply second round of filtering after insertion type inference 
        msg = '7. Apply second round of filtering after insertion type inference'
        log.header(msg)
        filters2Apply = ['PERC-RESOLVED']
        metaclustersPass, metaclustersFailed_round2 = filters.filter_metaclusters(metaclustersPass_round1, filters2Apply, self.confDict)

        # Remove variables to release memory
        del metaclustersPass_round1
        
        ## 8. Perform gene-based annotation with ANNOVAR of INS target region
        msg = '8. Perform gene-based annotation with ANNOVAR of INS target region'
        log.header(msg)

        # Do gene-based annotation step if enabled
        if ('INS' in metaclustersPass) and (self.confDict['annovarDir'] is not None):

            step = 'GENE-ANNOT'
            msg = 'Perform gene-based annotation with ANNOVAR of INS events'
            log.step(step, msg)

            ## Annotate
            annotDir = self.outDir + '/ANNOT/' 
            annotation.gene_annotation(metaclustersPass['INS'], self.confDict['annovarDir'], annotDir)

            # Cleanup
            unix.rm([annotDir])
        
        ### 9. Report SV calls into output files
        msg = '9. Report SV calls into output files'
        log.header(msg)

        ##  9.1 Report INS
        if 'INS' in metaclustersPass:
            outFileName = 'INS_MEIGA.PASS.tsv'
            output.write_INS(metaclustersPass['INS'], outFileName, self.outDir)

        if 'INS' in metaclustersFailed_round2:
            outFileName = 'INS_MEIGA.FAILED_2.tsv'
            output.write_INS(metaclustersFailed_round2['INS'], outFileName, self.outDir)
        

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
        metaclustersSVType, metaclustersSVTypeFailed = filters.filter_metaclusters(metaclustersSVType, filters2Apply, self.confDict)

        ## 8. Generate consensus event for SV metaclusters ##
        step = 'CONSENSUS'
        msg = 'Generate consensus event for SV metaclusters' 
        log.step(step, msg)

        targetSV = ['INS']
        outDir = binDir + '/CONSENSUS/' 
        clusters.create_consensus(metaclustersSVType, self.confDict, self.reference, targetSV, outDir)       

        # Do cleanup
        unix.rm([outDir, binDir])

        return metaclustersSVType, metaclustersSVTypeFailed


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
        rtAnnotBed = self.refDir + '/retrotransposons_repeatMasker.bed'
        rtAnnot = formats.BED()
        rtAnnot.read(rtAnnotBed, 'nestedDict')

        ## Create bin database
        refLengths = bamtools.get_ref_lengths(self.bam)
        self.repeatsBinDb = structures.create_bin_database(refLengths, rtAnnot.lines, 1)
        
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
            transducedBed.read(transducedPath, 'nestedDict')

            ## Create bin database
            self.transducedBinDb = structures.create_bin_database(refLengths, transducedBed.lines, 1)

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

