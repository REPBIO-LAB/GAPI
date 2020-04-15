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
import gRanges
import clustering
import sequences
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

        ### 1. Create SV clusters 
        msg = '1. Create SV clusters'
        log.header(msg)
        allMetaclusters = self.make_clusters()
        
        ### 2. Annotate SV clusters intervals  
        msg = '2. Annotate SV clusters intervals'
        log.header(msg)

        # Create output directory
        annotDir = self.outDir + '/ANNOT/'
        unix.mkdir(annotDir)

        # Reference lengths, needed for repeats annotation
        refLengths = bamtools.get_ref_lengths(self.bam)
        
        # Define annotation steps
        steps = ['REPEAT']

        if self.confDict['annovarDir'] is not None:
            steps.append('GENE')

        # For each cluster type
        for SV_type in allMetaclusters:
            
            metaclusters = allMetaclusters[SV_type]
            annotation.annotate(metaclusters, steps, refLengths, self.refDir, self.confDict['annovarDir'], self.confDict['processes'], annotDir)

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
            clusters.INS_type_metaclusters(allMetaclusters['INS'], self.reference, refLengths, self.refDir, self.confDict['transductionSearch'], 1, outDir)

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
            outDir = self.outDir + '/BRIDGES/'
            unix.mkdir(outDir)   

            allMetaclusters['BND'] = clusters.search4bridges_metaclusters_parallel(allMetaclusters['BND'], 10000, 80, self.confDict['minSupportingReads'], 25, refLengths, self.refDir, self.confDict['processes'], outDir)

            ### Search for BND junctions
            allJunctions = clusters.search4junctions_metaclusters(allMetaclusters['BND'], refLengths, self.confDict['processes'], self.confDict['minSupportingReads'], 25)
            
            # Remove output directory
            unix.rm([outDir])

        ### 6. Apply second round of filtering 
        msg = '6. Apply second round of filtering'
        log.header(msg)
        filters2Apply = ['PERC-RESOLVED']
        metaclustersPass, metaclustersFailed = filters.filter_metaclusters(allMetaclusters, filters2Apply, self.confDict)
                
        ### 7. Report SV calls into output files
        msg = '7. Report SV calls into output files'
        log.header(msg)
        
        ##  7.1 Report INS
        if 'INS' in metaclustersPass:
            outFileName = 'INS_MEIGA.PASS.tsv'
            output.write_INS(metaclustersPass['INS'], outFileName, self.outDir)

        if 'INS' in metaclustersFailed:
            outFileName = 'INS_MEIGA.FAILED.2.tsv'
            output.write_INS(metaclustersFailed['INS'], outFileName, self.outDir)

        ## 7.2 Report BND junctions
        if allJunctions:
            outFileName = 'BND_junctions_MEIGA.PASS.tsv'
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
        metaclustersPassList, metaclustersFailedList = zip(*pool.starmap(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        ### 3. Collapse metaclusters in a single dict and report metaclusters that failed filtering
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
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None, True)

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

        ## Compute reference lengths
        self.refLengths = bamtools.get_ref_lengths(self.bam)


    def minimap2_index(self):
        '''
        Return path to minimap2 index file
        '''
        index = os.path.splitext(self.reference)[0] + '.mmi' 

        return index

    def call(self):
        '''
        Search for integrations genome wide or in a set of target genomic regions
        '''
        # NOTE SR: if 'ME' in self.confDict['targetINT2Search'] annotRepeats is not and option
        if 'ME' in self.confDict['targetINT2Search']:
            annotDir = self.outDir + '/ANNOT/'
            self.annotations = annotation.load_annotations(['REPEATS', 'TRANSDUCTIONS'], self.refLengths, self.refDir, self.confDict['processes'], annotDir)
        else:
            self.annotations = {}
            self.annotations['REPEATS'], self.annotations['TRANSDUCTIONS'] = None, None

        ### 1. Create integration clusters 
        msg = '1. Create integration clusters. PID: ' + str(os.getpid())
        log.header(msg)      
        '''
        EXPLANATION
        Tuple of lists of lists -> metaclustersListofLists ([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], {}, {}, [], [], [], [],
        [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>',
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
            self.annotations = annotation.load_annotations(['REPEATS', 'TRANSDUCTIONS'], self.refLengths, self.refDir, self.confDict['processes'], annotDir)

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
            annotation.repeats_annotation(metaclustersList, self.annotations['REPEATS'], buffer)
            
            ## 6. Perform gene-based annotation with ANNOVAR of discordant read pair clusters ##
            # Do gene-based annotation step if enabled
            
            if self.confDict['annovarDir'] is not None:

                step = 'ANNOTATE'
                msg = 'Perform gene-based annotation with ANNOVAR of discordant read-pair clusters'
                log.step(step, msg)

                ## Annotate
                annotation.gene_annotation(metaclustersList, self.confDict['annovarDir'], annotDir)

            # Remove annotation directory
            unix.rm([annotDir])

        # TODO SR: Think if is worth it to make viral db inside MEIGA (headers, etc)
        # TODO SR: Index viral db
        ## 1.2 Create and index viral database
        #self.viralDb, self.viralDbIndex = databases.buildVirusDb(self.refDir, dbDir)
        
  
        # Report integrations calls into output files
        #output.write_DISCORDANT(discordantClusters, self.outDir)
        #metaclustersListofLists = list(allMetaclusters.values())
        outFileNameTSV = 'metaclusters.PASS.tsv'
        outFileName = 'metaclusters.PASS'
        #output.writeMetaclusters(metaclustersList, outFileNameTSV, self.outDir)
        output.INS2VCF_SR(metaclustersList, self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], self.confDict['VCFInfoFields'], self.confDict['VCFREF'], outFileName, self.outDir)


    def make_clusters(self):

        ### If viruses option is selected, collect read name and sequence of discordant low quality reads from all bam refs ##
        if 'VIRUS' in self.confDict['targetINT2Search']:
            # TEMP SR: DESILENCE
            '''
            # Make genomic bins
            bins = bamtools.makeGenomicBins(self.bam, self.confDict['binSize'], None)
            '''
            # Collect read name and sequence of discordant low quality reads from all bam refs
            l = mp.Lock()
            # TEMP SR: DESILENCE
            '''
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
                        
			# Filter by complexity (with komplexity)
            allFastas = sequences.komplexityFilter(self.confDict['komplexityThreshold'], 'allFastas_all.fasta', 'allFastas.fasta', self.outDir)
            
            # Align with bwa allFastas vs viralDb and filter resulting bam
            # TODO SR: bwa allFastas vs viralDb: check if bwa -T parameter does something that we need
            BAM = alignment.alignment_bwa_filtered(self.confDict['viralDb'], self.confDict['viralBamParcialMatch'], self.confDict['processes'], allFastas, 'viralAligment', self.outDir)

            # Index bam
            bamtools.samtools_index_bam(BAM, self.outDir)

            # TEMP SR: Remove allfastas
            #unix.rm([allFastas])
            '''
            #TEMP
            BAM = self.outDir + '/' + 'viralAligment' + '.bam'
            # Read bwa result and store in a dictionary
            self.viralSeqs = bamtools.BAM2FastaDict(BAM)
            
        ### 1. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 2. Search for SV clusters in each bin ##
        # Create output directory
        unix.mkdir(self.outDir + '/CLUSTER/')

        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        metaclustersPassListofLists, metaclustersFailedListofLists = zip(*pool.starmap(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        del self.viralSeqs
        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        '''
        EXPLANATION
        Tuple of lists of lists -> metaclustersFailedListofLists ([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], {}, {}, [], [], [], [],
        [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>',
        '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}]], [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}]], [], [], [])
        '''
        if metaclustersFailedListofLists:
            #metaclustersFailedListofLists = list(metaclustersFailed.values())
            outFileNameTSV = 'metaclusters.FAILED.tsv'
            outFileName = 'metaclusters.FAILED'
            '''
            EXPLANATION
            Lists of lists -> metaclustersFailedList [[CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, POS, '.', '<INS>', '.', 'PASS', {INFO_Dict}], [CHROM, 56715108, '.', '<INS>', '.', 'PASS', {INFO_Dict}]]
            '''
            # Flat metaclustersFailedList
            metaclustersFailedList = list(itertools.chain(*metaclustersFailedListofLists))
            #output.writeMetaclusters(metaclustersFailedList, outFileName, self.outDir)
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
            discordants = bamtools.collectDISCORDANT(ref, beg, end, self.bam, self.confDict, None, False)

        # b) Paired sample mode (tumour & matched normal)
        else:
            discordants = bamtools.collectDISCORDANT_paired(ref, beg, end, self.bam, self.normalBam, self.confDict, False)


        counts = str(len(discordants))
        
        step = 'COLLECT'
        msg = 'Number of DISCORDANT events in bin ' + binId + ': ' + counts + '. PID: ' + str(os.getpid())
        log.step(step, msg)
        
        ## 2. Discordant read pair identity ##
        ## Determine identity
        # NOTE SR: I think adapting for PAIRED mode is not needed
        #if self.mode == "SINGLE":
            # TODO SR: If we want to analyse RT, we should call determine_discordant_identity in another way, depending if we are analysing RT, virus or both.


        discordantsIdentity = events.determine_discordant_identity(discordants, None, None,self.bam, None, binDir, self.confDict['targetINT2Search'], self.viralSeqs)

        #discordantsIdentity = events.determine_discordant_identity(discordants, None, None,self.bam, None, binDir, self.confDict['targetINT2Search'])
        #else:
            #discordantsIdentity = events.determine_discordant_identity(discordantDict['DISCORDANT'], self.annotations['REPEATS'], self.annotations['TRANSDUCTIONS'],self.bam, None, binDir, self.confDict['viralDb'])
            #discordantsIdentity = events.determine_discordant_identity(discordantDict['DISCORDANT'], None, None,self.bam, None, binDir, self.confDict['viralDb'])
        
        del discordants

        step = 'IDENTITY'
        SV_types = sorted(discordantsIdentity.keys())
        counts = [str(len(discordantsIdentity[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events per identity in bin ' + binId + '(' + ','.join(SV_types) + '): ' + ','.join(counts) + '. PID: ' + str(os.getpid())
        log.step(step, msg)

        ## 3. Organize discordant read pairs into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize discordant read pairs into genomic bins prior clustering'
        log.step(step, msg)

        ## Define bin database sizes 
        minBinSize = min([self.confDict['maxInsDist'], self.confDict['maxBkpDist']])
        binSizes = [minBinSize, 1000, 10000, 100000, 1000000]

        ## Create bins
        discordantsBinDb = structures.create_bin_database_interval(ref, beg, end, discordantsIdentity, binSizes)

        del discordantsIdentity

        ## 4. Group discordant read pairs into clusters based on their mate identity ##
        buffer = 100
        discordantClustersDict = clusters.create_discordantClusters(discordantsBinDb, self.confDict['minClusterSize'], buffer)

        del discordantsBinDb
    
        step = 'DISCORDANT-CLUSTERING'
        SV_types = sorted(discordantClustersDict.keys())
        counts = [str(len(discordantClustersDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of created discordant clusters in bin ' + binId + '(' + ','.join(SV_types) + '): ' + ','.join(counts) + '. PID: ' + str(os.getpid())
        log.step(step, msg)

        # Return if no DISCODANT clusters found.
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
        annotation.repeats_annotation(allDiscordantClusters, self.annotations['REPEATS'], buffer)
        
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

        ### Do cleanup
        #unix.rm([binDir])

        ## 5. Filter discordant clusters ##
        filters2Apply = {}
        # Filters to apply to non-identified clusters
        filters2Apply['GENERIC'] = ['MAX-NBREADS', 'AREAMAPQ', 'AREASMS', 'IDENTITY']
        # Filters to apply to VIRUS clusters
        # NOTE SR: dont look at 'MIN-NBREADS' at cluster level
        filters2Apply['VIRUS'] = ['MAX-NBREADS', 'AREAMAPQ', 'AREASMS', 'IDENTITY']
        # TODO SR: Add here proper filters for MEs clusters
        filters2Apply['ME']  = ['MAX-NBREADS', 'AREAMAPQ', 'AREASMS', 'IDENTITY']
        discordantClustersDict, discordantClustersDictFailed = filters.filter_clusters(discordantClustersDict, filters2Apply, self.confDict, self.bam)

        step = 'FILTER-CLUSTERS'
        msg = 'Filter out clusters'  + '. PID: ' + str(os.getpid())
        log.step(step, msg)

        ## 6. Organize discordant clusters in bin database structure ##
        discordantClustersBinDb = structures.create_bin_database_interval(ref, beg, end, discordantClustersDict, binSizes)
        discordantClustersFailedBinDb = structures.create_bin_database_interval(ref, beg, end, discordantClustersDictFailed, binSizes)

        del discordantClustersDict
        del discordantClustersDictFailed

        ## 7. Make reciprocal clusters ##
        # TODO SR: Adjust parameters of clustering.reciprocal. Maybe it is worth it make them running arguments
        reciprocalClustersDict = clustering.reciprocal(discordantClustersBinDb, 1, 1, 300)
        reciprocalClustersFailedDict = clustering.reciprocal(discordantClustersFailedBinDb, 1, 1, 300)

        step = 'RECIPROCAL-CLUSTERING'
        msg = 'Performing reciprocal clustering. PID: ' + str(os.getpid())
        log.step(step, msg)

        del discordantClustersBinDb
        del discordantClustersFailedBinDb
     
        ## 8. Organize reciprocal and independent discordant clusters in bin database structure ##
        reciprocalClustersBinDb = structures.create_bin_database_interval(ref, beg, end, reciprocalClustersDict, binSizes)
        reciprocalClustersFailedBinDb = structures.create_bin_database_interval(ref, beg, end, reciprocalClustersFailedDict, binSizes)

        del reciprocalClustersDict
        del reciprocalClustersFailedDict

        ## 9. Get all identities ##
        identities = set([iden.split('-')[2] for iden in reciprocalClustersBinDb.eventTypes])
        identitiesFailed = set([iden.split('-')[2] for iden in reciprocalClustersFailedBinDb.eventTypes])

        metaclusters=[]
        metaclustersFailed=[]
        ## Create metaclusters from reciprocal and independent discordant clusters taking into account identities (eventType = those eventTypes (MINUS, PLUS and RECIPROCAL) that correspond to identity).
        ## By this way, only clusters with same identity are metaclustered toguether.
        # TODO SR: Review here, when calling create_discordant_metaclusters, maybe there are some parameters to adjust.
        for identity in identities:
            currentEventTypes = [eventType for eventType in reciprocalClustersBinDb.eventTypes if (identity in eventType)]
            metaclusters.extend(clusters.create_discordant_metaclusters(reciprocalClustersBinDb, currentEventTypes))

        for identityFailed in identitiesFailed:
            currentEventTypesFailed = [eventType for eventType in reciprocalClustersFailedBinDb.eventTypes if (identityFailed in eventType)]
            metaclustersFailed.extend(clusters.create_discordant_metaclusters(reciprocalClustersFailedBinDb, currentEventTypesFailed))

        step = 'META-CLUSTERING'
        msg = 'Number of created metaclusters: PASS -> ' + str(len(metaclusters)) + ' FILTERED -> ' + str(len(metaclustersFailed)) + '. PID: ' + str(os.getpid())
        log.step(step, msg)

        ## 10. Analyse metaclusters features and add supporting clipping reads ##
        # TODO SR: Now adding clipping step is better than before. Even though maybe it is not neccessary, it would be great to choose in a better way which clipping we should add to the metacluster.
        bkp.analyzeMetaclusters(metaclusters, self.confDict, self.bam, self.normalBam, self.mode, binDir)
        bkp.analyzeMetaclusters(metaclustersFailed, self.confDict, self.bam, self.normalBam, self.mode, binDir)

        step = 'BKP-ANALYSIS'
        msg = 'Analysing BKP and adding clipping supporting reads. PID: ' + str(os.getpid())
        log.step(step, msg)

        metaclustersSVTypeBfSecondFilter = {}
        metaclustersSVTypeFailed = {}
      
        metaclustersSVTypeBfSecondFilter['DISCORDANT'] = metaclusters
        metaclustersSVTypeFailed['DISCORDANT'] = metaclustersFailed

        del metaclusters
        del metaclustersFailed

        # NOTE SR: BY THIS WAY WE ARE ALSO ASSEGING MUT.ORIGIN TO THE METACLUSTER.
        # NOTE SR: All failed metaclusters are going to only one file.
        # TODO SR: I think this could return a list instead of a dictionary.

        filters2Apply = {}
        # Filters to apply to non-identified clusters
        filters2Apply['GENERIC'] = ['MIN-NBREADS', 'MAX-NBREADS']
        # Filters to apply to VIRUS clusters
        filters2Apply['VIRUS'] = ['MIN-NBREADS', 'MAX-NBREADS']
        # TODO SR: Add here proper filters for MEs clusters
        filters2Apply['ME']  = ['MIN-NBREADS', 'MAX-NBREADS']
        dictMetaclustersSecondFilter, dictMetaclustersSecondFilterFailedTemp = filters.filter_metaclusters(metaclustersSVTypeBfSecondFilter, filters2Apply, self.confDict)

        step = 'FILTER METACLUSTERS'
        msg = 'Filter out metaclusters' 
        log.step(step, msg)

        # Flat list
        #metaclustersList = list(itertools.chain(*dictMetaclustersSecondFilter.values()))
        metaclustersList = [item for sublist in dictMetaclustersSecondFilter.values() for item in sublist]

        # Get from metaclusters those fields that should be printed in VCF output
        metaclustersFields = output.VCFMetaclustersFields(metaclustersList)

        # Flat lists
        dictMetaclustersSecondFilterFailedTempList = list(itertools.chain(*dictMetaclustersSecondFilterFailedTemp.values()))
        metaclustersSVTypeFailedList = list(itertools.chain(*metaclustersSVTypeFailed.values()))
        # Merge list of failed clusters and list of failed metaclusters
        metaclustersFailedList = dictMetaclustersSecondFilterFailedTempList + metaclustersSVTypeFailedList
        # Get from metaclusters those fields that should be printed in VCF output
        metaclustersFailedFields = output.VCFMetaclustersFields(metaclustersFailedList)

        # TODO SR: Not sure if this dels are needed, because function is finishing and I think they get deleted automatically.
        del metaclustersList
        del dictMetaclustersSecondFilterFailedTempList
        del metaclustersSVTypeFailedList
        del metaclustersFailedList
        del dictMetaclustersSecondFilter
        del dictMetaclustersSecondFilterFailedTemp

        ### Do cleanup
        unix.rm([binDir])

        '''
        ## 11. Lighten up metaclusters  ##
        ## Metaclusters passing all the filters        
        clusters.lighten_up_metaclusters(metaclustersSVType)

        ## Filtered metaclusters
        clusters.lighten_up_metaclusters(metaclustersSVTypeFailed)
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
                
        ### 2. Define genomic bins to search for SV (will correspond to transduced areas)
        bins = bamtools.binning(transducedPath, None, None, None)

        ## Organize bins into a dictionary
        self.rangesDict = gRanges.rangeList2dict(bins)

        ### 3. Associate to each bin the src identifier
        BED = formats.BED()
        BED.read(transducedPath, 'List', None)   

        for index, coords in enumerate(bins):
            coords.append(BED.lines[index].optional['cytobandId'])
        
        unix.rm([tdDir])

        ### 4. Search for SV clusters in each bin 
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        clusterPerSrc = pool.starmap(self.make_clusters_bin, bins)
        pool.close()
        pool.join()

        # Convert into dictionary
        clusterPerSrcDict = {srcId:clusters for srcId,clusters in clusterPerSrc}

        ### 5. Write calls to file
        ## 5.1 Transduction counts per source element
        output.write_tdCounts_surelect(clusterPerSrcDict, self.outDir)

        ## 5.2 Transduction calls
        output.write_tdCalls_surelect(clusterPerSrcDict, self.outDir)

    def make_clusters_bin(self, ref, beg, end, srcId):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''

        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None, True)

        # b) Paired sample mode (tumour & matched normal)
        else:
            eventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        step = 'COLLECT'
        SV_types = sorted(eventsDict.keys())
        counts = [str(len(eventsDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)
        
        ## 2. Organize discordant read pairs into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize discordant read pairs into genomic bins prior clustering'
        log.step(step, msg)

        ## Define bin database sizes 
        ## Note: bigger window sizes are needed for SR (see comments, ask Eva where are the comments?)
        binSizes = [500, 1000, 10000, 100000, 1000000]

        ## Create bins
        discordantsBinDb = structures.create_bin_database_interval(ref, beg, end, eventsDict, binSizes)

        ## 3. Group discordant read pairs into clusters based on reciprocal overlap ##
        step = 'CLUSTERING'
        msg = 'Group discordant read pairs into clusters based on reciprocal overlap'
        log.step(step, msg)

        buffer = 150
        discordantClustersDict = clusters.create_discordantClusters(discordantsBinDb, self.confDict['minClusterSize'], buffer)

        ## 4. Do an extra clustering step based on mate position ##
        step = 'GROUP-BY-MATE'
        msg = 'Group discordant read pairs based on mate position'
        log.step(step, msg)

        ## Compute reference lengths
        refLengths = bamtools.get_ref_lengths(self.bam)

        ## Make groups
        discordantClustersDict = clusters.extra_clustering_by_matePos(discordantClustersDict, refLengths, self.confDict['minClusterSize'])

        ## 5. Filter out those clusters over NOT target reference ##
        step = 'FILTER-REF'
        msg = 'Filter out those clusters over NOT target reference'
        log.step(step, msg)
        filteredDiscordants = filters.filter_discordant_mate_ref(discordantClustersDict['DISCORDANT'], self.confDict['targetRefs'])

        ## 6. Filter out those clusters whose mates aligns over any source element downstream region ##
        step = 'FILTER-DOWSTREAM'
        msg = 'Filter out those clusters whose mates aligns over any source element downstream region'
        log.step(step, msg)
        filteredDiscordants = filters.filter_discordant_mate_position(filteredDiscordants, self.rangesDict, 10000)        
        
        ## 7. Filter out those clusters based on average MAPQ for mate alignments ##
        step = 'FILTER-MATE-MAPQ'
        msg = 'Filter out those clusters based on average MAPQ for mate alignments'
        log.step(step, msg)
        filteredDiscordants = filters.filter_discordant_mate_MAPQ(filteredDiscordants, 20, self.bam)

        return [srcId, filteredDiscordants]  

