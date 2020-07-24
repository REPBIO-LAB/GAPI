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
from VIGA_SR import bamtools_VIGASR
from VIGA_SR import virusSR


#from modules.callers import SV_caller

## CLASSES ##
# TODO: When everything is organised put this back in callers (as it is commented), and import it from there. I cant do it now because it gets in travel with MEIGA.py importing callers and formats importing call_VIGASR.
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

        # If ME analysis is done, this is loaded at the beggining
        # DESILENCE
        if self.confDict['annotRepeats']:
            annotDir = self.outDir + '/ANNOT/'
            self.annotations = annotation.load_annotations(['REPEATS', 'TRANSDUCTIONS'], self.refLengths, self.refDir, self.confDict['germlineMEI'], self.confDict['processes'], annotDir)

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
        self.viralSeqs, self.identDbPath = virusSR.find_virus_discordants(self.bam, self.normalBam, self.confDict['viralDb'], self.confDict['komplexityThreshold'], self.confDict['minTotalMatchVirus'], self.confDict['minParcialMatchVirus'], self.confDict['maxMatchCheckMAPQVirus'], self.confDict['minMAPQVirus'], self.confDict['maxBasePercVirus'], self.confDict['minLccVirus'], self.confDict['processes'], self.outDir)
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
        bamtools_VIGASR.collectDiscodantsLowMAPQSeq(ref, binBeg, binEnd, self.bam, self.confDict['discordantMatesMaxMAPQ'], self.confDict['discordantMatesCheckUnmapped'], self.confDict['discordantMatesSupplementary'], self.confDict['discordantMatesMaxBasePerc'], self.confDict['discordantMatesMinLcc'], self.outDir)

    def callCollectSeqNormal(self, ref, binBeg, binEnd):
        '''
        Wrapper to call function for collecting read name and sequence of discordant low quality reads from all NORMAL bam refs
        '''
        bamtools_VIGASR.collectDiscodantsLowMAPQSeq(ref, binBeg, binEnd, self.normalBam, self.confDict['discordantMatesMaxMAPQ'], self.confDict['discordantMatesCheckUnmapped'], self.confDict['discordantMatesSupplementary'], self.confDict['discordantMatesMaxBasePerc'], self.confDict['discordantMatesMinLcc'], self.outDir)

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
        discordantsIdentity = events.determine_discordant_identity(discordants, None, None, self.bam, None, binDir, ['VIRUS'], self.viralSeqs)

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
        # Pick viral identities db which is created internally.
        bkp.analyzeMetaclusters(metaclustersPass1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, self.identDbPath)
        if self.confDict['analyseFiltered']:
            bkp.analyzeMetaclusters(metaclustersFailed1, self.confDict, self.bam, self.normalBam, self.mode, binDir, binId, self.identDbPath)

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