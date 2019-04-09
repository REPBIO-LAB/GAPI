'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import multiprocessing as mp

# Internal
import log
import unix
import databases
import bamtools
import structures
import clustering
import clusters
import filters
import output

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
        self.retrotransposonDb = None
        self.retrotransposonDbIndex = None


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
        ### 1. Create and index reference databases prior SV calling ##
        dbDir = self.outDir + '/databases'
        unix.mkdir(dbDir)

        self.retrotransposonDb, self.retrotransposonDbIndex = databases.buildRetrotransposonDb(self.refDir, self.confDict['transductionSearch'], dbDir)

        ### 2. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 3. Search for SV clusters in each bin ##
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        INS_clusters, DEL_clusters, left_CLIPPING_clusters, right_CLIPPING_clusters = zip(*pool.map(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        ### 4. Report clusters into output file
        # Create dictionary containing clusters 
        clusters = {}
        clusters['INS-CLUSTER'] = INS_clusters
        clusters['DEL-CLUSTER'] = DEL_clusters
        clusters['LEFT-CLIPPING-CLUSTER'] = left_CLIPPING_clusters 
        clusters['RIGHT-CLIPPING-CLUSTER']= right_CLIPPING_clusters
    
        # Write clusters
        output.writeClusters(clusters, self.outDir)

        ### 5. Do cleanup
        unix.rm([dbDir])

    def make_clusters_bin(self, window):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''

        ## 0. Set bin id and create bin directory ##
        ref, beg, end = window
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        binDir = self.outDir + '/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
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

        ## 2. Organize all the SV events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior metaclustering'
        log.step(step, msg)

        ## Define bin database sizes 
        binSizes = [self.confDict['maxEventDist'], 1000, 10000, 100000, 1000000]

        ## Create bins
        eventsBinDb = structures.create_bin_database(ref, beg, end, eventsDict, binSizes)

        ## 3. Group events into SV metaclusters ##
        metaclustersBinDb = clusters.create_metaclusters(eventsBinDb, self.confDict)
        
        step = 'META-CLUSTERING'
        msg = 'Number of created metaclusters: ' + str(metaclustersBinDb.nbEvents()[0])
        log.step(step, msg)

        
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
        ### 1. Create and index reference databases prior SV calling ##
        #dbDir = self.outDir + '/databases'
        #unix.mkdir(dbDir)

        #self.retrotransposonDb, self.retrotransposonDbIndex = databases.buildRetrotransposonDb(self.refDir, self.confDict['transductionSearch'], dbDir)

        ### 2. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 3. Search for SV clusters in each bin ##
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        left_DISCORDANT_clusters, right_DISCORDANT_clusters = zip(*pool.map(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        ### 4. Report clusters into output file
        # Create dictionary containing clusters 
        clusters = {}
        clusters['LEFT-DISCORDANT-CLUSTER'] = left_DISCORDANT_clusters 
        clusters['RIGHT-DISCORDANT-CLUSTER']= right_DISCORDANT_clusters
    
        # Write clusters
        output.writeClusters(clusters, self.outDir)

        ### 5. Do cleanup
        #unix.rm([dbDir])

    def make_clusters_bin(self, window):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''
        ## 0. Set bin id and create bin directory ##
        ref, beg, end = window
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        binDir = self.outDir + '/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            discordantEventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:
            discordantEventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        step = 'COLLECT'
        SV_types = sorted(discordantEventsDict.keys())

        counts = [str(len(discordantEventsDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
    
        log.step(step, msg)

        ## 2. Organize discordant events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior metaclustering'
        log.step(step, msg)

        ## Define bin database sizes 
        ## Big window sizes are needed for SR (see comments)
        binSizes = [self.confDict['maxEventDist'], 10000, 100000, 1000000]

        ## Create bins
        eventsBinDb = structures.create_bin_database(ref, beg, end, discordantEventsDict, binSizes)

        ## 3. Group discordant events into metaclusters ##
        metaclustersBinDb = clusters.create_metaclusters(eventsBinDb, self.confDict)

        step = 'META-CLUSTERING'
        msg = 'Number of created metaclusters: ' + str(metaclustersBinDb.nbEvents()[0])
        log.step(step, msg)