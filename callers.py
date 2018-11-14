'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import multiprocessing as mp

# Internal
import log
import bamtools
import formats
import structures
import clustering

## FUNCTIONS ##

## CLASSES ##
class SVcaller_nano():
    '''
    Structural variation (SV) caller from single molecule sequencing data
    '''    
    def __init__(self, mode, bam, normalBam, confDict, outDir):

        self.mode = mode
        self.bam = bam
        self.normalBam = normalBam
        self.confDict = confDict
        self.outDir = outDir

    def callSV(self):
        '''
        Search for structural variants (SV) genome wide or in a set of target genomic regions
        '''
        ### 1. Define genomic bins to search for SV ##
        # a) Create bins the novo
        if self.confDict['targetBins'] == None:

            ## Split the reference genome into a set of genomic bins
            bins = bamtools.makeGenomicBins(self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        # b) Read bins from bed file
        else:
            bed = formats.bed()
            bed.read(self.confDict['targetBins'])
            bins = [ (line.ref, line.beg, line.end) for line in bed.lines]

        ### 2. Call SV in each bin ##
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        pool.map(self.callSV_bin, bins)
        pool.close()
        pool.join()


    def callSV_bin(self, window):
        '''
        Search for structural variants (SV) in a genomic bin/window
        '''
    
        ref, beg, end = window

        msg = "SV calling in bin: " + "_".join([str(ref), str(beg), str(end)]) 
        log.subHeader(msg)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            INS_events, DEL_events, CLIPPING_left_events, CLIPPING_right_events = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:

            ## Search for SV events in the tumour
            INS_events_T, DEL_events_T, CLIPPING_left_events_T, CLIPPING_right_events_T = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, 'TUMOUR')

            ## Search for SV events in the normal
            INS_events_N, DEL_events_N, CLIPPING_left_events_N, CLIPPING_right_events_N = bamtools.collectSV(ref, beg, end, self.normalBam, self.confDict, 'NORMAL')

            ## Join tumour and normal lists
            INS_events = INS_events_T + INS_events_N
            DEL_events = DEL_events_T + DEL_events_N
            CLIPPING_left_events = CLIPPING_left_events_T + CLIPPING_left_events_N
            CLIPPING_right_events = CLIPPING_right_events_T + CLIPPING_right_events_N

        step = 'COLLECT'
        msg = 'Number of SV events (INS, DEL, CLIPPING_left, CLIPPING_right): ' +  "\t".join([str(len(INS_events)), str(len(DEL_events)), str(len(CLIPPING_left_events)), str(len(CLIPPING_right_events))])    
        log.step(step, msg)

        ## 2. Organize all the SV events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior clustering'
        log.step(step, msg)

        ## 2.1 Insertions
        binSizes = [50]
        data = [(INS_events, 'INS')]
        INS_bins = structures.createBinDb(data, binSizes)

        ## 2.2 Deletions
        # missing events: 22_24000000_25000000, 22_23000000_2400000. I guess due to size longer than maximum bin size (1Mb)
        binSizes = [100, 1000, 10000, 100000, 1000000]
        data = [(DEL_events, 'DEL')]
        DEL_bins = structures.createBinDb(data, binSizes)
            
        ## 2.3 Left-clippings
        binSizes = [50]
        data = [(CLIPPING_left_events, 'LEFT-CLIPPING')]
        left_CLIPPING_bins = structures.createBinDb(data, binSizes)

        ## 2.4 Right-clippings
        binSizes = [50]
        data = [(CLIPPING_right_events, 'RIGHT-CLIPPING')]
        right_CLIPPING_bins = structures.createBinDb(data, binSizes)

        ## 3. Group events into SV clusters ##     
        ## 3.1 Cluster insertions       
        INS_clusters = clustering.clusterByDist1D(INS_bins, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'INS')
        step = 'CLUSTER-INS'
        msg = 'Number of INS clusters: ' +  str(INS_clusters.nbEvents()[0])  
        log.step(step, msg)    

        ## 2.2 Cluster deletions            
        DEL_clusters = clustering.clusterByRcplOverlap(DEL_bins, self.confDict['minPercRcplOverlap'], self.confDict['minRootClusterSize'], 'DEL')
        step = 'CLUSTER-DEL'
        msg = 'Number of DEL clusters: ' +  str(DEL_clusters.nbEvents()[0])  
        log.step(step, msg)

        ## 2.3 Cluster clippings  
        left_CLIPPING_clusters = clustering.clusterByDist1D(left_CLIPPING_bins, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'LEFT-CLIPPING')
        right_CLIPPING_clusters = clustering.clusterByDist1D(right_CLIPPING_bins, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'RIGHT-CLIPPING')
        step = 'CLUSTER-CLIPPING'
        msg = 'Number of CLIPPING clusters (left clipping, right clipping): ' +  "\t".join([str(left_CLIPPING_clusters.nbEvents()[0]), str(right_CLIPPING_clusters.nbEvents()[0])])    
        log.step(step, msg)
