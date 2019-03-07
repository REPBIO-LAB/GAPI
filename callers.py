'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import multiprocessing as mp

# Internal
import log
import unix
import bamtools
import formats
import structures
import clustering
import variants
import filters
import output

## FUNCTIONS ##

## CLASSES ##
class SVcaller_nano():
    '''
    Structural variation (SV) caller from single molecule sequencing data
    '''
    def __init__(self, mode, bam, normalBam, reference, dbDir, confDict, outDir):

        self.mode = mode
        self.bam = bam
        self.normalBam = normalBam
        self.reference = reference
        self.dbDir = dbDir
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
            BED = formats.BED()
            BED.read(self.confDict['targetBins'])
            bins = [ (line.ref, line.beg, line.end) for line in BED.lines]

        ### 2. Search for SV clusters in each bin ##
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        INS_clusters, DEL_clusters, left_CLIPPING_clusters, right_CLIPPING_clusters = zip(*pool.map(self.SV_clusters_bin, bins))
        pool.close()
        pool.join()

        # Create dictionary containing clusters 
        clusters = {}
        clusters['INS-CLUSTER'] = INS_clusters
        clusters['DEL-CLUSTER'] = DEL_clusters
        clusters['LEFT-CLIPPING-CLUSTER'] = left_CLIPPING_clusters 
        clusters['RIGHT-CLIPPING-CLUSTER']= right_CLIPPING_clusters
        
        # Write output
        output.writeClusters(clusters, self.outDir)

    def SV_clusters_bin(self, window):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''

        ## 0. Set bin id and create bin directory
        ref, beg, end = window
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        binDir = self.outDir + '/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            INS_events, DEL_events, CLIPPING_left_events, CLIPPING_right_events, FASTQ = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:

            ## Search for SV events in the tumour
            INS_events_T, DEL_events_T, CLIPPING_left_events_T, CLIPPING_right_events_T, FASTQ_T = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, 'TUMOUR')

            ## Search for SV events in the normal
            INS_events_N, DEL_events_N, CLIPPING_left_events_N, CLIPPING_right_events_N, FASTQ_N = bamtools.collectSV(ref, beg, end, self.normalBam, self.confDict, 'NORMAL')

            ## Join tumour and normal lists
            INS_events = INS_events_T + INS_events_N
            DEL_events = DEL_events_T + DEL_events_N
            CLIPPING_left_events = CLIPPING_left_events_T + CLIPPING_left_events_N
            CLIPPING_right_events = CLIPPING_right_events_T + CLIPPING_right_events_N

            # Cleanup
            del INS_events_T, DEL_events_T, CLIPPING_left_events_T, CLIPPING_right_events_T  
            del INS_events_N, DEL_events_N, CLIPPING_left_events_N, CLIPPING_right_events_N

            ## Merge tumour and normal FASTQ     
            FASTQ = formats.FASTQ()
            FASTQ.fastqDict = {**FASTQ_T.fastqDict, **FASTQ_N.fastqDict} 
            
            # Cleanup
            del FASTQ_T, FASTQ_N  

        step = 'COLLECT'
        msg = 'Number of SV events in bin (INS, DEL, CLIPPING_left, CLIPPING_right): ' +  "\t".join([binId, str(len(INS_events)), str(len(DEL_events)), str(len(CLIPPING_left_events)), str(len(CLIPPING_right_events))])
        log.step(step, msg)

        ## 2. Organize all the SV events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior clustering'
        log.step(step, msg)

        ## 2.1 Insertions
        binSizes = [50]
        data = [(INS_events, 'INS')]
        INS_bins = structures.createBinDb(ref, beg, end, data, binSizes)

        # Cleanup
        del INS_events, data  

        ## 2.2 Deletions
        binSizes = [100, 1000, 10000, 100000, 1000000]
        data = [(DEL_events, 'DEL')]
        DEL_bins = structures.createBinDb(ref, beg, end, data, binSizes)

        # Cleanup
        del DEL_events, data  

        ## 2.3 Left-clippings
        binSizes = [50]
        data = [(CLIPPING_left_events, 'LEFT-CLIPPING')]
        left_CLIPPING_bins = structures.createBinDb(ref, beg, end, data, binSizes)

        # Cleanup
        del CLIPPING_left_events, data  

        ## 2.4 Right-clippings
        binSizes = [50]
        data = [(CLIPPING_right_events, 'RIGHT-CLIPPING')]
        right_CLIPPING_bins = structures.createBinDb(ref, beg, end, data, binSizes)

        # Cleanup
        del CLIPPING_right_events, data  

        ## 3. Group events into SV clusters ##
        ## 3.1 Cluster insertions
        INS_clusters = clustering.clusterByDist1D(INS_bins, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'INS')

        # Cleanup
        del INS_bins

        ## 3.2 Cluster deletions
        DEL_clusters = clustering.clusterByRcplOverlap(DEL_bins, self.confDict['minPercRcplOverlap'], self.confDict['minRootClusterSize'], 'DEL')

        # Cleanup
        del DEL_bins

        ## 3.3 Cluster clippings
        left_CLIPPING_clusters = clustering.clusterByDist1D(left_CLIPPING_bins, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'LEFT-CLIPPING')
        right_CLIPPING_clusters = clustering.clusterByDist1D(right_CLIPPING_bins, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'RIGHT-CLIPPING')

        step = 'CLUSTERING'
        msg = 'Number of clusters (INS, DEL, CLIPPING_left, CLIPPING_right): ' +  "\t".join([binId, str(INS_clusters.nbEvents()[0]), str(DEL_clusters.nbEvents()[0]), str(left_CLIPPING_clusters.nbEvents()[0]), str(right_CLIPPING_clusters.nbEvents()[0])])
        log.step(step, msg)

        # Cleanup
        del left_CLIPPING_bins, right_CLIPPING_bins

        ## 4. Polish SV clusters ##
        step = 'POLISH'
        msg = 'Polish SV clusters'
        log.step(step, msg)

        ## 4.1 Polish insertions
        variants.polishClusters(INS_clusters, 'INS-CLUSTER')

        ## 4.2 Polish deletions
        variants.polishClusters(DEL_clusters, 'DEL-CLUSTER')

        ## 5. Filter SV clusters ##
        step = 'FILTER'
        msg = 'Filter SV clusters'
        log.step(step, msg)

        ## 5.1 Filter insertions
        filters.filterClusters(INS_clusters, 'INS-CLUSTER', self.confDict)

        ## 5.2 Filter deletions
        filters.filterClusters(DEL_clusters, 'DEL-CLUSTER', self.confDict)

        ## 5.3 Filter left-clippings
        filters.filterClusters(left_CLIPPING_clusters, 'LEFT-CLIPPING-CLUSTER', self.confDict)

        ## 5.4 Filter right-clippings
        filters.filterClusters(right_CLIPPING_clusters, 'RIGHT-CLIPPING-CLUSTER', self.confDict)

        ## 6. Make consensus sequence for SV clusters ##
        step = 'CONSENSUS'
        msg = 'Make consensus sequence for SV clusters '
        log.step(step, msg)

        ## 6.1 Consensus for insertions
        variants.consensusClusters(INS_clusters, 'INS-CLUSTER', FASTQ, self.reference, self.confDict, binDir)

        ## 6.2 Consensus for deletions
        variants.consensusClusters(DEL_clusters, 'DEL-CLUSTER', FASTQ, self.reference, self.confDict, binDir)

        ## 7. Determine what has been inserted (insertion type) ##
        step = 'INS-TYPE'
        msg = 'Determine what has been inserted (insertion type)'
        log.step(step, msg)

        variants.insTypeClusters(INS_clusters, self.dbDir, self.confDict, binDir)

        ## Do Cleanup once bin was processed
        #unix.rm([binDir])

        return INS_clusters, DEL_clusters, left_CLIPPING_clusters, right_CLIPPING_clusters
