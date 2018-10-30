'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import time
import random
import threading

# Internal
import log
import bamtools
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
        Search for structural variants (SV) genome wide
        '''
        ## 1. Split the reference genome into a set of genomic bins
        windows = bamtools.makeGenomicBins(self.bam, self.confDict['windowSize'], self.confDict['targetRefs'])

        ## 2. Split windows list into X evenly sized chunks. X = number of threads
        chunkSize = int(round(len(windows) / float(self.confDict['threads'])))
        chunks = [windows[i:i + chunkSize] for i in range(0, len(windows), chunkSize)]

        ## 3. Assign one thread per chunk
        threads = list()
        counter = 1

        for chunk in chunks:
            msg = "chunk" + str(counter) + ": " + str(len(chunk)) + " windows to process"
            log.info(msg)

            threadName = "THREAD-" + str(counter)
            thread = threading.Thread(target=self.callSV_windows, args=([chunk]), name=threadName)
            threads.append(thread)
            counter += 1

        # Launch threads
        [t.start() for t in threads]

        # Wait till threads have finished with MEI calling
        [t.join() for t in threads]

    def callSV_windows(self, windows):
        '''
        Search for structural variants (SV) in a set of genomic windows of interest
        '''
        threadId = threading.currentThread().getName()
        msg = threadId + ' launched'
        log.info(msg)
        time.sleep(random.uniform(0, 1))
    
        # print("THREAD: ", self.bam, self.normalBam, self.mode, threadId)

        ## For each window
        for window in windows:
    
            ref, beg, end = window

            msg = "SV calling in window: " + "_".join([str(ref), str(beg), str(end)]) 
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

            ## 2. Group events into SV clusters ##          
            ## 2.1 Cluster insertions 
            INS_clusters = clustering.clusterByPos(INS_events, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'INS')
            
            step = 'CLUSTER-INS'
            msg = 'Number of INS clusters: ' +  str(INS_clusters.nbEvents())  
            log.step(step, msg)

            #for windowIndex, clustersWindow in INS_clusters.data.items():
            #    print('INS-WINDOW', windowIndex, [(cluster.ref, cluster.beg, cluster.end, len(cluster.events)) for cluster in clustersWindow])
    
            ## 2.2 Cluster deletions
            # We should use a different clustering algorithm for deletions. Only rely on deletion begin makes no sense (POS)
            # I think a better strategy would be to rely on reciprocal overlap. I.E. those deletions with at least 80% with reciprocal 
            # overlap are clustered together. 
            
            ## 2.3 Cluster CLIPPINGS  
            CLIPPING_left_clusters = clustering.clusterByPos(CLIPPING_left_events, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'CLIPPING')
            CLIPPING_right_clusters = clustering.clusterByPos(CLIPPING_right_events, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'], 'CLIPPING')

            step = 'CLUSTER-CLIPPING'
            msg = 'Number of CLIPPING clusters (left clipping, right clipping): ' +  "\t".join([str(CLIPPING_left_clusters.nbEvents()), str(CLIPPING_right_clusters.nbEvents())])    
            log.step(step, msg)

            '''
            for windowIndex, clustersWindow in CLIPPING_left_clusters.data.items():
                print('CLIPPING-LEFT-WINDOW', windowIndex, [(cluster.ref, cluster.beg, cluster.end, len(cluster.events)) for cluster in clustersWindow])
           
            for windowIndex, clustersWindow in CLIPPING_right_clusters.data.items():
                print('CLIPPING-RIGHT-WINDOW', windowIndex, [(cluster.ref, cluster.beg, cluster.end, len(cluster.events)) for cluster in clustersWindow])
            '''
            
            ## 3. Group clusters into SV metaclusters ##
            #metaTypes = ['INS+CLIPPING']
            #clustering.buildMetaClusters(INS_clusters, None, CLIPPING_left_clusters, CLIPPING_right_clusters, metaTypes, self.confDict)
            
            
