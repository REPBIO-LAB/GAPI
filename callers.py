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
        windowsList = bamtools.makeGenomicBins(self.bam, self.confDict['windowSize'], self.confDict['targetRefs'])

        ## 2. Split windows list into X evenly sized chunks. X = number of threads
        chunkSize = int(round(len(windowsList) / float(self.confDict['threads'])))
        chunksList = [windowsList[i:i + chunkSize] for i in range(0, len(windowsList), chunkSize)]

        ## 3. Assign one thread per chunk
        threads = list()
        counter = 1

        for chunk in chunksList:
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

    def callSV_windows(self, windowsList):
        '''
        Search for structural variants (SV) in a set of genomic windows of interest
        '''
        threadId = threading.currentThread().getName()
        msg = threadId + ' launched'
        log.info(msg)
        time.sleep(random.uniform(0, 1))
    
        # print("THREAD: ", self.bam, self.normalBam, self.mode, threadId)

        ## For each window
        for window in windowsList:
    
            ref, beg, end = window

            msg = "SV calling in window: " + "_".join([str(ref), str(beg), str(end)]) 
            log.subHeader(msg)

            ## 1. Search for SV candidate events in the bam file/s ##
            # a) Single sample mode
            if self.mode == "SINGLE":
                INS_list, DEL_list, CLIPPING_left_list, CLIPPING_right_list = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

            # b) Paired sample mode (tumour & matched normal)
            else:

                ## Search for SV events in the tumour
                INS_list_T, DEL_list_T, CLIPPING_left_list_T, CLIPPING_right_list_T = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, 'TUMOUR')

                ## Search for SV events in the normal
                INS_list_N, DEL_list_N, CLIPPING_left_list_N, CLIPPING_right_list_N = bamtools.collectSV(ref, beg, end, self.normalBam, self.confDict, 'NORMAL')

                ## Join tumour and normal lists
                INS_list = INS_list_T + INS_list_N
                DEL_list = DEL_list_T + DEL_list_N
                CLIPPING_left_list = CLIPPING_left_list_T + CLIPPING_left_list_N
                CLIPPING_right_list = CLIPPING_right_list_T + CLIPPING_right_list_N

            step = 'COLLECT'
            msg = 'Number of SV events (INS, DEL, CLIPPING_left, CLIPPING_right): ' +  "\t".join([str(len(INS_list)), str(len(DEL_list)), str(len(CLIPPING_left_list)), str(len(CLIPPING_right_list))])    
            log.step(step, msg)

            ## 2. Group events into SV clusters ##
            ## 2.1 Cluster CLIPPINGS  
            CLIPPING_left_clusters, nbClustersLeft = clustering.clusterByPos(CLIPPING_left_list, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'])
            CLIPPING_right_clusters, nbClustersRight = clustering.clusterByPos(CLIPPING_right_list, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'])

            step = 'CLUSTER-CLIPPING'
            msg = 'Number of CLIPPING clusters (left clipping, right clipping): ' +  "\t".join([str(nbClustersLeft), str(nbClustersRight)])    
            log.step(step, msg)
            
            ## 2.2 Cluster insertions 
            INS_clusters, nbINS = clustering.clusterByPos(INS_list, self.confDict['maxBkpDist'], self.confDict['minRootClusterSize'])

            step = 'CLUSTER-INS'
            msg = 'Number of INS clusters: ' +  str(nbINS)    
            log.step(step, msg)

            for windowIndex, clustersWindow in INS_clusters.items():
                print('WINDOW', windowIndex, [(cluster.ref, cluster.beg, cluster.end, len(cluster.events)) for cluster in clustersWindow])
            
            ## 2.3 Cluster deletions
            # We should use a different clustering algorithm for deletions. Only rely on deletion begin makes no sense (POS)
            # I think a better strategy would be to rely on reciprocal overlap. I.E. those deletions with at least 80% with reciprocal 
            # overlap are clustered together. 

            ## 2. Group clusters into SV metaclusters ##






            
