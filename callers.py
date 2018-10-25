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
import structures
import variants

## FUNCTIONS ##
def clusterCLIPPING(CLIPPING_list, maxBkpDist):
    '''
    Cluster CLIPPING events based on breakpoint position coordinates

    Input:
        1. CLIPPING_list: list of CLIPPING objects
        2. maxDist: maximum distance between two breakpoints to include them into the same cluster

    Output:
        1. clusterDict: positions dictionary containing the clusters
    ''' 
    step = 'CLUSTER-CLIPPING'
    msg = 'Input (CLIPPING_list, maxBkpDist): ' +  "\t".join([str(len(CLIPPING_list)), str(maxBkpDist)])    
    log.step(step, msg)

    ## 1. Organize CLIPPING events by their breakpoint position ##
    posDict = structures.buildPosDict(CLIPPING_list, maxBkpDist)

    ## 2. Cluster CLIPPING events ##
    # Initialize list of already processed windows
    processedWindows = []
    clusterList = []

    # For each genomic window
    for windowIndex, events in posDict.items():
    
        # Skip already processed windows
        if windowIndex in processedWindows:
            continue

        processedWindows.append(windowIndex)

        ## Create cluster containing CLIPPINGS on the region
        cluster = variants.CLIPPING_cluster(events)
        clusterList.append(cluster)

        ## 2.1 Go backward from current window (*).  
        #       <---2--- <---1---
        # |---------|--------|----*----|---------
        # Go one window backward in each iteration. Extend the cluster 
        # if last clipping breakpoint in previous window within 
        # maximum cluster distance or end iteration, otherwise 
        backwardIndex = windowIndex - 1

        while True:
        
            # A) There are CLIPPINGS in the window 
            if backwardIndex in posDict:
                
                # Compute the distance between the right most CLIPPING in the new window
                # and the cluster begin 
                CLIPPINGS = posDict[backwardIndex]
                lastCLIPPING = CLIPPINGS[-1]
                bkpDist = cluster.pos - lastCLIPPING.pos

                # a) Last CLIPPING within maximum distance 
                if bkpDist <= maxBkpDist:

                    processedWindows.append(backwardIndex)

                    # Add CLIPPINGS within window to the cluster 
                    cluster.add(CLIPPINGS, 'left')
         
                # b) CLIPPING outside 
                else:
                    break

            # B) No CLIPPINGS in the window. Stop iterating
            else:
                break 
            
            backwardIndex -= 1
        
        ## 2.2 Go forward from current window (*).  
        #                ---1---> ---2--->
        # |---------|----*----|--------|---------
        # Go one window forward in each iteration. Extend the cluster 
        # if first clipping breakpoint in next window within 
        # maximum cluster distance or end iteration, otherwise 
        forwardIndex = windowIndex + 1

        while True:
        
            # A) There are CLIPPINGS in the window 
            if forwardIndex in posDict:
                
                # Compute the distance between the left most CLIPPING in the new window
                # and the cluster end 
                CLIPPINGS = posDict[forwardIndex]
                firstCLIPPING = CLIPPINGS[0]
                bkpDist = firstCLIPPING.pos - cluster.end 

                # a) Last CLIPPING within maximum distance 
                if bkpDist <= maxBkpDist:

                    processedWindows.append(forwardIndex)

                    # Add CLIPPINGS within window to the cluster 
                    cluster.add(CLIPPINGS, 'right')
         
                # b) CLIPPING outside 
                else:
                    break

            # B) No CLIPPINGS in the window. Stop iterating
            else:
                break 
            
            forwardIndex += 1        

    ## 2. Organize CLIPPING clusters into a dictionary ##
    clustersDict = structures.buildPosDict(clusterList, 1000)

    return clustersDict

## Pending, create these functions:
#clusterINS
#clusterDEL


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
                INS_list, DEL_list, CLIPPING_left_list, CLIPPING_right_list = bamtools.collectSV(ref, beg, end, self.bam, self.confDict['targetSV'], None)

            # b) Paired sample mode (tumour & matched normal)
            else:

                ## Search for SV events in the tumour
                INS_list_T, DEL_list_T, CLIPPING_left_list_T, CLIPPING_right_list_T = bamtools.collectSV(ref, beg, end, self.bam, self.confDict['targetSV'], 'TUMOUR')

                ## Search for SV events in the normal
                INS_list_N, DEL_list_N, CLIPPING_left_list_N, CLIPPING_right_list_N = bamtools.collectSV(ref, beg, end, self.normalBam, self.confDict['targetSV'], 'NORMAL')

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
            CLIPPING_left_clustersDict = clusterCLIPPING(CLIPPING_left_list, self.confDict['maxBkpDist'])
            CLIPPING_right_clustersDict = clusterCLIPPING(CLIPPING_right_list, self.confDict['maxBkpDist'])

            ## 2.2 Cluster insertions 

            ## 2.3 Cluster deletions




            
