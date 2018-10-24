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
    
        print("THREAD: ", self.bam, self.normalBam, self.mode, threadId)

        ## For each window
        for window in windowsList:
    
            ref, beg, end = window

            print("WINDOW: ", ref, beg, end, self.mode)

            ## 1. Search for SV in the bam file/s ##
            # a) Single sample mode
            if self.mode == "SINGLE":
                INS_list, DEL_list, CLIPPING_left_list, CLIPPING_right_list = bamtools.collectSV(ref, beg, end, self.bam, self.confDict['targetSV'], None)

            # b) Paired sample mode (tumour & matched normal)
            else:

                ## Search for SV in the tumour
                INS_list_T, DEL_list_T, CLIPPING_left_list_T, CLIPPING_right_list_T = bamtools.collectSV(ref, beg, end, self.bam, self.confDict['targetSV'], 'TUMOUR')

                ## Search for SV in the normal
                INS_list_N, DEL_list_N, CLIPPING_left_list_N, CLIPPING_right_list_N = bamtools.collectSV(ref, beg, end, self.normalBam, self.confDict['targetSV'], 'NORMAL')

                ## Join tumour and normal lists
                INS_list = INS_list_T + INS_list_N
                DEL_list = DEL_list_T + DEL_list_N
                CLIPPING_left_list = CLIPPING_left_list_T + CLIPPING_left_list_N
                CLIPPING_right_list = CLIPPING_right_list_T + CLIPPING_right_list_N

            

            ## 2. Organize events into a dictionary by genomic position ##
            event_list = INS_list + DEL_list + CLIPPING_left_list + CLIPPING_right_list

            print('TIOOO', self.mode, len(event_list), len(INS_list), len(DEL_list), len(CLIPPING_left_list), len(CLIPPING_right_list))

            bkpDictObj = structures.breakpointDict()

            bkpDictObj.buildBkpDict(event_list, self.confDict['clusterSize'])

