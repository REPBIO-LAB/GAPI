#!/usr/bin/env python
#coding: utf-8


def worker(bam, normalBam, windowsList, mode, outDir):
    '''
    '''
    threadId = threading.currentThread().getName()
    msg = threadId + ' launched'
    log.info(msg)
    time.sleep(random.uniform(0, 1))
    
    print("THREAD: ", bam, normalBam, mode, threadId)

    ## Process each window separately
    for window in windowsList:
    
        ref, beg, end = window

        print("WINDOW: ", ref, beg, end, bam, normalBam, window, mode, outDir)
        #targetSV = ['INS', 'DEL', 'CLIPPING']
        targetSV = ['INS', 'CLIPPING']
        bamtools.collectSV(ref, beg, end, bam, targetSV, 'TUMOUR')

        
        # a) Single sample mode
        #if mode == "SINGLE":
        #    callMEIWindow(chrom, beg, end, bam, outDir)
        
        # b) Paired sample mode (tumour & matched normal)
        #else:
        #    callMEIWindow_paired(chrom, beg, end, bam, normalBam, outDir)

#### MAIN ####

## Import modules ##
import argparse
import sys
import os
import time
import pysam
import threading
import random
import bamtools
import log

# Global variables:
global debugBool ## debug logging mode. Boolean.

# Environmental variables:

## Get user's input ##
parser = argparse.ArgumentParser(description= "Identify retrotransposon and viral insertions from Nanopore/Pacbio whole genome sequencing data. Two running modes: 1) SINGLE: individual sample; 2) PAIRED: tumour and matched normal samples")
parser.add_argument('bam', help='Input BAM file. Will correspond to the tumour sample in the PAIRED mode')
parser.add_argument('--normal-bam', default="NA", dest='normalBam', help='Matched normal bam file')
parser.add_argument('-t', '--threads', default=1, dest='threads', type=int, help='Number of threads. Default: 1' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
bam = args.bam
normalBam = args.normalBam
threads = args.threads
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Determine running mode:
mode = "SINGLE" if normalBam == "NA" else "PAIRED"


## Display configuration to standard output ##
print()
print("***** ", scriptName, " configuration *****")
print("mode: ", mode)
print("bam: ", bam)
print("normalBam: ", normalBam)
print("threads: ", threads)
print("outDir: ", outDir)
print()
print("***** Executing ", scriptName, ".... *****")
print()


## Start ## 

###  1. Split the reference genome into a set of genomic bins
##############################################################
windowSize = 1000000

targetRefList = ['22']
#targetRefList = list(range(1, 23, 1))
#targetRefList = None

windowsList = bamtools.makeGenomicBins(bam, windowSize, targetRefList)
  
### 2. Search for 
#######################################################
# Split windows list into X evenly sized chunks. X = number of threads
chunkSize = int(round(len(windowsList) / float(threads)))
chunksList = [windowsList[i:i + chunkSize] for i in range(0, len(windowsList), chunkSize)]

threads = list()
counter = 1

## Assing one thread per chunk
for chunk in chunksList:
    msg = "chunk" + str(counter) + ": " + str(len(chunk)) + " windows to process"
    log.info(msg)

    threadName = "THREAD-" + str(counter)
    thread = threading.Thread(target=worker, args=(bam, normalBam, chunk, mode, outDir), name=threadName)
    threads.append(thread)
    counter += 1

# Launch threads
[t.start() for t in threads]

# Wait till threads have finished with MEI calling
[t.join() for t in threads]

print("***** Finished! *****")
print()

