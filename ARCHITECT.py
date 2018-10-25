#!/usr/bin/env python
#coding: utf-8



#### MAIN ####

## Import modules ##
# External
import argparse
import sys
import os

# Internal
import callers
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

## 1. Create configuration dictionary
confDict = {}
windowSize = 1000000
maxBkpDist = 100
targetRefs = ['22']
#targetRefs = list(range(1, 23, 1))
#targetRefs = None
targetSV = ['INS', 'CLIPPING']

confDict['threads'] = threads
confDict['windowSize'] = windowSize
confDict['maxBkpDist'] = maxBkpDist
confDict['targetRefs'] = targetRefs
confDict['targetSV'] = targetSV

## 2. Launch structural variation (SV) caller
callerObj = callers.SVcaller_nano(mode, bam, normalBam, confDict, outDir)
callerObj.callSV()

print("***** Finished! *****")
print()

