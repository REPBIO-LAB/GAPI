#!/usr/bin/env python
#coding: utf-8

## Import modules ##
# External
import argparse
import sys
import os

# Internal
import callers
import log
import bamtools

# Global variables:
global debugBool ## debug logging mode. Boolean.

# Environmental variables:

## Get user's input ##
parser = argparse.ArgumentParser(description= "Identify retrotransposon and viral insertions from Nanopore/Pacbio whole genome sequencing data. Two running modes: 1) SINGLE: individual sample; 2) PAIRED: tumour and matched normal samples")
parser.add_argument('bam', help='Input BAM file. Will correspond to the tumour sample in the PAIRED mode')
parser.add_argument('--normal-bam', default="NA", dest='normalBam', help='Matched normal bam file')
parser.add_argument('-t', '--threads', default=1, dest='threads', type=int, help='Number of threads. Default: 1' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

parser.add_argument('-wS', '--windowSize', default=1000000, dest='windowSize', type=int, help='Input bams will be analised in windows of this size. Default: 1000000' )
parser.add_argument('-refs', '--refs', default="ALL", dest='refs', type=str, help='References to analyse from the bam file. Default: All references are analysed.')
parser.add_argument('-SV', '--SV', default="INS,DEL,CLIPPING", dest='SV', type=str, help='Structural variants to look for within reads. INS: insertion, DEL: deletion, CLIPPING: clip reads. Default: INS, DEL, CLIPPING')
parser.add_argument('-minMAPQ', '--minMAPQ', default=20, dest='minMAPQ', type=int, help='Minimum mapping quality requires for each read.')
parser.add_argument('-minINDELlen', '--minINDELlen', default=50, dest='minINDELlen', type=int, help='Minimum indel length.')
parser.add_argument('-minCLIPPINGlen', '--minCLIPPINGlen', default=500, dest='minCLIPPINGlen', type=int, help='Minimum length of clipping region for each read.')
parser.add_argument('-maxBkpDist', '--maxBkpDist', default=50, dest='maxBkpDist', type=int, help='Maximum distance bewteen two adjacent breakpoints of the same cluster (applies only for those SVs which only one breakpoint (i.e. INS, CLIPPING)).')
parser.add_argument('-minRootClusterSize', '--minRootClusterSize', default=2, dest='minRootClusterSize', type=int, help='Minimum number of reads of the first cluster found (before extending it).')
parser.add_argument('-maxClusterDist', '--maxClusterDist', default=100, dest='maxClusterDist', type=int, help='Maximum distance bewteen different clusters.')

args = parser.parse_args()
bam = args.bam
normalBam = args.normalBam
threads = args.threads
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])
windowSize = args.windowSize
refs = args.refs
SV = args.SV
minMAPQ = args.minMAPQ
minINDELlen = args.minINDELlen
minCLIPPINGlen = args.minCLIPPINGlen
maxBkpDist = args.maxBkpDist
minRootClusterSize = args.minRootClusterSize
maxClusterDist = args.maxClusterDist


# If no reference is specified, get all that are present in the bam file.
if refs == "ALL":
	refs = bamtools.getREFS(bam)

# Convert comma-separated string inputs into lists:
targetSV = SV.split(',')
targetRefs = refs.split(',')


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

confDict['threads'] = threads
confDict['windowSize'] = windowSize
confDict['targetRefs'] = targetRefs
confDict['targetSV'] = targetSV

confDict['minMAPQ'] = minMAPQ
confDict['minINDELlen'] = minINDELlen
confDict['minCLIPPINGlen'] = minCLIPPINGlen

confDict['maxBkpDist'] = maxBkpDist
confDict['minRootClusterSize'] = minRootClusterSize
confDict['maxClusterDist'] = maxClusterDist

## 2. Launch structural variation (SV) caller
callerObj = callers.SVcaller_nano(mode, bam, normalBam, confDict, outDir)
callerObj.callSV()

print("***** Finished! *****")
print()
