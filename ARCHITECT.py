#!/usr/bin/env python
#coding: utf-8

####################
## Import modules ##
####################

# External
import argparse
import sys
import os

# Internal
import callers
import log
import bamtools

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description= "Call structural variants (SV) from Nanopore/Pacbio whole genome sequencing data. Two running modes: 1) SINGLE: individual sample; 2) PAIRED: tumour and matched normal sample")
parser.add_argument('bam', help='Input bam file. Will correspond to the tumour sample in the PAIRED mode')

## General
parser.add_argument('--normalBam', default="NA", dest='normalBam', help='Matched normal bam file. If provided ARCHITECT will run in PAIRED mode')
parser.add_argument('-p', '--processes', default=1, dest='processes', type=int, help='Number of processes. Default: 1')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='Output directory. Default: current working directory')

## BAM processing
parser.add_argument('--targetBins', default=None, dest='targetBins', type=str, help='Bed file containing target genomic bins for SV calling. Overrides --binSize and --refs. Default: None')
parser.add_argument('-bS', '--binSize', default=1000000, dest='binSize', type=int, help='Input bams will be analised in genomic bins of this size. Default: 1000000')
parser.add_argument('--refs', default="ALL", dest='refs', type=str, help='Comma separated list of target references to call SV (i.e. 1,2,3,X). Default: All references included in the bam file')
parser.add_argument('--SV', default="INS,DEL,CLIPPING", dest='SV', type=str, help='Comma separated list of SV event types to collect (INS, DEL and CLIPPING). Default: INS,DEL,CLIPPING')

## Filtering thresholds
parser.add_argument('--minMAPQ', default=20, dest='minMAPQ', type=int, help='Minimum mapping quality required for each read. Default: 20')
parser.add_argument('--minINDELlen', default=50, dest='minINDELlen', type=int, help='Minimum indel length. Default: 50')
parser.add_argument('--minCLIPPINGlen', default=500, dest='minCLIPPINGlen', type=int, help='Minimum clipped sequence length for each read. Default: 500')
parser.add_argument('--minRootClusterSize', default=2, dest='minRootClusterSize', type=int, help='Minimum number of reads composing a root cluster (before extension step). Default: 2')
parser.add_argument('--maxBkpDist', default=50, dest='maxBkpDist', type=int, help='Maximum distance bewteen two adjacent breakpoints for INS and CLIPPING clustering. Default: 50')
parser.add_argument('--minPercOverlap', default=70, dest='minPercRcplOverlap', type=int, help='Minimum percentage of reciprocal overlap for DEL clustering. Default: 50')
parser.add_argument('--maxClusterDist', default=100, dest='maxClusterDist', type=int, help='Maximum distance bewteen different clusters for metaclustering. Default: 100')
parser.add_argument('--clusterFilters', default="NBREADS,CV,OUTLIERS", dest='clusterFilters', type=str, help='Comma separated list of cluster filters to apply (minimum number of reads, minimum Coefficient of Variation and minimum percentage of outliers). Default: NBREADS,CV,OUTLIERS')
parser.add_argument('--minClusterSize', default=2, dest='minClusterSize', type=int, help='Minimum number of reads composing a cluster (after extension step). Default: 2')
parser.add_argument('--maxClusterCV', default=15, dest='maxClusterCV', type=int, help='Maximum Coefficient of Variation of a cluster. Default: 15')
parser.add_argument('--maxOutliers', default=0.5, dest='maxOutliers', type=int, help='Maximum percentage of events supporting a cluster that have been removed during the polish step. Default: 0.5')




## 2. Parse user´s input and initialize variables ##
args = parser.parse_args()

## General
bam = args.bam
normalBam = args.normalBam
processes = args.processes
outDir = args.outDir

## BAM processing
targetBins = args.targetBins
binSize = args.binSize
refs = args.refs
SV = args.SV

## Filtering thresholds
minMAPQ = args.minMAPQ
minINDELlen = args.minINDELlen
minCLIPPINGlen = args.minCLIPPINGlen
maxBkpDist = args.maxBkpDist
minRootClusterSize = args.minRootClusterSize
minPercRcplOverlap = args.minPercRcplOverlap
maxClusterDist = args.maxClusterDist
clusterFilters = args.clusterFilters
minClusterSize = args.minClusterSize
maxClusterCV = args.maxClusterCV
maxOutliers = args.maxOutliers

# If no reference is specified, get all that are present in the bam file.
if refs == "ALL":
	refs = bamtools.getREFS(bam)

# Convert comma-separated string inputs into lists:
targetSV = SV.split(',')
targetRefs = refs.split(',')

## Determine running mode:
mode = "SINGLE" if normalBam == "NA" else "PAIRED"

##############################################
## Display configuration to standard output ##
##############################################
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]
version='0.0.2'

print()
print('***** ', scriptName, version, 'configuration *****')
print('** General **')
print('mode: ', mode)
print('bam: ', bam)
print('normalBam: ', normalBam)
print('processes: ', processes)
print('outDir: ', outDir, "\n")

print('** BAM processing **')
print('targetBins: ', targetBins)
print('binSize: ', binSize)
print('targetRefs: ', refs)
print('targetSVs: ', SV, "\n")

print('** Filtering thresholds **')
print('minMAPQ: ', minMAPQ)
print('minINDELlength: ', minINDELlen)
print('minCLIPPINGlength: ', minCLIPPINGlen)
print('minRootClusterSize: ', minRootClusterSize)
print('maxBkpDistance: ', maxBkpDist)
print('minPercOverlap: ', minPercRcplOverlap)
print('maxClusterDistance: ', maxClusterDist)
print('clusterFilters: ', clusterFilters)
print('maxClusterDistance: ', minClusterSize)
print('maxClusterCV: ', maxClusterCV)
print('maxOutliers: ', maxOutliers, "\n")



print('***** Executing ', scriptName, '.... *****', "\n")


##########
## CORE ## 
##########
## 1. Create configuration dictionary
confDict = {}

## General
confDict['processes'] = processes

## BAM processing
confDict['targetBins'] = targetBins
confDict['binSize'] = binSize
confDict['targetRefs'] = targetRefs
confDict['targetSV'] = targetSV

## Filtering thresholds
confDict['minMAPQ'] = minMAPQ
confDict['minINDELlen'] = minINDELlen
confDict['minCLIPPINGlen'] = minCLIPPINGlen
confDict['maxBkpDist'] = maxBkpDist
confDict['minRootClusterSize'] = minRootClusterSize
confDict['minPercRcplOverlap'] = minPercRcplOverlap
confDict['maxClusterDist'] = maxClusterDist
confDict['clusterFilters'] = clusterFilters
confDict['minClusterSize'] = minClusterSize
confDict['maxClusterCV'] = maxClusterCV
confDict['maxOutliers'] = maxOutliers

## 2. Launch structural variation (SV) caller
callerObj = callers.SVcaller_nano(mode, bam, normalBam, confDict, outDir)
callerObj.callSV()

print('***** Finished! *****')
print()
