## DEPENDENCIES ##
# External
import os
import sys
import argparse
import pandas as pd 

# Internal
import haplotype
import log

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('sequences', help='Path to FASTA containing input sequences')
parser.add_argument('reference', help='Path to FASTA containing reference sequence')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory')

## 2. Parse user input ##
args = parser.parse_args()
sequences = args.sequences
reference = args.reference
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]
version='0.1.0'

print()
print('***** ', scriptName, version, 'configuration *****')
print('sequences: ', sequences)
print('reference: ', reference)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 1. Generate haplotypes comparing with the reference
tmpDir = outDir + '/tmp/'
haplo = haplotype.compute_haplotypes(sequences, reference, tmpDir)

## 2. Write Haplotypes matrix into output file
outFile = outDir + '/haplotypes.tsv'
haplo.to_csv(outFile, sep='\t')
