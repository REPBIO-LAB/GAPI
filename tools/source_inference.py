## DEPENDENCIES ##
# External
import os
import sys
import argparse

# Internal
import haplotype
import log

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('source', help='Path to FASTA containing source sequences')
parser.add_argument('ancestor', help='Path to FASTA containing common ancestor for the sources')
parser.add_argument('target', help='Path to FASTA containing target sequences for source inference')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
source = args.source
ancestor = args.ancestor
target = args.target
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('source: ', source)
print('ancestor: ', ancestor)
print('target: ', target)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 1. Generate haplotypes using the ancestor as reference
###########################################################
msg = '1. Generate haplotypes using the ancestor as reference'
log.header(msg)

### 1.1. Source elements haplotype matrix
sourceDir = outDir + '/source/'
sourceHaplo = haplotype.compute_haplotypes(source, ancestor, sourceDir)

# Note: add option for using precomputed source matrix

### 1.2. Input target sequences haplotype matrix
targetDir = outDir + '/target/'
targetHaplo = haplotype.compute_haplotypes(target, ancestor, targetDir)

## 2. Generate diagnostic score matrix 
#######################################
msg = '2. Generate diagnostic score matrix'
log.header(msg)

# For each position at each source sequence compute its diagnostic score
scores = haplotype.diagnostic_scores(sourceHaplo)

# Note: add option for using precomputed scores matrix 

## 3. Source haplotypes assignation 
#####################################
msg = '3. Source haplotypes assignation'
log.header(msg)

# For each input haplotype make list with more to less likely matching source haplotypes
assignations = haplotype.assign_haplotypes(targetHaplo, sourceHaplo, scores)

## Report source element assignations into a tsv file
outFile = outDir + '/assignations.tsv'
assignations.to_csv(outFile, sep='\t')
