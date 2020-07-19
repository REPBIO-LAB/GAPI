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
parser.add_argument('source', help='Path to FASTA containing source sequences')
parser.add_argument('ancestor', help='Path to FASTA containing common ancestor for the sources')
parser.add_argument('target', help='Path to FASTA containing target sequences for source inference')
parser.add_argument('--source-haplo', default=None, dest='sourceHaplo', help='Table containing haplotype for database of source elements')
parser.add_argument('--scores', default=None, dest='scores', help='Table containing score for each source haplotype position')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory')

## 2. Parse user input ##
args = parser.parse_args()
source = args.source
ancestor = args.ancestor
target = args.target
sourceHaplo = args.sourceHaplo
scores = args.scores
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]
version='0.1.0'

print()
print('***** ', scriptName, version, 'configuration *****')
print('source: ', source)
print('ancestor: ', ancestor)
print('target: ', target)
print('sourceHaplo: ', sourceHaplo)
print('scores: ', scores)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 1. Generate haplotypes using the ancestor as reference
###########################################################
msg = '1. Generate haplotypes using the ancestor as reference'
log.header(msg)

### 1.1 Input target sequences haplotype matrix
step = 'HAPLO'
msg = 'Compute target sequences haplotype matrix'
log.step(step, msg)

targetDir = outDir + '/target/'
targetHaplo = haplotype.compute_haplotypes(target, ancestor, targetDir)

### 1.2 Source elements haplotype matrix
if sourceHaplo is None:

    step = 'HAPLO'
    msg = 'Compute source sequences haplotype matrix'
    log.step(step, msg)

    # Compute source haplotype
    sourceDir = outDir + '/source/'
    sourceHaplo = haplotype.compute_haplotypes(source, ancestor, sourceDir)

else:
    step = 'HAPLO'
    msg = 'Pre-computed haplotype matrix for source sequences provided as input'
    log.step(step, msg)
    sourceHaplo = pd.read_csv(sourceHaplo, header=0, index_col=0, sep='\t')


## 2. Generate diagnostic score matrix 
#######################################
msg = '2. Generate diagnostic score matrix'
log.header(msg)

if scores is None:

    step = 'SCORE'
    msg = 'Generate diagnostic score matrix'
    log.step(step, msg)

    # For each position at each source sequence compute its diagnostic score
    scores = haplotype.diagnostic_scores(sourceHaplo)

    # Remove non-informative positions
    sourceHaplo, scores = haplotype.remove_noninformative(sourceHaplo, scores)

    # Save filtered score and haplotype matrix
    outFile = outDir + '/haplo_ref.tsv'
    sourceHaplo.to_csv(outFile, sep='\t')

    outFile = outDir + '/scores_ref.tsv'
    scores.to_csv(outFile, sep='\t')

else:
    step = 'SCORE'
    msg = 'Pre-computed diagnostic score matrix provided as input'
    log.step(step, msg)
    scores = pd.read_csv(scores, header=0, index_col=0, sep='\t')

## 3. Source haplotypes assignation 
#####################################
msg = '3. Source haplotypes assignation'
log.header(msg)

# For each input haplotype make list with more to less likely matching source haplotypes
assignations = haplotype.assign_haplotypes(targetHaplo, sourceHaplo, scores)

# Report source element assignations into a tsv file
outFile = outDir + '/assignations.tsv'
assignations.to_csv(outFile, sep='\t')
