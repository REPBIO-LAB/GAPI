## DEPENDENCIES ##
# External
import argparse
import sys
import os

# Internal
import formats
import log

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('source', help='Path to FASTA containing source sequences')
parser.add_argument('length', type=int, help='Target length')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory')

## 2. Parse user input ##
args = parser.parse_args()
source = args.source
length = args.length
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('source: ', source)
print('length: ', length)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 1. Read input fasta file
############################
fasta = formats.FASTA()
fasta.read(source)

## 2. Trim sequences on the fasta file
#######################################
for seqId, seq in fasta.seqDict.items():
    pos = len(seq)- length
    fasta.seqDict[seqId] = seq[pos:]

## 3. Write sequence
#####################
outFile = outDir + '/subsequences.' + str(length) + '.tsv'
fasta.write(outFile)