## DEPENDENCIES ##
# External
import os
import sys
import argparse
import pandas as pd 

# Internal
import formats


######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('paths', help='Two columns text file: sampleId VCF_path')
parser.add_argument('outDir', help='Output directory')

## 2. Parse user input ##
args = parser.parse_args()
paths = args.paths
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('paths: ', paths)
print('outDir: ', outDir, "\n")


###############
## Count INS ##
###############

counts = {}
samples = []

paths = open(paths)

# For each sample
for line in paths:
            
    sampleId, VCF_path = line.split() 
    samples.append(sampleId)
    
    ## Read VCF
    VCF = formats.VCF()
    VCF.read(VCF_path)

    ## Initialize counters
    counts[sampleId] = {}
    counts[sampleId]['total'] = 0
    counts[sampleId]['nbAlu'] = 0
    counts[sampleId]['nbL1'] = 0
    counts[sampleId]['nbSVA'] = 0
    counts[sampleId]['nbERVK'] = 0
    counts[sampleId]['nbTd'] = 0
    counts[sampleId]['nbPseudo'] = 0
    counts[sampleId]['nbRNA'] = 0 

    ## Count event types
    for variant in VCF.variants:
        counts[sampleId]['total'] += 1

        ## a) Alu
        if ('FAM' in variant.info) and (variant.info['FAM'] == 'Alu'):
            counts[sampleId]['nbAlu'] += 1

        ## b) Solo L1
        elif ('FAM' in variant.info) and (variant.info['FAM'] == 'L1') and (variant.info['ITYPE'] == 'solo'):
            counts[sampleId]['nbL1'] += 1

        ## c) SVA
        elif ('FAM' in variant.info) and (variant.info['FAM'] == 'SVA'):
            counts[sampleId]['nbSVA'] += 1

        ## d) ERVK
        elif ('FAM' in variant.info) and (variant.info['FAM'] == 'ERVK'):
            counts[sampleId]['nbERVK'] += 1

        ## e) Transductions
        elif (variant.info['ITYPE'] == 'orphan') or (variant.info['ITYPE'] == 'partnered'):            
            counts[sampleId]['nbTd'] += 1

        ## f) Processed pseudogene
        elif (variant.info['ITYPE'] == 'pseudogene'): 
            counts[sampleId]['nbPseudo'] += 1

        ## g) rRNA
        elif ('FAM' in variant.info) and (variant.info['FAM'] == 'rRNA'):
            counts[sampleId]['nbRNA'] += 1

## Organize counts into a dataframe
data = list(counts.values())
indexes = list(counts.keys())
columns = ['total', 'nbAlu', 'nbL1', 'nbSVA', 'nbERVK', 'nbTd', 'nbPseudo', 'nbRNA']

df = pd.DataFrame(data, index=indexes, columns=columns) 

## Write dataframe into tsv file
outFile = outDir + '/counts.tsv' 
df.to_csv(outFile, sep='\t', index=True)