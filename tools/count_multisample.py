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
parser.add_argument('vcf', help='Path to multisample VCF file')
parser.add_argument('outDir', help='Output directory')

## 2. Parse user input ##
args = parser.parse_args()
vcf = args.vcf
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('vcf: ', vcf)
print('outDir: ', outDir, "\n")


###############
## Count INS ##
###############

## Read VCF
VCF = formats.VCF()
VCF.read(vcf)

## Initialize dictionary
types = ['total', 'Alu', 'L1', 'SVA', 'ERVK', 'transduction', 'pseudogene', 'full', 'competent']

counts = {}

for sampleId in VCF.header.sampleIds:
    counts[sampleId] = {}

    for iType in types:
        counts[sampleId][iType] = 0

## Update counts
for variant in VCF.variants:

    ### Determine insertion type
    ## a) Alu
    if ('FAM' in variant.info) and (variant.info['FAM'] == 'Alu'):
        iType = 'Alu'

    ## b) L1
    elif ('FAM' in variant.info) and (variant.info['FAM'] == 'L1'):
        iType = 'L1'

    ## c) SVA
    elif ('FAM' in variant.info) and (variant.info['FAM'] == 'SVA'):
        iType = 'SVA'

    ## d) ERVK
    elif ('FAM' in variant.info) and (variant.info['FAM'] == 'ERVK'):
        iType = 'ERVK'

    ## f) Processed pseudogene
    elif (variant.info['ITYPE'] == 'pseudogene'): 
        iType = 'pseudogene'

    ### Update sample counts
    for sampleId in VCF.header.sampleIds:
        
        ## Skip if not carrier:
        if variant.format[sampleId]['GT'] in ['0|0', '0/0']:
            continue

        counts[sampleId]['total'] += 1
        counts[sampleId][iType] += 1

        ## Transduction
        if (variant.info['ITYPE'] == 'orphan') or (variant.info['ITYPE'] == 'partnered'):  
            counts[sampleId]['transduction'] += 1

        ## Full length L1 insertion
        if ('FAM' in variant.info) and (variant.info['FAM'] == 'L1') and ('IS_FULL' in variant.info):
            counts[sampleId]['full'] += 1

        # Competent L1 insertion
        if ('FAM' in variant.info) and (variant.info['FAM'] == 'L1') and ('COMPETENT' in variant.info):
            counts[sampleId]['competent'] += 1

## Organize counts into a dataframe
data = list(counts.values())
indexes = list(counts.keys())
columns = ['total', 'Alu', 'L1', 'SVA', 'ERVK', 'transduction', 'pseudogene', 'full', 'competent']

df = pd.DataFrame(data, index=indexes, columns=columns) 

## Write dataframe into tsv file
outFile = outDir + '/counts.tsv' 
df.to_csv(outFile, sep='\t', index=True)
