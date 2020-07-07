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
alleleCounts = {}

## Update counts
for variant in VCF.variants:

    ## Collect basic info
    coord = variant.chrom + ':' + str(variant.pos) + '-' + str(variant.pos + 1)
    iType = variant.info['ITYPE']
    family = variant.info['FAM'] if 'FAM' in variant.info else 'None'

    ## Compute insertion allele count
    alleleCount = 0

    for sampleId in VCF.header.sampleIds:
        
        ## Not carrier:
        if variant.format[sampleId]['GT'] in ['0|0', '0/0']:
            continue

        ## a) Homozigous carrier
        elif variant.format[sampleId]['GT'] in ['1|1', '1/1']:
            alleleCount += 2 
            
        ## b) Heterozigous carrier or haploid chromosome
        elif variant.format[sampleId]['GT'] in ['1|0', '1/0', '0|1', '0/1', '1']:
            alleleCount += 1 

    alleleCounts[coord] = [iType, family, alleleCount]

## Organize counts into a dataframe
df = pd.DataFrame(alleleCounts).T
df.columns = ['type', 'fam', 'alleleCount']

## Write dataframe into tsv file
outFile = outDir + '/allele_counts.tsv' 
df.to_csv(outFile, sep='\t', index=True)
