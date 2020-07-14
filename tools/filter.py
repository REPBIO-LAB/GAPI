## DEPENDENCIES ##
# External
import os
import sys
import argparse


# Internal
import formats

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('vcf', help='Path to VCF to be filtered')
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


################
## Filter VCF ##
################

## 1. Load VCF ##
inVCF = formats.VCF()
inVCF.read(vcf)

## 2. Filter VCF ##
## Create output VCF
outVCF = formats.VCF()
outVCF.header = inVCF.header 

## Apply filters to each variant in input VCF
for variant in inVCF.variants:

    ## Filter based on the insertion type
    if variant.info['ITYPE'] not in ['solo', 'partnered', 'orphan', 'pseudogene']:
        continue

    ## Filter solo based on family
    if (variant.info['ITYPE'] == 'solo') and (variant.info['FAM'] not in ['L1', 'Alu', 'SVA', 'ERVK', 'rRNA']):
        continue

    ## Filter retrotranspositions based on mechanism (filter out if not TPRT insertions)
    if (('FAM' in variant.info and variant.info['FAM'] in ['L1', 'Alu', 'SVA']) or variant.info['ITYPE'] in ['partnered', 'orphan']) and ('MECHANISM' not in variant.info or variant.info['MECHANISM'] != 'TPRT'):
        continue

    ## Filter ERVK if insertion over an annotated ERVK on the reference
    if ('FAM' in variant.info and variant.info['FAM'] == 'ERVK') and ('REP' in variant.info):

        filtered = False

        for rep, dist in zip(variant.info['REP'].split(','), variant.info['DIST'].split(',')):
            
            if rep == 'ERVK' and dist == '0':
                filtered = True
                break

        if filtered:
            continue
            
    ## Filter Solo Alu, L1 and SVA insertions based on their length (max length threshold)
    if (variant.info['ITYPE'] == 'solo') and (variant.info['FAM'] in ['L1', 'Alu', 'SVA']):

        if (variant.info['FAM'] == 'L1') and (int(variant.info['LEN']) > 6500):
            continue
        
        elif (variant.info['FAM'] == 'Alu') and (int(variant.info['LEN']) > 400): 
            continue

        elif (variant.info['FAM'] == 'SVA') and (int(variant.info['LEN']) > 4500): 
            continue

    ## Add variant passing all the filters to the output VCF
    outVCF.variants.append(variant)

## Write output VCF
fileName = os.path.basename(vcf) 
outName = os.path.splitext(fileName)[0] + '.FILTERED'
outVCF.write(inVCF.info_order, outName, outDir)
