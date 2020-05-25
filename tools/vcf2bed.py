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
parser.add_argument('VCF_path', help='Path to VCF file')
parser.add_argument('outName', help='Output file name')
parser.add_argument('outDir', help='Output directory')

## 2. Parse user input ##
args = parser.parse_args()
VCF_path = args.VCF_path
outName = args.outName
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('VCF_path: ', VCF_path)
print('outName: ', outName)
print('outDir: ', outDir, "\n")


###########################
## Generate bed from vcf ##
###########################

## Initialize BED
BED = formats.BED()
BED.lines = []
BED.structure = 'List'

## Read VCF
VCF = formats.VCF()
VCF.read(VCF_path)

## Add variants to BED
for variant in VCF.variants:

    if 'FAM' in variant.info:
    
        ## Create bed entry 
        fields = [variant.chrom, variant.pos, variant.pos + 1, variant.info['FAM']]
        header = ['ref', 'beg', 'end', 'name']
        entry = formats.BED_entry(fields, header)
        entry.name = variant.info['FAM']
        BED.lines.append(entry)


## Write BED
outPath = outDir + '/' + outName + '.bed'
BED.write(outPath)
