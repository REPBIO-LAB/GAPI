## DEPENDENCIES ##
# External
import os
import sys
import argparse

# Internal
import formats
import clustering


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


################
## Merge VCFs ##
################

## 1. Read VCFs ##
paths = open(paths)
VCFs = []
samples = []

# For line in the file
for line in paths:
            
    sampleId, VCF_path = line.split() 

    ## Add sample to the list
    samples.append(sampleId)

    ## Read VCF
    VCF = formats.VCF()
    VCF.read(VCF_path)

    ## Set Donor id
    VCF.sampleId = sampleId

    for variant in VCF.variants:
        variant.sampleId = sampleId

    ## Add VCF to the list
    VCFs.append(VCF)

## 2. Load VCFs into a bin database ##
wgBinDb = formats.INS2binDb(VCFs, VCFs[0].header.refLengths, 1)

## 3. Do events clustering ##
clusters = []

# For each ref
for ref in wgBinDb:

    # For each event type
    for eventType in wgBinDb[ref].eventTypes:
        
        # Do clustering
        clusters = clusters + clustering.reciprocal_overlap_clustering(wgBinDb[ref], 1, 1, [eventType], 50, 'INS_VCF')

## 4. Create consensus events ##
cVariants = []

for cluster in clusters:
    cVariants.append(cluster.consensus())

## 5. Create VCF containing consensus variants ##
outVCF = formats.VCF()

## Define header
outVCF.header = VCFs[0].header
outVCF.header.info['SAMPLES'] = ['.', 'String', 'Comma separated list of samples where the variant was identified']

## Add consensus variants
outVCF.variants = cVariants

## Write output VCF
order = VCFs[0].info_order + ['SAMPLES']
outName = 'MERGED'
outVCF.write(order, outName, outDir)
