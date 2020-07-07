## DEPENDENCIES ##
# External
import os
import sys
import argparse
import pandas as pd 

# Internal
import formats
import structures
import clustering

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('paths', help='')
parser.add_argument('refLengths', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
paths = args.paths
refLengths = args.refLengths
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('paths: ', paths)
print('refLengths: ', refLengths)
print('outDir: ', outDir, "\n")

##########
## MAIN ##
##########

## 0. Read reference lengths
lenDict = {}
refLengths = open(refLengths)

for line in refLengths:
            
    ref, length = line.split() 
    lenDict[ref] = int(length)

## 1. Read Bed files
allCallers = []

paths = open(paths)
beds = []

# For line in the file
for line in paths:
            
    caller, bedPath = line.split() 
    allCallers.append(caller)

    BED = formats.BED()
    BED.read(bedPath, 'nestedDict', None)

    ## Add caller info to each bed entry (temporal, improve later)
    for ref in BED.lines:
        for name in BED.lines[ref]:
            for entry in BED.lines[ref][name]:
                entry.caller = caller

    ## Add bed to the list
    beds.append(BED)

## 2. Merge Bed files
mergedDict = {}

for ref in lenDict.keys():

    bed2merge = []

    for BED in beds:
        if ref in BED.lines:
            bed2merge.append(BED.lines[ref])

    mergedDict[ref] = structures.merge_dictionaries(bed2merge)
    
## 3. Organize Bed entries into a bin database
binDb = structures.create_bin_database(lenDict, mergedDict)

## 4. Cluster Bed entries
clusters = []

# For each ref
for ref in binDb:

    # For each event type
    for eventType in binDb[ref].eventTypes:
        
        # Do clustering
        clusters = clusters + clustering.reciprocal_overlap_clustering(binDb[ref], 1, 1, [eventType], 100, 'INS_VCF')

## 5. Generate matrix
matrix = [['id', 'name'] + allCallers]

for cluster in clusters:

    fields = []
    variant = cluster.variants[0]
    ID =  variant.ref + '_' + str(variant.beg) + '_' + str(variant.end)
    fields.append(ID)
    fields.append(variant.optional['name'])

    callers = [ variant.caller for variant in cluster.variants]

    for caller in allCallers:
        if caller in callers:
            fields.append(1)
        else:
            fields.append(0)
    
    matrix.append(fields)

## Convert to pandas dataframe
df = pd.DataFrame(matrix[1:], columns=matrix[0])
df = df.set_index('id')

## Write haplotype table into tsv file as output
outFile = outDir + '/intersection.tsv' 
df.to_csv(outFile, sep='\t', index=True)