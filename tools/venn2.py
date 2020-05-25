## DEPENDENCIES ##
# External
import os
import sys
import argparse
import pandas as pd 
from matplotlib import pyplot as plt
from matplotlib_venn import venn2


def count_intersection_types(intersections):
    ''' 
    '''
    ## Count number of 10
    nb10 = len(intersections[(intersections.iloc[:, 0] == 1) & (intersections.iloc[:, 1] == 0)].index)

    ## Count number of 01
    nb01 = len(intersections[(intersections.iloc[:, 0] == 0) & (intersections.iloc[:, 1] == 1)].index)

    ## Count number of 11
    nb11 = len(intersections[(intersections.iloc[:, 0] == 1) & (intersections.iloc[:, 1] == 1)].index)
    
    return nb10, nb01, nb11

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('intersections', help='')
parser.add_argument('outName', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
intersections = args.intersections
outName = args.outName
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('intersections: ', intersections)
print('outName: ', outName)
print('outDir: ', outDir, "\n")

##########
## CORE ##
##########

## Read intersection matrix
intersections = pd.read_csv(intersections, sep='\t', index_col=0)

## Count intersection types (10, 01, 11)
nb10, nb01, nb11 = count_intersection_types(intersections)

## Make Venn diagram
labels = list(intersections)
fig = plt.figure(figsize=(5,5))
venn2(subsets = (nb10, nb01, nb11), set_labels = labels)

## Save figure
fileName = outDir + '/' + outName + '_venn2.pdf'
plt.savefig(fileName)

fileName = outDir + '/' + outName + '_venn2.svg'
plt.savefig(fileName)

## Make summary matrix
total = float(nb10 + nb01 + nb11)
totalA = float(nb10 + nb11)
totalB = float(nb01 + nb11)

matrix = [['caller', 'perShared', 'percTotal'], 
        [labels[0], nb11/totalA*100, totalA/total*100], 
        [labels[1], nb11/totalB*100, totalB/total*100]]

## Convert matrix to dataframe
df = pd.DataFrame(matrix[1:], columns=matrix[0])
df = df.set_index('caller')

outFile = outDir + '/' + outName + '_summary.tsv'
df.to_csv(outFile, sep='\t', index=True)
