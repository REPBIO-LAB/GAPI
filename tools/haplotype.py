## DEPENDENCIES ##
# External
import os
import sys
import argparse
import pandas as pd 

# Internal
import formats
import alignment
import bamtools
import unix


def VCF2haplotype(vcf, lenght):
    '''
    Generate haplotype vector from a VCF containing SNPs calls.

    Input: 
        1. vcf: SNP and INDEL calls in VCF format
        2. length: length of reference sequence

    Output:
        1. haplotype: haplotype vector containing the status for each reference position 
                      (0, no change; 1-12: one of the 12 possible single base substitutions)
    '''
    ## 1. Create dictionary with substitution type codes:
    codes = {
            'A->T': 1,
            'A->G': 2,
            'A->C': 3,
            'T->A': 4,
            'T->G': 5,
            'T->C': 6,
            'G->A': 7,
            'G->T': 8,
            'G->C': 9,
            'C->A': 10,
            'C->T': 11,
            'C->G': 12           
            }

    ## 2. Initialize haplotype vector as always reference (0s)
    haplotype = [0] * lenght

    ## 3. Update haplotype vector with identified SNPs
    # Read VCF
    VCF = formats.VCF()
    VCF.read(vcf)

    # Iterate over each identified variant
    for variant in VCF.variants:

        ## Skip reference calls and indels
        if variant.alt == '.' or 'INDEL' in variant.info:
            continue
    
        ## Add SNP call to the vector
        key = variant.ref + '->' + variant.alt
        code = codes[key]
        index = variant.pos - 1 # Convert 1-based (vcf) to 0-based (vector)

        haplotype[index] = code

    return haplotype


######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('vcf', help='Path to VCF file')
parser.add_argument('fam', help='Target MEI family')
parser.add_argument('consensus', help='MEI family consensus sequence')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
vcf = args.vcf
fam = args.fam
consensus = args.consensus
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('vcf: ', vcf)
print('fam: ', fam)
print('consensus: ', consensus)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########
## Compute consensus length
fasta = formats.FASTA()             
fasta.read(consensus)
length = len(list(fasta.seqDict.values())[0])

## Read VCF
VCF = formats.VCF()
VCF.read(vcf)

## Compute variant haplotypes
haplotypes = {}

for variant in VCF.variants:

    # Skip not target insertions
    if (variant.info['ITYPE'] != 'solo') or ('FAM' not in variant.info) or (variant.info['FAM'] != fam):
        continue

    ## 1. Create tmp output directory
    seqId = variant.chrom + '_' + str(variant.pos) + '_' + variant.info['FAM'] 
    tmpDir = outDir + '/' + seqId
    unix.mkdir(tmpDir)

    ## 2. Create Fasta containing inserted sequence
    fasta = formats.FASTA()             
    fasta.seqDict[seqId] = variant.info['INSEQ']
    fastaPath = tmpDir + '/' + seqId + '.fa'
    fasta.write(fastaPath)

    ## 3. Map inserted sequence against consensus
    SAM = alignment.alignment_bwa(fastaPath, consensus, seqId, 1, tmpDir)

    ## 4. SAM to BAM conversion and sorting
    BAM = bamtools.SAM2BAM(SAM, tmpDir)

    ## 5. Call SNPs and INDELS
    calls = tmpDir + '/' + seqId + '.vcf' 
    command = 'samtools mpileup -f ' + consensus + ' -g ' + BAM + ' | bcftools call --ploidy 1 -c - > ' + calls
    os.system(command) 

    ## 6. Convert SNP calls into vector representation
    haplotype = VCF2haplotype(calls, length)
    haplotypes[seqId] = haplotype

    ## Cleanup
    unix.rm([tmpDir])

## Generate dataframe containing haplotypes
colNames = list(range(1, length + 1))
haplotypesDf = pd.DataFrame.from_dict(haplotypes, orient='index', columns=colNames)

## Write haplotype table into tsv file as output
outFile = outDir + '/' + fam + '_haplotypes.tsv' 
haplotypesDf.to_csv(outFile, sep='\t', index=True)