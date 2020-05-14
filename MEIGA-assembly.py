## DEPENDENCIES ##
# External
import os
import sys
import argparse

# Internal
import formats
import alignment
import unix
import sequences

###############
## Functions ##
###############


def search4polyA(sequence):
    '''
    Search for poly(A) at sequence end 
    '''
    ## Configuration for monomere search:
    windowSize = 8
    maxWindowDist = 2
    minMonomerSize = 10
    minPurity = 80  
    maxDist2Ends = 10 

    ## Seach poly(A) monomers
    targetMonomer = 'A'
    monomersA = sequences.find_monomers(sequence, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    if not monomersA:
        return False, None
        
    ## Select monomer closest to sequence end
    candidate = monomersA[-1]

    ## Filter out monomer if more than Xbp from end
    seqLen = len(sequence)
    dist2end = seqLen - candidate.end

    if dist2end <= 30:
        polyA = True

    else:
        polyA = False 

    return polyA, candidate

def search4polyT(sequence):
    '''
    Search for poly(T) at sequence begin 
    '''
    ## Configuration for monomere search:
    windowSize = 8
    maxWindowDist = 2
    minMonomerSize = 10
    minPurity = 80  
    maxDist2Ends = 10 

    ## Seach poly(T) monomers
    targetMonomer = 'T'
    monomersT = sequences.find_monomers(sequence, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    if not monomersT:
        return False, None

    ## Select monomer closest to sequence beg        
    candidate = monomersT[0]

    ## Filter out monomer if more than Xbp from beg
    dist2beg = candidate.beg

    if dist2beg <= 30:
        polyT = True

    else:
        polyT = False 

    return polyT, candidate

def candidate_MEI(VCF):
    '''
    '''
    tails = {}

    ## Create VCF with filtered calls
    filteredVCF = formats.VCF()
    filteredVCF.header = VCF.header

    ## For each variant
    for variant in VCF.variants:
        
        ## Filter calls distinct from INS
        if variant.info['SVTYPE'] != 'INS':
            continue

        ## Search for poly(A) tail at inserted sequence
        polyA, monomerA = search4polyA(variant.alt)

        ## Search for poly(T) tail at inserted sequence
        polyT, monomerT = search4polyT(variant.alt)

        ## Filter calls if poly(A) nor poly(T) found
        if not polyA and not polyT:
            continue

        ## Add variant passing all the filters
        filteredVCF.add(variant)

        ## Add polyA/T information to the dict
        insId = variant.chrom + '_' + str(variant.pos) 
        tails[insId] = {}

        # a) PolyA and polyT tracts identified
        if polyA and polyT:
            tails[insId]['polyA'] = monomerA 
            tails[insId]['polyT'] = monomerT

        # b) PolyA identified
        elif polyA:
            tails[insId]['polyA'] = monomerA

        # c) PolyT identified
        else:
            tails[insId]['polyT'] = monomerT

    return filteredVCF, tails

def call_MEI(vcf, tails):
    '''
    '''
    # Write inserted sequences into fasta file
    fastaPath = outDir + '/insertions.fa'
    fasta = ins2fasta(vcf, outDir)
    fasta.write(fastaPath)

    ## 4. Realign inserted sequence against retrotransposon consensus database
    ## 4.1 Create index or consensus sequences
    fileName = 'consensus'  
    consensusIndex = alignment.index_minimap2(consensus, fileName, outDir)

    ## 4.2 Realign inserted sequences against consensus:
    pafPath = alignment.alignment_minimap2(fastaPath, consensusIndex, 'hits2consensus', 1, outDir)
    paf = formats.PAF()
    paf.read(pafPath)

    ## 4.2 Generate a single paf object per inserted sequence:
    pafDict = {}

    for hit in paf.alignments:

        if hit.qName not in pafDict:
            pafDict[hit.qName] = formats.PAF()
    
        pafDict[hit.qName].alignments.append(hit)

    ## Add PolyA/T calls to the paf
    for insId in pafDict.keys():

        insLen = len(fasta.seqDict[insId])

        if 'polyA' in tails[insId]:
            polyA = tails[insId]['polyA']
            fields = [insId, insLen, polyA.beg, polyA.end, None, 'polyA', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            pafDict[hit.qName].alignments.append(hit) 

        if 'polyT' in tails[insId]:
            polyT = tails[insId]['polyT']

            fields = [insId, insLen, polyT.beg, polyT.end, None, 'polyT', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            pafDict[hit.qName].alignments.append(hit) 

    for insId in pafDict:
        chain = pafDict[insId].chain(20, 50)
        alignments = ';'.join([str(hit.qBeg) + '_' + str(hit.qEnd) + '_' + hit.tName for hit in chain.alignments])
        print('chain: ', chain.interval(), chain.perc_query_covered(), alignments, len(fasta.seqDict[insId]), fasta.seqDict[insId])



def ins2fasta(vcf, outDir):
    '''
    Write inserted sequences into a fasta

    Input:
        1. vcf: Path to VCF file
        2. outDir: Output directory

    Output: 
        1. fastaPath: path to fasta object inserted sequences
    '''     
    ## 1. Initialize fasta object
    fasta = formats.FASTA() 

    ## 2. Collect inserted sequences
    for variant in vcf.variants:
 
        insId = variant.chrom + '_' + str(variant.pos) 
        seq = variant.alt

        fasta.seqDict[insId] = seq

    return fasta

######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('vcf', help='Path to VCF file with assembly-based SV calls')
parser.add_argument('consensus', help='Path to FASTA file containing consensus retrotransposon sequences')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
vcf = args.vcf
consensus = args.consensus
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('vcf: ', vcf)
print('consensus: ', consensus)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 1. Read VCF
VCF = formats.VCF()
VCF.read(vcf)

## 2. Filter VCF by selecting retrotransposition insertion candidates 
# (inserted sequences with polyA/T tails at their ends)
filteredVCF, tails = candidate_MEI(VCF)
print('filteredVCF: ', filteredVCF, len(filteredVCF.variants))

## 3. Do MEI calling for candidate insertions
call_MEI(filteredVCF, tails)
 