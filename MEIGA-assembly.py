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

def call_MEI_candidate(VCF):
    '''
    '''
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

    return filteredVCF

def call_MEI(vcf):
    '''
    '''
    ## 1. Write inserted sequences into fasta file
    fastaPath = outDir + '/insertions.fa'
    fasta = ins2fasta(vcf, outDir)
    fasta.write(fastaPath)

    ## 2. Create index or consensus sequences
    fileName = 'consensus'  
    consensusIndex = alignment.index_minimap2(consensus, fileName, outDir)

    ## 3. Realign inserted sequences against consensus:
    pafPath = alignment.alignment_minimap2(fastaPath, consensusIndex, 'hits2consensus', 1, outDir)
    paf = formats.PAF()
    paf.read(pafPath)

    ## 4. Generate a single paf object per inserted sequence:
    pafDict = group_alignments(paf)

    ## 5. Resolve structure for each insertion with matches on retrotransposon consensus sequences
    structures = {}

    for insId in pafDict:
        structures[insId] = MEI_structure(pafDict[insId], fasta.seqDict[insId])

    ## 6. Realign full or bits of the inserted sequence into the reference genome to:
    # - Determine partnered transduction procedence
    # - Identify orphan transductions
    # - Identify processed pseudogene insertions

    ## 7. Generate output VCF containing MEI calls

    ## Add MEI specific fields to the VCF header

    ## For each INS calls a

def MEI_structure(paf, insertSeq):
    '''
    '''
    structure = {}
    seqLen = len(insertSeq)

    ## 1. Chain alignments
    structure['CHAIN'] = paf.chain(20, 50)

    ## 2. Determine insertion family
    families = list(set([hit.tName.split('|')[1] for hit in structure['CHAIN'].alignments]))

    # a) L1 insertion
    if 'L1' in families:
        structure['FAM'] = 'L1'

    # b) SVA insertion
    elif 'SVA' in families:
        structure['FAM'] = 'SVA'

    # c) Alu insertion
    elif 'Alu' in families:
        structure['FAM'] = 'Alu'

    ## 3. Search for polyA/T tails at unresolved insert ends and determine insertion type
    structure['ITYPE'] = 'solo'
    rtBeg, rtEnd = structure['CHAIN'].interval()

    ## Set parameters for monomer search
    windowSize = 8
    maxWindowDist = 2
    minMonomerSize = 10
    minPurity = 80  

    ### 3.1 PolyA search
    ## Search for monomers
    targetSeq = insertSeq[rtEnd:]
    monomersA = sequences.find_monomers(targetSeq, 'A', windowSize, maxWindowDist, minMonomerSize, minPurity)

    ## Map to insert sequence coordinates
    for monomer in monomersA:
        monomer.beg = monomer.beg + rtEnd
        monomer.end = monomer.end + rtEnd

    ## Make polyA calls
    structure['POLYA'] = 0
    structure['STRAND'] = None

    # a) Single polyA
    if (len(monomersA) == 1):
        dist2rt = monomersA[0].beg - rtEnd
        dist2end = seqLen - monomersA[0].end

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'solo'
            structure['POLYA'] = 1
            structure['STRAND'] = '+'

    # b) Multiple polyA (Transduction candidate)
    # AAAAAAAAAAAAAA-------------AAAAAAAAAAAAAA
    elif (len(monomersA) > 1):
        dist2rt = monomersA[0].beg - rtEnd
        dist2end = seqLen - monomersA[-1].end

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'Partnered'
            structure['POLYA'] = len(monomersA)
            structure['STRAND'] = '+'

    ## 3.2 PolyT search 
    ## Search for monomers
    targetSeq = insertSeq[:rtBeg]
    monomersT = sequences.find_monomers(targetSeq, 'T', windowSize, maxWindowDist, minMonomerSize, minPurity)

    ## Make polyT calls
    structure['POLYT'] = 0

    # a) Single polyT
    if (len(monomersT) == 1):
        dist2end = monomersT[0].beg
        dist2rt = rtBeg - monomersT[0].end 

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'solo'
            structure['POLYT'] = 1
            structure['STRAND'] = '-'

    # b) Multiple polyT (Partnered transduction candidate)
    # TTTTTTTTTTTTT-------------TTTTTTTTTTTTT
    elif (len(monomersT) == 2):
        dist2end = monomersT[0].beg
        dist2rt = rtBeg - monomersT[-1].end 

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'Partnered'
            structure['POLYT'] = len(monomersT)
            structure['STRAND'] = '-'
        
    ## 3.3 Add polyA/T annotation to the chain
    # A) PolyA found
    if (structure['POLYA'] != 0) and (structure['POLYT'] == 0):

        # a) Solo
        if (structure['POLYA'] == 1):
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersA[0].beg, monomersA[0].end, None, 'polyA', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.append(hit) 

        # b) Partnered
        else:
            # First polyA
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersA[0].beg, monomersA[0].end, None, 'polyA', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.append(hit)

            # Partnered 
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersA[0].end, monomersA[-1].beg, None, 'partnered', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.append(hit)       

            # Second polyA
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersA[-1].beg, monomersA[-1].end, None, 'polyA', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.append(hit)            

    # B) PolyT found
    elif (structure['POLYT'] != 0) and (structure['POLYA'] == 0): 

        # a) Solo
        if (structure['POLYT'] == 1):
  
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersT[0].beg, monomersT[0].end, None, 'polyT', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.insert(0, hit)

        # b) Partnered
        else:
            # Second polyT
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersT[-1].beg, monomersT[-1].end, None, 'polyT', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.insert(0, hit)

            # Partnered 
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersT[0].end, monomersT[-1].beg, None, 'partnered', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.insert(0, hit)

            # First polyT
            fields = [structure['CHAIN'].alignments[0].qName, seqLen, monomersT[0].beg, monomersT[0].end, None, 'polyT', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            structure['CHAIN'].alignments.insert(0, hit)

    # Compute % of inserted resolved
    structure['PERC-RESOLVED'] = structure['CHAIN'].perc_query_covered()

    ## 4. Apply filters 
    failed = []

    # 4.1 Percentage resolved filter

    if structure['PERC-RESOLVED'] < 60:
        failed.append('PERC-RESOLVED')

    # 4.2 Length filtering for solo insertions
    if structure['ITYPE'] == 'solo':

        if (structure['FAM'] == 'L1') and (seqLen > 6500):
            failed.append('LEN')
        
        elif (structure['FAM'] == 'Alu') and (seqLen > 400): 
            failed.append('LEN')

        elif (structure['FAM'] == 'SVA') and (seqLen > 4500): 
            failed.append('LEN')

    # a) Insertions passes all the filters
    if not failed:
        structure['PASS'] = True

    # b) At least one failed filter
    else:
        structure['PASS'] = False

    print('RESULT: ', structure['FAM'], structure['ITYPE'], structure['POLYA'], structure['POLYT'], structure['PERC-RESOLVED'], structure['PASS'], failed, insertSeq)

    return structure


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

def group_alignments(paf):
    '''

    '''     
    pafDict = {}

    ## For each hit
    for hit in paf.alignments:

        # Initialize paf object for this inserted sequence
        if hit.qName not in pafDict:
            pafDict[hit.qName] = formats.PAF()
    
        # Add hit to the corresponding paf
        pafDict[hit.qName].alignments.append(hit)

    return pafDict

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
filteredVCF = call_MEI_candidate(VCF)
print('filteredVCF: ', filteredVCF, len(filteredVCF.variants))

## 3. Do MEI calling for candidate insertions
call_MEI(filteredVCF)
 

