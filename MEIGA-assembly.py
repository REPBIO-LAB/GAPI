## DEPENDENCIES ##
# External
import os
import sys
import argparse
import copy
import re

# Internal
import formats
import alignment
import unix
import sequences
import gRanges
import annotation

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

def call_MEI(vcf, consensus, reference, sourceDb, outDir):
    '''
    '''    
    ## 0. Create temporary folder
    tmpDir = outDir + '/tmp'
    unix.mkdir(tmpDir)

    ## 1. Write inserted sequences into fasta file
    fastaPath = tmpDir + '/insertions.fa'
    fasta = ins2fasta(vcf, tmpDir)
    fasta.write(fastaPath)

    ## 2. Create index for consensus sequences
    fileName = 'consensus'  
    consensusIndex = alignment.index_minimap2(consensus, fileName, tmpDir)

    ## 3. Align inserted sequences against consensus:
    PAF_path = alignment.alignment_minimap2(fastaPath, consensusIndex, 'hits2consensus', 1, tmpDir)
    #PAF_path = '/Users/brodriguez/Research/Projects/MEIGA/MEIGA/test/0.20.0/assembly_caller/tmp/hits2consensus.paf'
    PAF_consensus = formats.PAF()
    PAF_consensus.read(PAF_path)

    ## Align inserted sequences against the reference genome
    SAM_path = alignment.alignment_bwa(fastaPath, reference, 'hits2genome', 1, tmpDir)
    #SAM_path = '/Users/brodriguez/Research/Projects/MEIGA/MEIGA/test/0.20.0/assembly_caller/tmp/hits2genome.sam'
    PAF_path = alignment.sam2paf(SAM_path, 'hits2genome', tmpDir)
    PAF_genome = formats.PAF()
    PAF_genome.read(PAF_path)

    ## 4. Generate single PAF objects per inserted sequence:
    PAFs_consensus = group_alignments(PAF_consensus)
    PAFs_genome = group_alignments(PAF_genome)

    ## 5. Resolve structure for each insertion with matches on retrotransposon consensus sequences
    structures = {}

    for insId in PAFs_consensus:
        structures[insId] = MEI_structure(PAFs_consensus[insId], fasta.seqDict[insId])
        seqBeg, seqEnd = structures[insId]['CHAIN'].interval()

    ## 6. Resolve 3' partnered transductions
    structures = resolve_partnered_3prime(structures, fasta, reference, sourceDb, tmpDir)

    ## 6. Search for 5' partnered transductions
    structures = search4partnered_5prime(structures, fasta, reference, tmpDir)

    ## 7. Search for orphan transductions
    ## Remove resolved insertions
    for insId in structures:
        if structures[insId]['PASS']:
            del PAFs_genome[insId]

    ## Do orphan transduction search
    search4orphan(PAFs_genome, sourceDb, fasta) # TO FINISH LATER (Only two L1 orphan transductions so far..)

    ## 8. Generate output VCF containing MEI calls
    ## Create header for output dictionary
    outVCF = formats.VCF()
    outVCF.header = vcf.header

    ## Add MEI specific fields to the VCF header
    info2add = {'ITYPE': ['.', 'String', 'Type of insertion (solo, partnered or orphan)'], \
                '3PRIME': ['0', 'Flag', 'Partnered 3-prime transduction'], \
                '5PRIME': ['0', 'Flag', 'Partnered 5-prime transduction'], \
                'FAM': ['.', 'String', 'Repeat family'], \
                'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
                'TDCOORD': ['1', 'Integer', 'Transduced sequence coordinates'], \
                'TDLEN': ['1', 'Integer', 'Transduction length'], \
                'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
                }

    outVCF.header.info.update(info2add)

    ## Select INS corresponding to MEI calls and add update info field with MEI features
    for variant in vcf.variants:
        insId = variant.chrom + ':' + str(variant.pos) 
                            
        # Discard unresolved inserted sequences
        if (insId not in structures) or ((insId in structures) and (structures[insId]['PASS'] is False)):
            continue

        variant2add = copy.deepcopy(variant)
        variant2add.info.update(structures[insId])
        outVCF.add(variant2add)

    ## 9. Do cleanup
    #unix.rm([tmpDir])

    return outVCF

def search4orphan(hits, sourceDb, fasta):
    '''
    TO FINISH LATER
    '''
    
    # For each partnered event
    for insId in hits:
    
        # For each hit
        for hit in hits[insId].alignments:

            hit.tName = 'chr' + hit.tName

            #print('CANDIDATE: ', insId, hit.qLen, hit.qBeg, hit.qEnd, hit.alignmentPerc(), hit.MAPQ, fasta.seqDict[insId], sourceDb[hit.tName].collect_interval(hit.tBeg, hit.tEnd, 'ALL'))

            ## Filter out hits 
            if (hit.alignmentPerc() < 75) or (hit.MAPQ < 30) or (hit.tName not in sourceDb):
                continue            
            
            #if sourceDb[hit.tName].collect_interval(hit.tBeg, hit.tEnd, 'ALL'):
                #print('CALL: ', insId, hit.qLen, hit.qBeg, hit.qEnd, hit.alignmentPerc(), hit.MAPQ, fasta.seqDict[insId], sourceDb[hit.tName].collect_interval(hit.tBeg, hit.tEnd, 'ALL'))

def search4partnered_5prime(structures, fasta, reference, outDir):
    '''
    '''
    ## 1. Create Fasta with sequences to realign
    seq2realign = formats.FASTA()

    for insId in structures:

        # Discard if strand not determined
        if structures[insId]['STRAND'] is None:
            continue

        ## Extract unresolved 5' sequence if any
        qBeg, qEnd = structures[insId]['CHAIN'].interval()

        if structures[insId]['STRAND'] == '+':
            seq2realign.seqDict[insId] = fasta.seqDict[insId][:qBeg]

        else:
            seq2realign.seqDict[insId] = fasta.seqDict[insId][qEnd:]

    fastaPath = outDir + '/seq2realign.5prime.fasta'
    seq2realign.write(fastaPath)

    ## 2. Realign sequences on the reference with BWA-mem
    SAM_path = alignment.alignment_bwa(fastaPath, reference, 'hits2genome.5prime', 1, outDir)
    #SAM_path = '/Users/brodriguez/Research/Projects/MEIGA/MEIGA/test/0.20.0/assembly_caller/tmp/hits2genome.5prime.sam'
    PAF_path = alignment.sam2paf(SAM_path, 'hits2genome.5prime', outDir)
    
    PAF = formats.PAF()
    PAF.read(PAF_path)

    ## 3. Make 5' transduction calls
    # For each hit
    for hit in PAF.alignments:
        
        hit.tName = 'chr' + hit.tName
        iRef, coord = hit.qName.split(':')
        iBeg = int(coord) - 500
        iEnd = int(coord) + 500

        ## Filter out hits 
        if (hit.alignmentPerc() < 75) or (hit.MAPQ < 30) or (iRef == hit.tName and gRanges.overlap(iBeg, iEnd, hit.tBeg, hit.tEnd)[0]):
            continue

        ## Make call
        structures[hit.qName]['ITYPE'] = 'partnered'
        structures[hit.qName]['5PRIME'] = True
        structures[hit.qName]['TDCOORD'] = hit.tName + ':' + str(hit.tBeg) + '-' +  str(hit.tEnd)
        structures[hit.qName]['TDLEN'] = hit.tEnd - hit.tBeg
        
    return structures
        

def resolve_partnered_3prime(structures, fasta, reference, sourceDb, outDir):
    '''
    '''
    ## 1. Create Fasta with sequences to realign
    seq2realign = formats.FASTA()
    pattern = re.compile("Partnered_[0-9]+")
    partneredDict = {}

    for insId in structures:

        # Discard solo
        if structures[insId]['ITYPE'] != 'partnered':
            continue

        ## Initialize partnered dict for the ins
        partneredDict[insId] = {}
        partneredDict[insId]['NB_PARTNERED'] = 0
        partneredDict[insId]['NB_RESOLVED'] = 0

        # For each hit
        for hit in structures[insId]['CHAIN'].alignments:

            # Discard if not partnered
            if not pattern.match(hit.tName):
                continue

            # Add candidate partnered sequence to the fasta
            seqId = insId + '|' + hit.tName
            seq = fasta.seqDict[insId][hit.qBeg:hit.qEnd]
            seq2realign.seqDict[seqId] = seq

            ## Update partnered dictionary
            partneredDict[insId]['NB_PARTNERED'] += 1       
            partneredDict[insId][hit.tName] = hit


    fastaPath = outDir + '/seq2realign.3prime.fasta'
    seq2realign.write(fastaPath)

    ## 2. Realign sequences on the reference with BWA-mem
    SAM_path = alignment.alignment_bwa(fastaPath, reference, 'hits2genome.3prime', 1, outDir)
    #SAM_path = '/Users/brodriguez/Research/Projects/MEIGA/MEIGA/test/0.20.0/assembly_caller/tmp/hits2genome.3prime.sam'
    PAF_path = alignment.sam2paf(SAM_path, 'hits2genome.3prime', outDir)
    PAF = formats.PAF()
    PAF.read(PAF_path)

    ## 3. Add hit information to partnered transduction candidates
    hits = PAF.hits2dict()

    # For each partnered event
    for ID in hits:
    
        insId, tdId = ID.split('|')
        partneredDict[insId]['CYTOID'] = None

        # For each hit
        for hit in hits[ID]:

            hit.tName = 'chr' + hit.tName

            ## Check if it´s a partnered transduction from a known source element
            if (hit.tName in sourceDb) and sourceDb[hit.tName].collect_interval(hit.tBeg, hit.tEnd, 'ALL'):
                source = sourceDb[hit.tName].collect_interval(hit.tBeg, hit.tEnd, 'ALL')[0][0]
                partneredDict[insId]['CYTOID'] = source.optional['cytobandId']

            ## Filter out hits 
            if (hit.alignmentPerc() < 75 or hit.MAPQ < 30) and (partneredDict[insId]['CYTOID'] is None):
                continue
            
            ## Add hit information  
            partneredDict[insId]['NB_RESOLVED'] += 1
            partneredDict[insId][tdId].tName = hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd)
            
    ## 4. Filter partnered transduction candidates
    for insId in partneredDict:

        # a) Make transduction call
        if (partneredDict[insId]['NB_PARTNERED'] > 0) and (partneredDict[insId]['NB_PARTNERED'] == partneredDict[insId]['NB_RESOLVED']):

            tdIds = [key for key in partneredDict[insId].keys() if key not in ['NB_PARTNERED', 'NB_RESOLVED', 'CYTOID']]
            structures[insId]['TDCOORD'] = ','.join([partneredDict[insId][tdId].tName for tdId in tdIds])
            structures[insId]['CYTOID'] = partneredDict[insId]['CYTOID']
            structures[insId]['TDLEN'] = ','.join([str(partneredDict[insId][tdId].qEnd - partneredDict[insId][tdId].qBeg) for tdId in tdIds])

        # b) Filter out event
        else:
            structures[insId]['PASS'] = False
    
    return structures


def MEI_structure(PAF, insertSeq):
    '''
    '''
    structure = {}
    structure['LEN'] = len(insertSeq)

    ## 1. Chain alignments
    structure['CHAIN'] = PAF.chain(20, 50)

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
        dist2end = structure['LEN'] - monomersA[0].end

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'solo'
            structure['POLYA'] = 1
            structure['STRAND'] = '+'

    # b) Multiple polyA (Transduction candidate)
    # AAAAAAAAAAAAAA-------------AAAAAAAAAAAAAA
    elif (len(monomersA) > 1):
        dist2rt = monomersA[0].beg - rtEnd
        dist2end = structure['LEN'] - monomersA[-1].end

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'partnered'
            structure['3PRIME'] = True
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
            structure['ITYPE'] = 'partnered'
            structure['3PRIME'] = True
            structure['POLYT'] = len(monomersT)
            structure['STRAND'] = '-'
        
    ## 3.3 Determine if polyA or polyT found
    # a) PolyA found
    if (structure['POLYA'] != 0) and (structure['POLYT'] == 0):

        tail = 'polyA'
        structure['NBPOLY'] = structure['POLYA']
        monomers = monomersA

    # b) PolyT found
    elif (structure['POLYT'] != 0) and (structure['POLYA'] == 0):
        tail = 'polyT'
        structure['NBPOLY'] = structure['POLYT']
        monomers = monomersT

    # c) No tail or ambiguous
    else:
        tail = None
        structure['NBPOLY'] = 0

    ## 3.4 Determine candidate insertion type based on the number of polyA/T tails found
    hits2add = []

    # a) Solo
    if structure['NBPOLY'] == 1:
        fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomers[0].beg, monomers[0].end, None, 'PolyA/T', 0, 0, 0, 0, 0, 0]
        hit = formats.PAF_alignment(fields)
        hits2add.append(hit)

    # b) Partnered
    elif structure['NBPOLY'] > 1:

        ## First polyA/T
        fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomers[0].beg, monomers[0].end, None, 'PolyA/T', 0, 0, 0, 0, 0, 0]
        hit = formats.PAF_alignment(fields)
        hits2add.append(hit)

        ## Add partnered region/s plus polyAT/s
        counter = 1

        for monomer1, monomer2 in zip(monomers, monomers[1:]):

            # Partnered 
            fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomer1.end, monomer2.beg, None, 'Partnered_' + str(counter), 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            hits2add.append(hit)

            # Next polyA/T
            fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomer2.beg, monomer2.end, None, 'PolyA/T', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            hits2add.append(hit)

            counter += 1        
      
    ## 3.5 Add polyA/T plus transduced annotation to the chain
    if tail == 'polyA':
        structure['CHAIN'].alignments = structure['CHAIN'].alignments + hits2add

    elif tail == 'polyT':
        structure['CHAIN'].alignments = hits2add + structure['CHAIN'].alignments 

    # Compute % of inserted resolved
    structure['PERC-RESOLVED'] = structure['CHAIN'].perc_query_covered()

    ## 4. Apply filters 
    failed = []

    # 4.1 Percentage resolved filter
    if structure['PERC-RESOLVED'] < 60:
        failed.append('PERC-RESOLVED')

    # 4.2 Length filtering for solo insertions
    if structure['ITYPE'] == 'solo':

        if (structure['FAM'] == 'L1') and (structure['LEN'] > 6500):
            failed.append('LEN')
        
        elif (structure['FAM'] == 'Alu') and (structure['LEN'] > 400): 
            failed.append('LEN')

        elif (structure['FAM'] == 'SVA') and (structure['LEN'] > 4500): 
            failed.append('LEN')

    # a) Insertions passes all the filters
    if not failed:
        structure['PASS'] = True

    # b) At least one failed filter
    else:
        structure['PASS'] = False

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
 
        insId = variant.chrom + ':' + str(variant.pos) 
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
parser.add_argument('reference', help='Path to FASTA file containing the reference genome')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
vcf = args.vcf
consensus = args.consensus
reference = args.reference
fileName = args.fileName
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('vcf: ', vcf)
print('consensus: ', consensus)
print('reference: ', reference)
print('fileName: ', fileName)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 1. Read VCF
VCF = formats.VCF()
VCF.read(vcf)


## 2. Load source elements database
annotDir = '/Users/brodriguez/Research/Projects/MEIGA/MEIGA/databases/H.Sapiens/hg38/'
annotations = annotation.load_annotations(['TRANSDUCTIONS'], VCF.header.refLengths, annotDir, None, 1, outDir)

print('sourceDb: ', annotations['TRANSDUCTIONS'])

## 2. Filter VCF by selecting retrotransposition insertion candidates 
# (inserted sequences with polyA/T tails at their ends)
filteredVCF = call_MEI_candidate(VCF)

## 3. Do MEI calling for candidate insertions
outVCF = call_MEI(filteredVCF, consensus, reference, annotations['TRANSDUCTIONS'], outDir)

## 4. Write VCF containing MEI calls
IDS = ['VARTYPE', 'SVTYPE', 'SVLEN', 'CONTIG', 'CONTIG_COORD', 'CONTIG_STRAND', \
       'ITYPE', '3PRIME', '5PRIME', 'FAM', 'CYTOID', 'TDCOORD', 'TDLEN', 'STRAND']

outVCF.write(IDS, fileName, outDir)
