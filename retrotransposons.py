'''
Module 'retrotransposons' - Contains functions for the identification and characterization of retrotransposon sequences
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
import log
import unix
import formats
import alignment
import sequences
import annotation

## FUNCTIONS ##
def retrotransposon_structure(FASTA_file, index, outDir):
    '''    
    Infer the insertion size, structure, poly-A, target site duplication length and other insertion structural features

    Input:
        1. FASTA_file: Path to FASTA file containing the sequence
        2. index: Minimap2 index for consensus retrotransposon sequences database
        3. outDir: Output directory
        
    Output:
        1. structure: dictionary containing insertion structure information
    '''         
    structure = {}

    ## 0. Create logs directory ##
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align the sequence into the retrotransposon sequences database ##
    PAF_file = alignment.alignment_minimap2(FASTA_file, index, 1, outDir)

    ## 2. Read PAF alignments ##
    PAF = formats.PAF()
    PAF.read(PAF_file)

    # Exit function if no hit on the retrotransposons database
    if not PAF.alignments:
        return structure

    ## 3. Chain complementary alignments ##
    chain = PAF.chain(100, 20)

    ## 4. Infer insertion features ##
    ## Retrieve inserted seq
    FASTA = formats.FASTA()
    FASTA.read(FASTA_file)
    sequence = list(FASTA.seqDict.values())[0]

    ## 4.1 Insertion type
    structure['INS_TYPE'], structure['FAMILY'], structure['CYTOBAND'] = insertion_type(chain)

    ## 4.2 Insertion strand
    structure['STRAND'], structure['POLYA'] = infer_strand(structure['INS_TYPE'], sequence, chain)

    ## 4.3 Sequence lengths 
    lengths = infer_lengths(structure['INS_TYPE'], chain, structure['STRAND'])
    structure.update(lengths)

    ## 4.4 Insertion mechanism (TPRT or EI)
    structure['MECHANISM'] = infer_integration_mechanism(chain, structure['TRUNCATION_3_LEN'], structure['POLYA'])

    ## 4.5 Target site duplication (TO DO LATER...)
    #search4tsd()
    
    ## 4.6 Percentage resolved
    structure['PERC_RESOLVED'] = chain.perc_query_covered()

    return structure
    

def insertion_type(chain):
    '''
    Scan alignments chain to determine the type of insertion (solo, transduction...)

    Input:
        1. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        
    Output:
        1. insType: Insertion type (solo, nested, orphan, partnered or None)
        2. family: List of retrotransposon families
        3. srcId: List of source element ids
    ''' 
    ## Make list containing all the different templates the sequence aligns into
    templateTypes = list(set([alignment.tName.split("|")[0] for alignment in chain.alignments]))
    nbTemplateTypes = len(templateTypes)

    ## Make list containing the id for the source element transduced regions the sequence aligns into
    sourceElements = list(set([alignment.tName.split("|")[2] for alignment in chain.alignments if ('transduced' in alignment.tName)]))
    nbSource = len(sourceElements)

    ## Make list containing the families the sequence aligns into 
    families = list(set([alignment.tName.split("|")[1] for alignment in chain.alignments if ('consensus' in alignment.tName)]))
    nbFamilies = len(families)

    ## a) Solo insertion
    # //////RT//////     
    if (nbTemplateTypes == 1) and ('consensus' in templateTypes) and (nbFamilies == 1):
        insType = 'solo'
        family = families
        srcId = []

    ## b) Nested insertion (Insertion composed by multiple retrotransposons from different families)
    # //////RT_1//////\\\\\\\\RT_2\\\\\\\     
    elif (nbTemplateTypes == 1) and ('consensus' in templateTypes):
        insType = 'nested'
        family = families
        srcId = []
        
    ## c) Orphan (inserted sequence only matching one transduced region)
    # SOURCE//////TD//////    
    elif (nbTemplateTypes == 1) and ('transduced' in templateTypes) and (nbSource == 1):
        insType = 'orphan'
        family = []
        srcId = sourceElements

    ## d) Partnered (inserted sequence matching consensus and one transduced sequence)
    # SOURCE >>>>>>L1>>>>>>//////TD///////    
    elif (nbTemplateTypes == 2) and (set(['consensus', 'transduced']).issubset(templateTypes)) and (nbFamilies == 1) and (nbSource == 1):
        insType = 'partnered'
        family = families
        srcId = sourceElements

    ## e) Unknown insertion type
    else:
        insType = None
        family = []
        srcId = [] 

    return insType, family, srcId
    

def infer_strand(insType, sequence, chain):
    '''
    Infer insertion strand based on two criteria:
        1) Alignment orientation for the insert 3' end over the template sequence
        2) Location of polyA/T tail at sequence ends

    Input:
        1. insType: Insertion type (solo, nested, orphan, partnered or None)
        2. sequence: consensus inserted sequence
        3. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions

    Output:
        1. strand: Insertion strand (+, - or None) 
        2. polyA: boolean specifying if polyA/T sequence was found
    '''

    ## 1. Strand based on polyA/T presence 
    strandPolyA, polyA = infer_strand_polyA(sequence, chain)

    ## 2. Strand based on alignment orientation
    strandAlignment = infer_strand_alignment(insType, chain)

    ## 3. Define consensus strand
    # a) PolyA/T has preference over the alignment orientation
    if strandPolyA is not None:
        strand = strandPolyA
    
    # b) Strand based on the alignment orientation
    elif strandAlignment is not None:
        strand = strandAlignment
    
    # c) Unknown strand
    else:
        strand = None
    
    return strand, polyA
    

def infer_strand_polyA(sequence, chain):
    '''
    Infer insertion strand based on two criteria:
        1) Location of polyA/T tail at sequence ends
        2) Alignment strand for the insert 3' end over the template sequence

    Input: 
        1. sequence: consensus inserted sequence
        2. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions

    Output:
        1. strand: Insertion strand (+, - or None) 
        2. polyA: boolean specifying if polyA/T sequence was found
    '''
    
    ### Set up configuration parameters
    windowSize = 8
    maxWindowDist = 2
    minMonomerSize = 10
    minPurity = 80 

    maxDist2Ends = 10 
    minInternalMonomerSize = 20

    ## 1. Search for polyA at the insert 3' end ##
    # 1.1 Extract unaligned 3' end of the inserted sequence
    lastHit = chain.alignments[-1]
    targetSeq = sequence[lastHit.qEnd:]

    # 1.2 Search for poly(A) monomers on the 3' end 
    targetMonomer = 'A'
    monomers3end = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    # 1.3 Filter internal monomers
    monomers3end = sequences.filter_internal_monomers(monomers3end, targetSeq, maxDist2Ends, minInternalMonomerSize)

    ## 2. Search for polyT at the insert 5' end ##
    # 2.1 Extract unaligned 5' end of the inserted sequence
    firstHit = chain.alignments[0]
    targetSeq = sequence[:firstHit.qBeg]

    # 2.2 Search for poly(T) monomers on the 5' end 
    targetMonomer = 'T'
    monomers5end = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    # 2.3 Filter internal monomers
    monomers5end = sequences.filter_internal_monomers(monomers5end, targetSeq, maxDist2Ends, minInternalMonomerSize)

    ## 3. Determine strand ##
    # 3.1 Compute 3' monomers accumulative length
    monomers3endLengths = [monomer.length() for monomer in monomers3end]
    accumulative3prime = sum(monomers3endLengths)
     
    # 3.2 Compute 5' monomers accumulative length
    monomers5endLengths = [monomer.length() for monomer in monomers5end]
    accumulative5prime = sum(monomers5endLengths)

    # 3.3 Determine if polyA/T at 5' or 3' end (indicative of strand orientation) 
    # a) Unknown strand if:
    # - No monomer found in any end OR
    # - Ambiguity as 3' and 5' monomers with equal size
    if ((accumulative3prime == 0) and (accumulative5prime == 0)) or (accumulative3prime == accumulative5prime) :
        monomers = None
        strand = None
        polyA = False

    # b) Positive strand
    elif accumulative3prime > accumulative5prime:
        monomers = monomers3end
        strand = '+'
        polyA = True

    # c) Negative strand
    else:
        monomers = monomers5end
        strand = '-'
        polyA = True

    ## 4. Convert monomer coordinates to inserted sequence space ##
    ## a) + strand
    # -----insert------**unaligned**
    if strand == '+':

        for monomer in monomers:
            monomer.beg = monomer.beg + lastHit.qEnd
            monomer.end = monomer.end + lastHit.qEnd 

    ## b) - strand
    # NOT NEEDED as the unaligned sequence correspond to the leftmost end of the insert 
    # **unaligned**-----insert------ 

    ## 5. Add polyA/T to the chain of alignments
    firstAlignment = chain.alignments[0]

    ## a) + strand
    if strand == '+':

        # For each monomer
        for monomer in monomers:

            # Create PAF line containing poly(A) info
            fields = [firstAlignment.qName, firstAlignment.qLen, monomer.beg, monomer.end, None, 'poly(A/T)', 0, 0, 0, 0, 0, 0]            
            alignment = formats.PAF_line(fields)

            # Add to the end of the chain
            chain.alignments.append(alignment) 

    ## b) - strand
    elif strand == '-':

        # For each monomer
        for monomer in monomers[::-1]:

            # Create PAF line containing poly(T) info
            fields = [firstAlignment.qName, firstAlignment.qLen, monomer.beg, monomer.end, None, 'poly(A/T)', 0, 0, 0, 0, 0, 0]
            alignment = formats.PAF_line(fields)

            # Add to the begin of the chain
            chain.alignments.insert(0, alignment) 

    return strand, polyA


def infer_strand_alignment(insType, chain):
    '''
    Determine insertion strand (+ or -). The insertion strand can be infered from the alignment
    strand for the insert 3' end over the template sequence
    
    Input:
        1. insType: Insertion type (solo, nested, orphan, partnered or None)
        2. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        
    Output:
        1. strand: Insertion strand (+, - or None) 
    ''' 
    ## a) Solo or orphan transduction    
    if insType in ['solo', 'orphan']:

        ## Sort hits from 5' to 3'
        sortedHits = sorted(chain.alignments, key=lambda alignment: alignment.tEnd, reverse=False)

        ## Select 3' end hit as its strand will be == insertion strand
        strand = sortedHits[-1].strand
        
    ## b) Partnered transduction
    elif insType == 'partnered':

        ## Select hits over the transduced region
        transductionHits = [alignment for alignment in chain.alignments if 'transduced' in alignment.tName]

        ## Sort hits from 5' to 3'
        sortedHits = sorted(transductionHits, key=lambda alignment: alignment.tEnd, reverse=False)

        ## Select 3' end hit as its strand will be == insertion strand
        strand = sortedHits[-1].strand

    ## c) Other
    else:
        strand = None

    return strand


def infer_lengths(insType, chain, strand):
    '''
    Determine the length of each type of sequence composing the insertion (transduction, retrotransposon, insertion, inversion)
    
    Input:
        1. insType: Insertion type (solo, nested, orphan, partnered or None)
        2. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        3. strand: Insertion strand (+ or -)

    Output:
        1. lengths: dictionary containing length information
    ''' 
    ### Initialize dictionary
    lengths = {}

    for feature in ['RETRO_COORD', 'RETRO_LEN', 'IS_FULL', 'TRUNCATION_5_LEN', 'TRUNCATION_3_LEN', 'TRANSDUCTION_COORD', 'TRANSDUCTION_LEN', 'INVERSION_LEN']:
        lengths[feature] = None

    ### 1. Compute the length of each type of sequence composing the insertion
    # 1.1 Retroelement length
    if insType in ['solo', 'partnered']:

        ## Pick only those hits over retrotransposon consensus sequence
        retroHits = [hit for hit in chain.alignments if 'consensus' in hit.tName]

        ## Determine piece of consensus sequence that has been integrated
        ref = retroHits[0].tName.split('|')[1]
        retroBeg = min([hit.tBeg for hit in retroHits])
        retroEnd = max([hit.tEnd for hit in retroHits])
        lengths['RETRO_COORD'] = str(ref) + ':' + str(retroBeg) + '-' + str(retroEnd)

        ## Compute length
        lengths['RETRO_LEN'] = retroEnd - retroBeg

        ## Assess if full length retrotransposon insertion
        consensusLen = retroHits[0].tLen 
        percConsensus = float(lengths['RETRO_LEN']) / consensusLen * 100
        lengths['IS_FULL'] = True if percConsensus >= 95 else False

        ## Compute truncation length at both ends
        lengths['TRUNCATION_5_LEN'] = retroBeg   
        lengths['TRUNCATION_3_LEN'] = consensusLen - retroEnd

    # 1.2 Transduction length
    if insType in ['partnered', 'orphan']:

        ## Pick only those hits over transduced region
        transductionHits = [hit for hit in chain.alignments if 'transduced' in hit.tName]
 
        ## Compute offset to translate from interval to genomic coordinates
        interval = transductionHits[0].tName.split('|')[3]
        ref, coord = interval.split(':')
        offset = int(coord.split('-')[0])
        
        ## Determine piece of transduced area that has been integrated
        transductionBeg = min([hit.tBeg for hit in transductionHits]) + offset
        transductionEnd = max([hit.tEnd for hit in transductionHits]) + offset
        lengths['TRANSDUCTION_COORD'] = (ref, transductionBeg, transductionEnd)

        ## Compute length
        lengths['TRANSDUCTION_LEN'] = transductionEnd - transductionBeg

    # 1.3 Inversion length
    inversionHits = [hit for hit in chain.alignments if ((hit.strand != 'None') and (hit.strand != strand))]

    # a) 5' inversion
    if inversionHits:
        lengths['INVERSION_LEN'] = sum([hit.tEnd - hit.tBeg for hit in inversionHits])

    # b) No inversion
    else:
        lengths['INVERSION_LEN'] = 0   
     
    return lengths


def infer_integration_mechanism(chain, truncation3len, polyA):
    '''
    Determine the mechanism of integration (TPRT: Target Primed Reversed Transcription; EI: Endonuclease Independent)
    
    Input:
        1. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        2. truncation3len: number of base pairs the inserted sequence has been truncated at its 3'
        3. polyA: boolean specifying if polyA/T sequence was found

    Output:
        1. mechanism: TPRT, EI or unknown
    '''  
    ## A) TPRT hallmarks: 
    # - 3' truncation <= 100bp OR None (if orphan transduction)
    # - polyA 
    # - Incorporate TSD later as an additional evidence...
    if ((truncation3len is None) or (truncation3len <= 100)) and polyA:
        mechanism = 'TPRT'

    ## B) EI hallmarks:
    # - % resolved > 95%
    # - 3' truncation > 100bp 
    # - no polyA
    elif (chain.perc_query_covered() >= 95) and ((truncation3len is not None) and (truncation3len > 100)) and not polyA:
        mechanism = 'EI'

    ## C) Unknown mechanism:
    else:
        mechanism = 'unknown'

    return mechanism


def is_interspersed_ins(sequence, PAF, repeatsDb, transducedDb):
    '''
    Determine if input sequence corresponds to a interspersed insertion

    Input:
        1. sequence: Input sequence
        2. PAF: PAF object containing input sequence alignments on the reference genome
        3. repeatsDb: bin database containing annotated repeats in the reference. None if not available
        4. transducedDb: bin database containing regions transduced by source elements. None if not available

    Output:
        1. INTERSPERSED: Boolean specifying if inserted sequence corresponds to an intersersed repeat (True) or not (False)
        2. INS_features: dictionary containing interspersed repeat insertion features
        3. chain: alignments chain on the reference. None if sequence does not align on the reference
    '''
    INS_features = {}

    ## 0. Sequence does not align on the reference ##
    if not PAF.alignments:
        INTERSPERSED = False
        INS_features['INS_TYPE'] = 'unknown'
        INS_features['PERC_RESOLVED'] = 0
        return INTERSPERSED, INS_features, None

    ## 1. Create chain of alignments ##
    chain = PAF.chain(300, 20)

    ## 2. Annotate each hit in the chain ##
    repeatMatch = False
    transducedMatch = False

    ## For each hit
    for hit in chain.alignments:

        hit.annot = {}

        ## 2.1 Does the hit corresponds to a polyA/T? 
        polyA = is_polyA(sequence, 80)

        # Skip next annotation steps if hit corresponds to a polyA
        if polyA:
            continue

        ## 2.2 Hit matches an annotated repeat 
        overlaps = annotation.annotate_interval(hit.tName, hit.tBeg, hit.tEnd, repeatsDb)

        if overlaps:
            repeatMatch = True
            hit.annot['REPEAT'] = overlaps[0][0] # Select repeat with longest overlap

        ## 2.3 Hit matches a transduced region
        overlaps = annotation.annotate_interval(hit.tName, hit.tBeg, hit.tEnd, transducedDb)

        if overlaps:
            transducedMatch = True
            hit.annot['SOURCE_ELEMENT'] = overlaps[0][0] # Select trandsduced with longest overlap

    ## 3. Make list of distinct features matching the inserted sequence ##
    features = {}
    features['SOURCE_ELEMENT'] = []   
    features['REPEATS'] = {}
    features['REPEATS']['FAMILIES'] = []
    features['REPEATS']['SUBFAMILIES'] = []

    for hit in chain.alignments:

        ## Hit overlaps a repeat
        if 'REPEAT' in hit.annot:
            family = hit.annot['REPEAT'].optional['family']
            subfamily = hit.annot['REPEAT'].optional['subfamily']

            if family not in features['REPEATS']['FAMILIES']:
                features['REPEATS']['FAMILIES'].append(family)

            if subfamily not in features['REPEATS']['SUBFAMILIES']:
                features['REPEATS']['SUBFAMILIES'].append(subfamily)

        ## Hit overlaps a source element
        if 'SOURCE_ELEMENT' in hit.annot:
            cytobandId = hit.annot['SOURCE_ELEMENT'].optional['cytobandId']
            family = hit.annot['SOURCE_ELEMENT'].optional['family']

            if cytobandId not in features['SOURCE_ELEMENT']:
                features['SOURCE_ELEMENT'].append(cytobandId)

            if family not in features['REPEATS']['FAMILIES']:
                features['REPEATS']['FAMILIES'].append(family)

    ## 4. Determine insertion type based on hits annotation ##
    # A) Partnered transduction
    if repeatMatch and transducedMatch:

        INTERSPERSED = True
        INS_features['INS_TYPE'] = 'partnered'

        ## Repeat info
        INS_features['FAMILY'] = features['REPEATS']['FAMILIES'] 
        INS_features['SUBFAMILY'] = features['REPEATS']['SUBFAMILIES']

        ## Transduction info
        INS_features['CYTOBAND'] = features['SOURCE_ELEMENT']

    # B) Orphan
    elif transducedMatch:
        INTERSPERSED = True
        INS_features['INS_TYPE'] = 'orphan'

        INS_features['FAMILY'] = features['REPEATS']['FAMILIES'] 
        INS_features['CYTOBAND'] = features['SOURCE_ELEMENT']

    # C) Solo
    elif repeatMatch:
        INTERSPERSED = True
        INS_features['INS_TYPE'] = 'solo'

        INS_features['FAMILY'] = features['REPEATS']['FAMILIES'] 
        INS_features['SUBFAMILY'] = features['REPEATS']['SUBFAMILIES']

    # D) PolyA/T tail
    elif polyA:
        INTERSPERSED = True
        INS_features['INS_TYPE'] = 'poly(A/T)'

    # E) Unknown 
    else:       
        INTERSPERSED = False
        INS_features['PERC_RESOLVED'] = 0
        INS_features['INS_TYPE'] = 'unknown'

    return INTERSPERSED, INS_features, chain


def is_polyA(sequence, minPerc):
    '''
    Determine if input sequence corresponds to a polyA/T tail or not

    Input:
        1. sequence: Input sequence
        2. minPerc: Minimum percentage of base pairs corresponding to A or T to call a polyA/T tail

    Output:
        1. polyA: boolean specifying if input sequence corresponds to a polyA/T tail
    '''
    ## 1. Assess input sequence base composition
    baseCounts, basePercs = sequences.baseComposition(sequence)

    ## 2. Determine if polyA/T
    # a) PolyA/T
    if (basePercs['A'] >= minPerc) or (basePercs['T'] >= minPerc):
        polyA = True
    
    # b) Not PolyA/T
    else:
        polyA = False
    
    return polyA


