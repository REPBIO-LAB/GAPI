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

## FUNCTIONS ##
def retrotransposon_structure(FASTA_file, index, outDir):
    '''
    Determine if an input sequence correspond to a retrotransposition event (solo, partnered or orphan transduction from a known source element). 
    
    Infer the insertion size, structure, poly-A and target site duplication length

    Input:
        1. FASTA_file: Path to FASTA file containing the sequence
        2. index: Minimap2 index for consensus retrotransposon sequences database
        3. outDir: Output directory
        
    Output:
        1. insType: Insertion type (solo, nested, orphan, partnered or None)
        2. family: List of retrotransposon families
        3. srcId: List of source element ids
        4. strand. Insertion strand (+ or -)
        5. polyA: boolean specifying if polyA/T sequence was found
        6. structure: dictionary containing insertion structure information
        7. mechanism: TPRT, EI or unknown
    '''         
    ## 0. Create logs directory ##
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align the sequence into the retrotransposon sequences database ##
    PAF_file = alignment.alignment_minimap2(FASTA_file, index, 1, outDir)

    ## 2. Read PAF alignments ##
    PAF = formats.PAF()
    PAF.read(PAF_file)

    # Exit function if no hit on the retrotransposons database
    if not PAF.lines:
        percResolved, insType, family, srcId, strand, polyA, structure, mechanism = [0, None, [], [], None, False, {}, 'unknown']
        return percResolved, insType, family, srcId, strand, polyA, structure, mechanism

    ## 3. Chain complementary alignments ##
    chain = PAF.chain(100, 20)

    ## 4. Infer insertion features ##
    ## Retrieve inserted seq
    FASTA = formats.FASTA()
    FASTA.read(FASTA_file)
    sequence = list(FASTA.seqDict.values())[0]

    ## 4.1 Insertion type
    insType, family, srcId = insertion_type(chain)

    ## 4.2 Insertion strand
    strand, polyA = infer_strand(insType, sequence, chain)

    ## 4.3 Sequence structure 
    structure = infer_structure(insType, chain, strand)

    ## 4.4 Insertion mechanism (TPRT or EI)
    mechanism = infer_integration_mechanism(chain, structure['truncation3len'], polyA)

    ## 4.5 Target site duplication (TO DO LATER...)
    #search4tsd()
    
    return chain.perc_query_covered(), insType, family, srcId, strand, polyA, structure, mechanism
    

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


def infer_structure(insType, chain, strand):
    '''
    Infer inserted sequence structural features
    
    Input:
        1. insType: Insertion type (solo, nested, orphan, partnered or None)
        2. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        3. strand: Insertion strand (+ or -)

    Output:
        1. structure: dictionary containing insertion structure information
    ''' 
    ### Initialize dictionary
    structure = {}

    for feature in ['retroCoord', 'retroLen', 'isFull', 'truncation5len', 'truncation3len', 'transductionCoord', 'transductionLen', 'inversionLen']:
        structure[feature] = None

    ### 1. Compute the length of each type of sequence composing the insertion
    # 1.1 Retroelement length
    if insType in ['solo', 'partnered']:

        ## Pick only those hits over retrotransposon consensus sequence
        retroHits = [hit for hit in chain.alignments if 'consensus' in hit.tName]

        ## Determine piece of consensus sequence that has been integrated
        ref = retroHits[0].tName.split('|')[1]
        retroBeg = min([hit.tBeg for hit in retroHits])
        retroEnd = max([hit.tEnd for hit in retroHits])
        structure['retroCoord'] = str(ref) + ':' + str(retroBeg) + '-' + str(retroEnd)

        ## Compute length
        structure['retroLen'] = retroEnd - retroBeg

        ## Assess if full length retrotransposon insertion
        consensusLen = retroHits[0].tLen 
        percConsensus = float(structure['retroLen']) / consensusLen * 100
        structure['isFull'] = True if percConsensus >= 95 else False

        ## Compute truncation length at both ends
        structure['truncation5len'] = retroBeg   
        structure['truncation3len'] = consensusLen - retroEnd

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
        structure['transductionCoord'] = (ref, transductionBeg, transductionEnd)

        ## Compute length
        structure['transductionLen'] = transductionEnd - transductionBeg

    # 1.3 Inversion length
    inversionHits = [hit for hit in chain.alignments if ((hit.strand != 'None') and (hit.strand != strand))]

    # a) 5' inversion
    if inversionHits:
        structure['inversionLen'] = sum([hit.tEnd - hit.tBeg for hit in inversionHits])

    # b) No inversion
    else:
        structure['inversionLen'] = 0   
     
    return structure


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
