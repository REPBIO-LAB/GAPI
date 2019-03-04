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
def is_retrotransposition(FASTA_file, index, outDir):
    '''
    Determine if an input sequence correspond to a retrotransposition event (solo, partnered or orphan transduction from a known source element). 
    
    Infer the insertion size, structure, poly-A and target site duplication length

    Input:
        1. FASTA_file: Path to FASTA file containing the sequence
        2. index: Minimap2 index for consensus retrotransposon sequences database
        3. outDir: Output directory
        
    Output:
    ''' 
    ## 0. Create logs directory ##
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align the sequence into the retrotransposon sequences database ##
    PAF_file = outDir + '/alignments.paf'
    err = open(logDir + '/align.err', 'w') 
    command = 'minimap2 ' + index + ' ' + FASTA_file + ' > ' + PAF_file
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN-INSERT'
        msg = 'Insert alignment failed' 
        log.step(step, msg)

    ## 2. Read PAF alignments ##
    PAF = formats.PAF()
    PAF.read(PAF_file)

    # Exit function if no hit on the retrotransposons database
    if not PAF.lines:
        insType, family, srcId, status, percCovered, strand, hits = [None, None, None, None, None, None, None]
        return insType, family, srcId, status, percCovered, strand, hits

    ## 3. Chain complementary alignments ##
    chain = PAF.chain()
    
    ## 4. Infer insertion features ##
    ## 4.1 Insertion type
    insType, family, srcId = insertion_type(chain)

    ## 4.2 Insertion strand
    strand = infer_strand(insType, chain)
    
    ## 4.3 Identify PolyA tail
    ## Retrieve inserted seq
    FASTA = formats.FASTA()
    FASTA.read(FASTA_file)
    sequence = list(FASTA.fastaDict.values())[0]

    ## Search for polyA at 3' insert end
    search4polyA(sequence, chain, strand)

    ## 4.4 Sequence structure (TO DO LATER...)
    # infer_structure()
    
    ## 4.5 Target site duplication (TO DO LATER...)
    #search4tsd()

    ## 5. Determine status accoding to how much % of its sequence has been resolved ## 
    percCovered = chain.perc_query_covered()

    if (percCovered >= 70):
        status = 'resolved'

    elif (percCovered >= 40):
        status = 'partially_resolved'        

    else:
        status = 'unresolved'      

    ## Make hits string
    hits = ['_'.join([str(hit.qBeg), str(hit.qEnd), hit.strand, hit.tName]) for hit in chain.alignments]
    hits = ';'.join(hits)

    return insType, family, srcId, status, percCovered, strand, hits

def insertion_type(chain):
    '''
    Scan alignments chain to determine the type of insertion (solo, transduction...)

    Input:
        1. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        
    Output:
        1. insType: Insertion type 
    ''' 
    
    ## Make list containing all the different templates the sequence aligns into
    templateTypes = list(set([alignment.tName.split("|")[0] for alignment in chain.alignments]))
    nbTemplateTypes = len(templateTypes)

    ## Make list containing the id for the source element transduced regions the sequence aligns into
    sourceElements = list(set([alignment.tName.split("|")[1] for alignment in chain.alignments if ('Transduced5prime' in alignment.tName) or ('Transduced3prime' in alignment.tName)]))

    ## Make list containing the families the sequence aligns into 
    families = list(set([alignment.tName.split("|")[2] for alignment in chain.alignments if ('consensus' in alignment.tName) or ('consensus' in alignment.tName)]))
    nbFamilies = len(families)

    ## a) Unknown insertion type
    if not templateTypes:
        insType = None
        family = None
        srcId = None

    ## b) Solo insertion
    # //////RT//////     
    elif (nbTemplateTypes == 1) and ('consensus' in templateTypes):
        insType = 'Solo'
        family = ';'.join(families)
        srcId = None

    ## c) Nested insertion (Insertion composed by multiple retrotransposons from different families)
    # //////RT_1//////\\\\\\\\RT_2\\\\\\\     
    elif (nbTemplateTypes == 1) and ('consensus' in templateTypes) and (nbFamilies > 1):
        insType = 'Nested'
        family = ';'.join(families)
        srcId = None

    ## d) Orphan 5' 
    # //////TD//////SOURCE    
    elif (nbTemplateTypes == 1) and ('Transduced5prime' in templateTypes):
        insType = 'Orphan5prime'
        family = None
        srcId = ';'.join(sourceElements)
        
    ## e) Orphan 3' 
    # SOURCE//////TD//////    
    elif (nbTemplateTypes == 1) and ('Transduced3prime' in templateTypes):
        insType = 'Orphan3prime'
        family = None
        srcId = ';'.join(sourceElements)

    ## f) Partnered 5' 
    # SOURCE //////TD///////>>>>>>L1/SVA>>>>>>    
    elif (nbTemplateTypes == 2) and (set(['consensus', 'Transduced5prime']).issubset(templateTypes)):
        insType = 'Partnered5prime'
        family = ';'.join(families)
        srcId = ';'.join(sourceElements)

    ## g) Partnered 3'
    # SOURCE >>>>>>L1>>>>>>//////TD///////    
    elif (nbTemplateTypes == 2) and set(['consensus', 'Transduced3prime']).issubset(templateTypes):
        insType = 'Partnered3prime'
        family = ';'.join(families)
        srcId = ';'.join(sourceElements)

    ## h) Partnered transduction containing both 5' and 3' transduction
    else:
        insType = 'Partnered5,3prime'
        family = ';'.join(families)
        srcId = ';'.join(sourceElements)

    return insType, family, srcId


def infer_strand(insType, chain):
    '''
    Determine insertion strand (+ or -). The insertion strand can be infered from the alignment
    strand for the insert 3' end over the template sequence
    
    Input:
        1. insType: Insertion type (solo, transduction, ...)
        2. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        
    Output:
        1. strand. Insertion strand (+ or -)
    ''' 

    ## a) Solo or orphan transduction    
    if insType in ['Solo', 'Orphan5prime', 'Orphan3prime']:

        ## Sort hits from 5' to 3'
        sortedHits = sorted(chain.alignments, key=lambda alignment: alignment.tEnd, reverse=False)

        ## Select 3' end hit as its strand will be == insertion strand
        strand = sortedHits[-1].strand
        
    ## b) Partnered
    elif insType in ['Partnered5prime', 'Partnered3prime']:

        ## Select hits over the transduced region
        transductionHits = [alignment for alignment in chain.alignments if 'Transduced' in alignment.tName]

        ## Sort hits from 5' to 3'
        sortedHits = sorted(transductionHits, key=lambda alignment: alignment.tEnd, reverse=False)

        ## Select 3' end hit as its strand will be == insertion strand
        strand = sortedHits[-1].strand

    ## c) 5' and 3' partnered transductions
    elif insType == 'Partnered5,3prime':

        ## Select hits over the transduced region
        transductionHits = [alignment for alignment in chain.alignments if 'Transduced3prime' in alignment.tName]

        ## Sort hits from 5' to 3'
        sortedHits = sorted(transductionHits, key=lambda alignment: alignment.tEnd, reverse=False)

        ## Select 3' end hit as its strand will be == insertion strand
        strand = sortedHits[-1].strand

    ## d) Nested
    else:
        strand = None

    return strand


def search4polyA(sequence, chain, strand):
    '''
    Search for polyA/T sequence at 3' insert end

    Input:
        1. sequence: Inserted fragment 
        2. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        3. strand: Insertion strand (+ or -)
        
    Output:
        Update alignments chain object including poly(A/T) hits
    '''
    ### Set up configuration parameters
    windowSize = 5
    maxWindowDist = 2
    minMonomerSize = 8
    minPurity = 70 

    maxDist2Ends = 10 
    minInternalMonomerSize = 20

    ### Search for poly A/T tracts in the 3' ends of the insert 
    ## A) + strand
    if strand == '+':

        ## 1. Extract unaligned 3' end of the inserted sequence
        lastHit = chain.alignments[-1]
        targetSeq = sequence[lastHit.qEnd:]

        ## 2. Search for poly(A) monomers on the 3' end 
        targetMonomer = 'A'
        monomers = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

        ## 3. Request a larger size for internal monomers
        ## An internal monomer is located at more than Xbp from the inserted sequence ends
        filteredMonomers = []

        for monomer in monomers:
            dist2Beg = monomer.beg
            dist2End = len(targetSeq) - monomer.end

            ## a) Select external monomers located at less than or equal to X bp from insert ends 
            if (dist2Beg <= maxDist2Ends) or (dist2End <= maxDist2Ends):
                filteredMonomers.append(monomer)
                            
            ## b) Select internal monomers longer or equal than X bp
            elif (monomer.length() >= minInternalMonomerSize):
                filteredMonomers.append(monomer)
            
            ## c) Discard internal monomers smaller than X bp

        ## 4. Convert monomer coordinates to inserted sequence space
        for monomer in filteredMonomers:
            monomer.beg = monomer.beg + lastHit.qEnd
            monomer.end = monomer.end + lastHit.qEnd

        ## 5. Add poly(A) info to the chain of alignments         
        nbMonomers = len(filteredMonomers)
        
        ## a) Single poly(A) 
        if nbMonomers == 1:

            monomer = filteredMonomers[0]

            ## Create PAF line containing poly(A) info
            firstAlignment = chain.alignments[0]
            fields = [firstAlignment.qName, firstAlignment.qLen, monomer.beg, monomer.end, None, 'poly(A/T)', 0, 0, 0, 0, 0, 0]
            alignment = formats.PAF_line(fields)

            chain.alignments.append(alignment) ## Add to the chain
        
        ## b) Two poly(A) tracts. The most likely explanation of this structure is a short transduction
        # First poly(A) corresponds to the source L1 element
        # Second poly(A) corresponds to the transcript polyA tail
        # Sequence in between corresponds to a short transduction (TD)
        # AAAAAAAAAA-----TD-----AAAAAAAAAA
        elif nbMonomers == 2:
            firstMonomer, secondMonomer = filteredMonomers

            ## Create PAF line containing TD info
            firstAlignment = chain.alignments[0]
            fields = [firstAlignment.qName, firstAlignment.qLen, firstMonomer.beg, secondMonomer.beg, None, 'Transduced3prime_UNK_UNK', 0, 0, 0, 0, 0, 0]
            alignment = formats.PAF_line(fields)
            chain.alignments.append(alignment) ## Add to the chain

            ## Create PAF line containing poly(T) info
            fields = [firstAlignment.qName, firstAlignment.qLen, secondMonomer.beg, secondMonomer.end, None, 'poly(A/T)', 0, 0, 0, 0, 0, 0]
            alignment = formats.PAF_line(fields)
            chain.alignments.append(alignment) ## Add to the chain

    ## B) - strand
    else:
        
        ## 1. Extract unaligned 3' end of the inserted sequence
        firstHit = chain.alignments[0]
        targetSeq = sequence[:firstHit.qBeg]

        ## 2. Search for poly(T) monomers on the 3' end
        targetMonomer = 'T'
        monomers = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

        ## 3. Request a larger size for internal monomers
        ## An internal monomer is located at more than Xbp from the inserted sequence ends
        filteredMonomers = []

        for monomer in monomers:
            dist2Beg = monomer.beg
            dist2End = len(targetSeq) - monomer.end

            ## a) Select external monomers located at less than or equal to X bp from insert ends 
            if (dist2Beg <= maxDist2Ends) or (dist2End <= maxDist2Ends):
                filteredMonomers.append(monomer)
                            
            ## b) Select internal monomers longer or equal than X bp
            elif (monomer.length() >= minInternalMonomerSize):
                filteredMonomers.append(monomer)
            
            ## c) Discard internal monomers smaller than X bp

        ## 4. Convert monomer coordinates to inserted sequence space
        # NOT NEEDED as the unaligned sequence correspond to the leftmost end of the insert 
        # **unaligned**-----insert------ 

        ## 5. Add poly(T) info to the chain of alignments 
        nbMonomers = len(filteredMonomers)
        
        ## a) Single poly(T) 
        if nbMonomers == 1:

            monomer = filteredMonomers[0]

            ## Create PAF line containing poly(T) info
            firstAlignment = chain.alignments[0]
            fields = [firstAlignment.qName, firstAlignment.qLen, monomer.beg, monomer.end, None, 'poly(A/T)', 0, 0, 0, 0, 0, 0]
            alignment = formats.PAF_line(fields)
            chain.alignments.insert(0, alignment) ## Add to the chain
        
        ## b) Two poly(T) tracts. The most likely explanation of this structure is a short transduction
        # First poly(T) corresponds to the transcript polyA tail 
        # Second poly(T) corresponds to the source L1 element
        # Sequence in between corresponds to a short transduction (TD)
        # TTTTTTTTT-----TD-----TTTTTTTTT
        elif nbMonomers == 2:
            firstMonomer, secondMonomer = filteredMonomers
            
            ## Create PAF line containing TD info
            firstAlignment = chain.alignments[0]
            fields = [firstAlignment.qName, firstAlignment.qLen, firstMonomer.end, secondMonomer.end, None, 'Transduced3prime_UNK_UNK', 0, 0, 0, 0, 0, 0]
            alignment = formats.PAF_line(fields)
            chain.alignments.insert(0, alignment) ## Add to the chain

            ## Create PAF line containing poly(T) info
            fields = [firstAlignment.qName, firstAlignment.qLen, firstMonomer.beg, firstMonomer.end, None, 'poly(A/T)', 0, 0, 0, 0, 0, 0]
            alignment = formats.PAF_line(fields)
            chain.alignments.insert(0, alignment) ## Add to the chain


def infer_structure(chain):
    '''
    WORK ON PROGRESS
    
    Input:
        1. chain: Sequence chain of alignments over retrotranposon consensus sequences and/or transduced regions
        
    Output:
    ''' 
    #### 3. Collect transduction information  
    ## a) Not transduction
    if (insType == 'Solo') or (insType == 'Nested'):
        srcId = None 
        srcType = None
        srcCoord = None
        tdCoord = None
        tdLen = None
 
    ## b) 5' and 3' transductions
    #elif (insType == 'Partnered5,3prime'):

    ## c) Single transduction
    else:
        ## Select hits over the transduced region
        transductionHits = [alignment for alignment in chain.alignments if 'Transduced' in alignment.tName ]
        
        ## Collect src element id and coordinates from template name
        srcType = 'GERMLINE'   # Temporary
        srcId, srcCoord = transductionHits[0].tName.split(':')[1:]

        chrom, coords = srcCoord.split(':')
        beg, end = coords.split('-')

        ## 
        tdBeg = min([alignment.tBeg for alignment in transductionHits]) + int(beg)
        tdEnd = max([alignment.tEnd for alignment in transductionHits]) + int(beg)

        tdCoord = None
        tdLen = None
