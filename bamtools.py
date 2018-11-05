'''
Module 'bamtools' - Contains functions for extracting data from bam files
'''

## DEPENDENCIES ##
# External
import pysam

# Internal
import log
import variants

## FUNCTIONS ##
def makeGenomicBins(bam, windowSize, targetRefs):
    '''
    Split the genome into a set of non overlapping windows of 'windowSize' bp. 

    Input:
        1. bam: BAM file
        2. windowSize: size of the windows
        3. targetRefs: list of target references. None if all the references are considered

    Output:
        1. windows: List of non overlapping windows. Each list item corresponds to a tuple (ref, beg, end)
    '''    

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, 'rb')

    ## Make dictionary with the length for each reference
    refLengths  = dict(list(zip(bamFile.references, bamFile.lengths)))

    ## Select target references
    if targetRefs != None:
        targetRefs = [str(i) for i in targetRefs]
        refLengths  = {ref: refLengths [ref] for ref in targetRefs} 

    ## Split each reference into evenly sized windows
    windows = []

    # For each reference
    for ref, length in refLengths .items():
    
        ## Define window boundaries
        boundaries = [boundary for boundary in range(0, length, windowSize)]    
        boundaries = boundaries + [length]
 
        ## Make windows
        for idx, beg in enumerate(boundaries):

            ## Skip last element from the list
            if beg < boundaries[-1]:
                end = boundaries[idx + 1]
                window = (ref, beg, end)
                windows.append(window)

    ## Close bam file
    bamFile.close()
    
    return windows


def collectSV(ref, beg, end, bam, confDict, sample):
    '''
    Collect structural variant (SV) candidates from a genomic region. 

    Input:
        1. ref: target referenge
        2. beg: target interval begin position
        3. end: target interval end position
        4. bam: indexed BAM file
        5. confDict: 
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght

        6. sample: type of sample (TUMOUR, NORMAL or None)

    Output:
        1. INS_events: list of INS objects
        2. DEL_events: list of DEL objects   
        3. CLIPPING_left_events: list of CLIPPING objects on the left
        4. CLIPPING_left_events: list of CLIPPING objects on the right
        5. [TO DO] readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read
    
    * include secondary alignment filter???
    '''
    
    ## Initialize lists with SV
    INS_events = []    
    DEL_events = []    
    CLIPPING_left_events = []
    CLIPPING_right_events = []
    
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, beg, end)
            
    # For each read alignment
    for alignmentObj in iterator:
        
        ### Select good quality alignments not having any of these properties:
        # - Unmapped 
        # - Secondary alignment (disabled, explore the possibility of including this filter)
        # - PCR or optical duplicate
        # - Read sequence not available
        # - Mapping quality < threshold
        MAPQ = int(alignmentObj.mapping_quality)   
 
        if (alignmentObj.is_unmapped == False) and (alignmentObj.is_duplicate == False) and (alignmentObj.query_sequence != None) and (MAPQ >= confDict['minMAPQ']):
            
            ## 1. Collect CLIPPINGS
            if 'CLIPPING' in confDict['targetSV']:

                clippingLeftObj, clippingRightObj = collectCLIPPING(alignmentObj, confDict, sample)

                if clippingLeftObj != None:
                    CLIPPING_left_events.append(clippingLeftObj)

                if clippingRightObj != None:
                    CLIPPING_right_events.append(clippingRightObj)
    
            ## 2. Collect INDELS
            if ('INS' in confDict['targetSV']) or ('DEL' in confDict['targetSV']):
                INS_events_tmp, DEL_events_tmp = collectINDELS(alignmentObj, confDict, sample)

                INS_events = INS_events + INS_events_tmp
                DEL_events = DEL_events + DEL_events_tmp            

    # return sv candidates
    return INS_events, DEL_events, CLIPPING_left_events, CLIPPING_right_events

def collectCLIPPING(alignmentObj, confDict, sample):
    '''
    For a read alignment check if the read is clipped on each side and return the corresponding clipping objects 

    Input: 
        1. alignmentObj: pysam read alignment object 
        2. confDict: configuration dictionary containing the following mandatory key value pairs:
            * minCLIPPINGlen -> minimum clipping lenght
        3. sample: type of sample (TUMOUR, NORMAL or None). Move to confDict
        
    Output:
        1. clippingLeftObj: CLIPPING object for left clipping (None if no clipping found) 
        2. clippingRightObj: CLIPPING object for right clipping (None if no clipping found)
        3. [TO DO] readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read

    Note: include filter to discard reads clipped at both their begin and end (useful with illumina data)
    '''
    # Initialize as None
    clippingLeftObj, clippingRightObj = [None, None]
    
    # Select first and last operation from cigar to search for clipping
    firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
    lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]

    ## Clipping >= X bp at the left
    #  Note: soft (Operation=4) or hard clipped (Operation=5)     
    if ((firstOperation == 4) or (firstOperation == 5)) and (firstOperationLen >= confDict['minCLIPPINGlen']):
        
        clippingLeftObj = variants.CLIPPING(alignmentObj, 'left', sample)
        
    ## Clipping > X bp at the right
    if ((lastOperation == 4) or (lastOperation == 5)) and (lastOperationLen >= confDict['minCLIPPINGlen']):
        clippingRightObj = variants.CLIPPING(alignmentObj, 'right', sample) 
   
    return clippingLeftObj, clippingRightObj


def collectINDELS(alignmentObj, confDict, sample):
    '''
    Collect insertions and deletions longer than a threshold that are completely spanned within an input read alignment

    Input: 
        1. alignmentObj: pysam read alignment object instance
        2. confDict: configuration dictionary containing the following mandatory key value pairs:
            * targetSV    -> list with target SV (INS: insertion; DEL: deletion)
            * minINDELlen -> minimum INS and DEL lenght
        3. sample: type of sample (TUMOUR, NORMAL or None). Move to confDict

    Output:
        1. INS_events: list of INS objects
        2. DEL_events: list of DEL objects
        3. [TO DO] readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read
    '''
    INS_events = []    
    DEL_events = []    

    ## Set read id
    mate = '/1' if alignmentObj.is_read1 else '/2'
    readId = alignmentObj.query_name + mate  

    ## Initialize positions at query and ref
    posQuery = 0        
    posRef = alignmentObj.reference_start

    # Iterate over the CIGAR
    for cigarTuple in alignmentObj.cigartuples:
                
        operation = int(cigarTuple[0])
        length = int(cigarTuple[1])

        ## a) INSERTION to the reference >= Xbp 
        if ('INS' in confDict['targetSV']) and (operation == 1) and (length >= confDict['minINDELlen']):
       
            beg = posRef 
            end = posRef 
            insertBeg = posQuery
            insertEnd = posQuery + length
            insertSeq = alignmentObj.query_sequence[insertBeg:insertEnd]
            insertLength = len(insertSeq)

            insObj = variants.INS(alignmentObj.reference_name, beg, end, insertLength, insertSeq, readId, sample)
            INS_events.append(insObj)   

        ## b) DELETION to the reference >= Xbp 
        if ('DEL' in confDict['targetSV']) and (operation == 2) and (length >= confDict['minINDELlen']):

            beg = posRef
            end = posRef + length

            delObj = variants.DEL(alignmentObj.reference_name, beg, end, length, readId, sample)
            DEL_events.append(delObj)   

        #### Update position over reference and read sequence
        ### a) Operations consuming query and reference
        # - Op M, tag 0, alignment match (can be a sequence match or mismatch)
        # - Op =, tag 7, sequence match
        # - Op X, tag 8, sequence mismatch
        if (operation == 0) or (operation == 7) or (operation == 8):
            posQuery += length
            posRef += length

        ### b) Operations only consuming query 
        # - Op I, tag 1, insertion to the reference
        # - Op S, tag 4, soft clipping (clipped sequences present in SEQ)
        elif (operation == 1) or (operation == 4):
            posQuery += length

        ### c) Operations only consuming reference 
        # - Op D, tag 2, deletion from the reference
        # - Op N, tag 3, skipped region from the reference
        elif (operation == 2) or (operation == 3):
            posRef += length

        ### d) Operations not consuming query nor reference 
        # - Op H, tag 5, hard clipping (clipped sequences NOT present in SEQ)
        # - Op P, tag 6, padding (silent deletion from padded reference)
        # Do not do anything

    return INS_events, DEL_events 


 
