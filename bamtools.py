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
def makeGenomicBins(bam, windowSize, targetRefList):
    '''
    Split the genome into a set of non overlapping windows of 'windowSize' bp. 

    Input:
        1. bam: BAM file
        2. windowSize: size of the windows
        3. targetRefList: list of target references. None if all the references are considered

    Output:
        1. windowsList: List of non overlapping windows. Each list item corresponds to a tuple (ref, beg, end)
    '''    

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, 'rb')

    ## Make dictionary with the length for each reference
    refLengthsDict = dict(list(zip(bamFile.references, bamFile.lengths)))

    ## Select target references
    if targetRefList != None:
        targetRefList = [str(i) for i in targetRefList]
        refLengthsDict = {ref: refLengthsDict[ref] for ref in targetRefList} 

    ## Split each reference into evenly sized windows
    windowsList = []

    # For each reference
    for ref, length in refLengthsDict.items():
    
        ## Define window boundaries
        boundariesList = [boundary for boundary in range(0, length, windowSize)]    
        boundariesList = boundariesList + [length]
 
        ## Make windows
        for idx, beg in enumerate(boundariesList):

            ## Skip last element from the list
            if beg < boundariesList[-1]:
                end = boundariesList[idx + 1]
                windowTuple = (ref, beg, end)
                windowsList.append(windowTuple)

    ## Close bam file
    bamFile.close()
    
    return windowsList


def collectSV(ref, beg, end, bam, targetSV, sample):
    '''
    Collect structural variant (SV) candidates from a genomic region. 

    Input:
        1. ref: target referenge
        2. beg: target interval begin position
        3. end: target interval end position
        4. bam: indexed BAM file
        5. targetSV: list of target SV to search for (3 possible SV types: insertions (INS), deletions (DEL) and clippings (CLIPPING))
        6. sample: type of sample (TUMOUR, NORMAL or None)

    Output:
        1. INS_list: list of INS objects
        2. DEL_list: list of DEL objects   
        3. CLIPPING_left_list: list of CLIPPING objects on the left
        4. CLIPPING_left_list: list of CLIPPING objects on the right
        5. [TO DO] readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read


    Note: some filters are hardcoded (i.e. minimum lengths). At one point improve to take a dictionary as input specifying the filters configuration.
    '''
    
    ## Initialize lists with SV
    INS_list = []    
    DEL_list = []    
    CLIPPING_left_list = []
    CLIPPING_right_list = []
    
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, beg, end)
            
    # For each read alignment
    for alignmentObj in iterator:
        
        ### Discard unmapped reads, PCR duplicates and reads with sequence not available
        if (alignmentObj.is_unmapped == False) and (alignmentObj.is_duplicate == False) and (alignmentObj.query_sequence != None):
            
            ## 1. Collect CLIPPINGS
            if 'CLIPPING' in targetSV:

                clippingLeftObj, clippingRightObj = collectCLIPPING(alignmentObj, sample)

                if clippingLeftObj != None:
                    CLIPPING_left_list.append(clippingLeftObj)

                if clippingRightObj != None:
                    CLIPPING_right_list.append(clippingRightObj)
    
            ## 2. Collect INDELS
            if ('INS' in targetSV) or ('DEL' in targetSV):
                INS_list_tmp, DEL_list_tmp = collectINDELS(alignmentObj, targetSV, sample)

                INS_list = INS_list + INS_list_tmp
                DEL_list = DEL_list + DEL_list_tmp            

    # return sv candidates
    return INS_list, DEL_list, CLIPPING_left_list, CLIPPING_right_list

def collectCLIPPING(alignmentObj, sample):
    '''
    Check for a read alignment if the read is clipped on both sides and return the corresponding clipping objects 

    Input: 
        1. alignmentObj: pysam read alignment object 
        2. sample: type of sample (TUMOUR, NORMAL or None)
        
    Output:
        1. clippingLeftObj: CLIPPING object for left clipping (None if no clipping found) 
        2. clippingRightObj: CLIPPING object for right clipping (None if no clipping found)
        3. [TO DO] readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read

    Note: some filters are hardcoded (i.e. minimum lengths). At one point improve to take a dictionary as input specifying the filters configuration.
    Note: consider to include a minimum mapping quality filter
    Note: include filter to discard reads clipped at both their begin and end (useful with illumina data)
    '''
    # Initialize as None
    clippingLeftObj, clippingRightObj = [None, None]
    
    # Select first and last operation from cigar to search for clipping
    firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
    lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]

    ## Clipping >= X bp at the left
    #  Note: soft (Operation=4) or hard clipped (Operation=5)     
    if ((firstOperation == 4) or (firstOperation == 5)) and (firstOperationLen >= 500):
        clippingLeftObj = variants.CLIPPING(alignmentObj, 'left', sample)
        
    ## Clipping > X bp at the right
    if ((lastOperation == 4) or (lastOperation == 5)) and (lastOperationLen >= 500):
        clippingRightObj = variants.CLIPPING(alignmentObj, 'right', sample) 
   
    return clippingLeftObj, clippingRightObj


def collectINDELS(alignmentObj, targetSV, sample):
    '''
    Collect insertions and deletions longer than a threshold that are completely spanned within an input read alignment

    Input: 
        1. alignmentObj: pysam read alignment object instance
        2. sample: type of sample (TUMOUR, NORMAL or None)

    Output:
        1. INS_list: list of INS objects
        2. DEL_list: list of DEL objects
        3. [TO DO] readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read

    Note: some filters are hardcoded (i.e. minimum lengths). At one point improve to take a dictionary as input specifying the filters configuration.
    Note: consider to include a minimum mapping quality filter
    '''
    INS_list = []    
    DEL_list = []    

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
        if ('INS' in targetSV) and (operation == 1) and (length >= 50):
       
            insertBeg = posQuery
            insertEnd = posQuery + length
            insertSeq = alignmentObj.query_sequence[insertBeg:insertEnd]
            insertLength = len(insertSeq)

            insObj = variants.INS(alignmentObj.reference_name, posRef, insertLength, insertSeq, readId, sample)
            INS_list.append(insObj)   

        ## b) DELETION to the reference >= Xbp
        if ('DEL' in targetSV) and (operation == 2) and (length >= 50):

            beg = posRef
            end = posRef + length

            delObj = variants.DEL(alignmentObj.reference_name, beg, end, length, readId, sample)
            DEL_list.append(delObj)   

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

    return INS_list, DEL_list 


 
