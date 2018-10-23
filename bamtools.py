

'''
Module 'bamtools' - Contains Functions for extracting info froma a bam file
'''

## DEPENDENCIES ##
import pysam

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

    Output:
        1. INS_List: list of insertion objects    
        2. DEL_List: list of deletion objects    
        3. CLIP_List: list of clipped sequence objects
        4. readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read
    '''
    
    ## Initialize lists with SV
    INS_List = []    
    DEL_List = []    
    CLIPPED_left_List = []
    CLIPPED_right_List = []
    
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, beg, end)
            
    # Iterate over the read alignments
    for alignmentObj in iterator:

        ### Discard unmapped reads, PCR duplicates and reads with sequence not available
        if (alignmentObj.is_unmapped == False) and (alignmentObj.is_duplicate == False) and (alignmentObj.query_sequence != None):
            
            ## 1. Collect clipped read alignments
            if 'CLIPPING' in targetSV:

                clippingLeftObj, clippingRightObj = collectCLIPPING(alignmentObj, sample)

                if clippingLeftObj != None:

                    CLIPPED_left_List.append(clippingLeftObj)

                if clippingRightObj != None:
                    
                    CLIPPED_right_List.append(clippingRightObj)
    
            ## 2. Collect INDELS
            if ('INS' in targetSV) or ('DEL' in targetSV):
                INS_List_tmp, DEL_List_tmp = collectINDELS(alignmentObj, targetSV, sample)

                INS_List = INS_List + INS_List_tmp
                DEL_List = DEL_List + DEL_List_tmp                
                

def collectCLIPPING(alignmentObj, sample):
    '''

    Output:
        1. CLIP_List: list of clipped sequence objects
        2. readDict: dictionary containing the read alignments supporting the SV events. Format:
                key   -> read identifier
                value -> list of alignments for this read


    '''
    # Initialize as None
    clippingLeftObj, clippingRightObj = [None, None]
    
    # Select first and last operation from cigar to search for clipping
    firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
    lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]

    ## Clipping >= X bp at the left
    #  Note: soft (Operation=4) or hard clipped (Operation=5)     
    if ((firstOperation == 4) or (firstOperation == 5)) and (firstOperationLen >= 500):
        clippingLeftObj = CLIPPING(alignmentObj, 'left', sample)
        
    ## Clipping > X bp at the right
    if ((lastOperation == 4) or (lastOperation == 5)) and (lastOperationLen >= 500):
        clippingRightObj = CLIPPING(alignmentObj, 'right', sample) 
   
    return clippingLeftObj, clippingRightObj



def collectINDELS(alignmentObj, targetSV, sample):
    '''
    '''
    INS_List = []    
    DEL_List = []    

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

            insObj = INS(alignmentObj.reference_name, posRef, insertLength, insertSeq, readId, sample)
            INS_List.append(insObj)   

        ## b) DELETION to the reference >= Xbp
        if ('DEL' in targetSV) and (operation == 2) and (length >= 50):

            beg = posRef
            end = posRef + length

            delObj = DEL(alignmentObj.reference_name, beg, end, length, readId, sample)
            DEL_List.append(delObj)   

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

    return INS_List, DEL_List 


## CLASSES ##
class CLIPPING():
    '''
    '''
    def __init__(self, alignmentObj, clippedSide, sample):
        '''
        '''
        self.ref = alignmentObj.reference_name
        self.pos = alignmentObj.reference_start if clippedSide == 'left' else alignmentObj.reference_end
        self.clippedSide = clippedSide
        self.sample = sample

        ## Set supporting read id
        mate = '/1' if alignmentObj.is_read1 else '/2'
        self.readId = alignmentObj.query_name + mate  

        ## call functions
        self.type = ''
        self.clippingType(alignmentObj, clippedSide)

    def clippingType(self, alignmentObj, clippedSide):
        '''
        '''
        ### Extract operation
        ## a) Clipping at the begin 
        if clippedSide == 'left':

            operation =  alignmentObj.cigartuples[0][0]  

        ## b) Clipping at the end
        else:
            
            operation = alignmentObj.cigartuples[-1][0]

        ### Define if soft or had clipping
        self.type = 'soft' if (operation == 4) else 'hard'
       

class INS():
    """
    """
    def __init__(self, ref, pos, length, seq, readId, sample):

        """
        """
        self.ref = ref
        self.pos = pos
        self.length = length
        self.seq = seq
        self.readId = readId   
        self.sample = sample 

class DEL():
    """
    """
    def __init__(self, ref, beg, end, length, readId, sample):

        """
        """
        self.ref = ref
        self.beg = beg
        self.end = end
        self.length = length
        self.readId = readId   
        self.sample = sample 


 
