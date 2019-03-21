'''
Module 'bamtools' - Contains functions for extracting data from bam files
'''

## DEPENDENCIES ##
# External
import pysam
import subprocess

# Internal
import log
import unix
import events
import gRanges
import formats
import sequences

## FUNCTIONS ##
def getREFS(bam):
    '''
    Get all references present in the bam file.

	Input:
		1. bam: indexed BAM file
	
	Output:
		1. refs: String containing all references from the bam file, separated by commas.
    '''
    bamFile = pysam.AlignmentFile(bam, 'rb')
    refs  = ','.join(bamFile.references)
    return refs

def SAM2BAM(SAM, outDir):
    '''
    Convert SAM file into sorted BAM and make BAM index

	Input:
		1. SAM: File containing alignments in SAM format

	Output:
		1. BAM_sorted: Sorted and indexed BAM file. BAM index located in the same directory with the extension '.bai'
    '''
    ## 0. Create logs directory 
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Convert SAM into BAM 
    BAM = outDir + '/alignments.bam'
    err = open(logDir + '/SAM2BAM.err', 'w') 
    command = 'samtools view -Sb ' + SAM + ' > ' + BAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'SAM2BAM'
        msg = 'SAM to BAM conversion failed' 
        log.step(step, msg)

    ## 2. Sort bam 
    BAM_sorted = outDir + '/alignments.sorted.bam'
    err = open(logDir + '/sort.err', 'w') 
    command = 'samtools sort ' + BAM + ' > ' + BAM_sorted
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'SORT'
        msg = 'BAM sorting failed' 
        log.step(step, msg)    

    ## 3. Index bam 
    BAM_index = outDir + '/alignments.sorted.bam.bai'
    err = open(logDir + '/index.err', 'w') 
    command = 'samtools index ' + BAM_sorted + ' > ' + BAM_index
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'INDEX'
        msg = 'BAM indexing failed' 
        log.step(step, msg)

    return BAM_sorted

def phred2ASCII(phred):
    '''
    Convert Phred quality scores into ASCII (Sanger format used in FASTQ)

	Input:
		1. phred: List containing per base standard Phred quality scores (from 0 to 93) as provided by pysam.AlignmentFile.query_qualities attribute.
	
	Output:
		1. ASCII: List containing per base phred quality scores encoded by ASCII characters 33 to 126 
    '''
    ASCII = [chr(x + 33) for x in phred]

    return ASCII

def BAM2FASTQ_entry(alignmentObj):
    '''
    Transform a BAM alignment into a FASTQ_entry object. 

	Input:
		1. alignmentObj: pysam.AlignedSegment object.
	
	Output:
		1. FASTQ_entry: formats.FASTQ_entry object
    '''
    ## 1. Convert Phred quality scores to ASCII 
    if not alignmentObj.query_qualities == None:
        ASCII = phred2ASCII(alignmentObj.query_qualities)
        ASCII = "".join(ASCII)
    else:
        qual = None

    ## 2. Obtain raw read and quality strings (Prior alignment)
    # a) Read mapped in reverse -> Make complementary reverse of the sequence and the reverse of the quality 
    if alignmentObj.is_reverse:
        seq = sequences.rev_complement(alignmentObj.query_sequence)

        if not alignmentObj.query_qualities == None:
            qual = ASCII[::-1]
                
    # b) Read mapped in forward
    else:
        seq = alignmentObj.query_sequence

        if not alignmentObj.query_qualities == None:
            qual = ASCII

    ## 3. Create FASTQ_entry object
    FASTQ_entry = formats.FASTQ_entry(alignmentObj.query_name, seq, '', qual)
    return FASTQ_entry

def binning(targetBins, bam, binSize, targetRefs):
    '''
    Split the genome into a set of genomic bins. Two possible binning approaches:
    1) Use predefined bins if bed file provided (targetBins)
    2) Non overlapping bins of a given size (binSize) for a set of target references (targetRefs). Reference length extracted from the input 'bam' file

    Input:
        1. targetBins: Bed file containing predefined bins OR None (in this case bins will be created de novo)
        2. bam: BAM file used to know the length of the target references 
        3. binSize: Binning size
        3. targetRefs: Comma separated list of target references

    Output:
        1. bins: List of bins. Each list item corresponds to a tuple (ref, beg, end)
    '''

    # A) Create bins de novo
    if targetBins == None:

        ## Split the reference genome into a set of genomic bins
        bins = makeGenomicBins(bam, binSize, targetRefs)

    # B) Read bins from bed file
    else:
        BED = formats.BED()
        BED.read(targetBins)
        bins = [ (line.ref, line.beg, line.end) for line in BED.lines]
    
    return bins


def makeGenomicBins(bam, binSize, targetRefs):
    '''
    Split the genome into a set of non overlapping bins of 'binSize' bp.

    Input:
        1. bam: BAM file
        2. binSize: size of the bins
        3. targetRefs: list of target references. None if all the references are considered

    Output:
        1. bins: List of non overlapping bins. Each list item corresponds to a tuple (ref, beg, end)
    '''

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, 'rb')

    ## Make dictionary with the length for each reference
    refLengths  = dict(list(zip(bamFile.references, bamFile.lengths)))

    ## Select target references
    if targetRefs != None:
        targetRefs = [str(i) for i in targetRefs]
        refLengths  = {ref: refLengths [ref] for ref in targetRefs}

    ## Split each reference into evenly sized bins
    bins = []

    # For each reference
    for ref, length in refLengths .items():

        ## Define bins boundaries
        boundaries = [boundary for boundary in range(0, length, binSize)]
        boundaries = boundaries + [length]

        ## Make bins
        for idx, beg in enumerate(boundaries):

            ## Skip last element from the list
            if beg < boundaries[-1]:
                end = boundaries[idx + 1]
                window = (ref, beg, end)
                bins.append(window)

    ## Close bam file
    bamFile.close()

    return bins


def collectSV(ref, binBeg, binEnd, bam, confDict, sample):
    '''
    Collect structural variant (SV) events from a genomic region.

    Input:
        1. ref: target referenge
        2. binBeg: target interval begin position
        3. binEnd: target interval end position
        4. bam: indexed BAM file
        5. confDict:
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * quality        -> True (sequence qualities available) or False (not available). 

        6. sample: type of sample (TUMOUR, NORMAL or None)

    Output:
        1. INS_events: list of INS objects
        2. DEL_events: list of DEL objects
        3. CLIPPING_left_events: list of CLIPPING objects on the left
        4. CLIPPING_left_events: list of CLIPPING objects on the right
        5. supportingReads: FASTQ (if qualities available) or FASTA (if qualities NOT available) object containing the reads supporting the SV events

        Note: * include secondary alignment filter???
    '''
    
    ## Initialize lists with SV
    INS_events = []
    DEL_events = []
    CLIPPING_left_events = []
    CLIPPING_right_events = []

    ## Initialize FASTA/FASTQ object containing supporting reads    
    # a) Quality available  
    if confDict['quality']:
        supportingReads = formats.FASTQ()

    # b) Quality not available
    else:
        supportingReads = formats.FASTA()

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
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
                                
            informative = False

            ## 1. Collect CLIPPINGS
            if 'CLIPPING' in confDict['targetSV']:

                CLIPPING_left, CLIPPING_right = collectCLIPPING(alignmentObj, confDict, sample)

                ## Select CLIPPING left breakpoints within the target genomic bin
                if CLIPPING_left != None:
                    overlapLen = gRanges.overlap(CLIPPING_left.beg, CLIPPING_left.end, binBeg, binEnd)

                    if overlapLen > 0:
                        informative = True
                        CLIPPING_left_events.append(CLIPPING_left)

                ## Select CLIPPING right breakpoints within the target genomic bin
                if CLIPPING_right != None:
                    overlapLen = gRanges.overlap(CLIPPING_right.beg, CLIPPING_right.end, binBeg, binEnd)

                    if overlapLen > 0:
                        informative = True
                        CLIPPING_right_events.append(CLIPPING_right)

            ## 2. Collect INDELS
            if ('INS' in confDict['targetSV']) or ('DEL' in confDict['targetSV']):
                INS_events_tmp, DEL_events_tmp = collectINDELS(alignmentObj, confDict, sample)
                
                ## Select INS events within the target genomic bin
                for INS in INS_events_tmp:
                    overlapLen = gRanges.overlap(INS.beg, INS.end, binBeg, binEnd)

                    if overlapLen > 0:
                        informative = True
                        INS_events.append(INS)

                ## Select DEL events within the target genomic bin
                for DEL in DEL_events_tmp:
                    overlapLen = gRanges.overlap(DEL.beg, DEL.end, binBeg, binEnd)
                    
                    if overlapLen > 0:
                        informative = True
                        DEL_events.append(DEL)

            ## 3. Add read to the FASTQ object if supports a SV event (DEL, INS or CLIPPING)
            if informative:
 
                # a) Quality available -> Add supporting read to FASTQ
                if confDict['quality']:

                    ## Transform the BAM alignment into a FASTQ_entry object.
                    entry = BAM2FASTQ_entry(alignmentObj)

                    ## Add entry to FASTQ
                    supportingReads.add(entry)

                # b) Quality not available -> Add supporting read to FASTA
                else:

                    ## Obtain raw read and quality strings (Prior alignment)
                    # a) Read mapped in reverse 
                    if alignmentObj.is_reverse:
                        seq = sequences.rev_complement(alignmentObj.query_sequence)

                    # b) Read mapped in forward
                    else:
                        seq = alignmentObj.query_sequence

                    ## Add sequence to FASTA
                    supportingReads.seqDict[alignmentObj.query_name] = seq

    ## Close 
    bamFile.close()

    # return sv candidates
    return INS_events, DEL_events, CLIPPING_left_events, CLIPPING_right_events, supportingReads


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
        clippingLeftObj = events.CLIPPING(alignmentObj, 'left', sample)

    ## Clipping > X bp at the right
    if ((lastOperation == 4) or (lastOperation == 5)) and (lastOperationLen >= confDict['minCLIPPINGlen']):
        clippingRightObj = events.CLIPPING(alignmentObj, 'right', sample)

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
    '''
    INS_events = []
    DEL_events = []

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
            insObj = events.INS(alignmentObj.reference_name, beg, end, insertLength, insertSeq, alignmentObj.query_name, sample)
            INS_events.append(insObj)

        ## b) DELETION to the reference >= Xbp
        if ('DEL' in confDict['targetSV']) and (operation == 2) and (length >= confDict['minINDELlen']):

            beg = posRef 
            end = posRef + length
            delObj = events.DEL(alignmentObj.reference_name, beg, end, length, alignmentObj.query_name, sample)
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


