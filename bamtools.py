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


def collectSV_paired(ref, binBeg, binEnd, tumourBam, normalBam, confDict):
    '''
    Collect structural variant (SV) events in a genomic bin from tumour and matched normal bam files

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. tumourBam: indexed tumour BAM file
        5. normalBam: indexed normal BAM file
        6. confDict:
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * quality        -> True (sequence qualities available) or False (not available). 

    Output:
        1. eventsDict: dictionary containing list of SV events grouped according to the SV type:
            * INS -> list of INS objects
            * DEL -> list of DEL objects
            * LEFT-CLIPPING -> list of left CLIPPING objects
            * RIGHT-CLIPPING -> list of right CLIPPING objects
    '''
    ## Search for SV events in the tumour
    eventsDict_T = collectSV(ref, binBeg, binEnd, tumourBam, confDict, 'TUMOUR')

    ## Search for SV events in the normal
    eventsDict_N = collectSV(ref, binBeg, binEnd, normalBam, confDict, 'NORMAL')

    ## Join tumour and normal lists
    eventsDict = {}

    for SV_type in eventsDict_T:        
        eventsDict[SV_type] = eventsDict_T[SV_type] + eventsDict_N[SV_type]

    return eventsDict


def collectSV(ref, binBeg, binEnd, bam, confDict, sample):
    '''
    Collect structural variant (SV) events in a genomic bin from a bam file

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. confDict:
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * overhang       -> Number of flanking base pairs around the SV event to be collected from the supporting read. If 'None' the complete read sequence will be collected)
        6. sample: type of sample (TUMOUR, NORMAL or None)

    Output:
        1. eventsDict: dictionary containing list of SV events grouped according to the SV type:
            * INS -> list of INS objects
            * DEL -> list of DEL objects
            * LEFT-CLIPPING -> list of left CLIPPING objects
            * RIGHT-CLIPPING -> list of right CLIPPING objects
        Note: * include secondary alignment filter???
    '''
    # Define target interval
    targetInterval = (binBeg, binEnd)

    ## Initialize dictionary to store SV events
    eventsDict = {}

    for SV_type in confDict['targetSV']:

        # a) Dividide clipping events into left and right clippings
        if (SV_type == 'CLIPPING'):
            eventsDict['LEFT-CLIPPING'] = []
            eventsDict['RIGHT-CLIPPING'] = []

        # b) Other types of events
        else:
            eventsDict[SV_type] = []

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
                                
            ## 1. Collect CLIPPINGS
            if 'CLIPPING' in confDict['targetSV']:

                left_CLIPPING, right_CLIPPING = collectCLIPPING(alignmentObj, confDict['minCLIPPINGlen'], targetInterval, confDict['overhang'], sample)

                # Left CLIPPING found
                if left_CLIPPING != None:
                    eventsDict['LEFT-CLIPPING'].append(left_CLIPPING)
        
                # Right CLIPPING found
                if right_CLIPPING != None:
                    eventsDict['RIGHT-CLIPPING'].append(right_CLIPPING)

            ## 2. Collect INDELS
            if ('INS' in confDict['targetSV']) or ('DEL' in confDict['targetSV']):

                INS_events, DEL_events = collectINDELS(alignmentObj, confDict['targetSV'], confDict['minINDELlen'], targetInterval, confDict['overhang'], sample)
                
                # Add events to the pre-existing lists
                eventsDict['INS'] = eventsDict['INS'] + INS_events
                eventsDict['DEL'] = eventsDict['DEL'] + DEL_events

    ## Close 
    bamFile.close()

    # return sv candidates
    return eventsDict


def collectCLIPPING(alignmentObj, minCLIPPINGlen, targetInterval, overhang, sample):
    '''
    For a read alignment check if the read is clipped on each side and return the corresponding clipping objects

    Input:
        1. alignmentObj: pysam read alignment object
        2. minCLIPPINGlen: minimum clipping lenght
        3. targetInterval: tuple containing begin and end position of the target genomic interval to extract events from. If 'None' all clippings will be reported
        4. overhang: Number of flanking base pairs around the SV event to be collected from the supporting read. If 'None' the complete read sequence will be collected)
        5. sample: type of sample (TUMOUR, NORMAL or None). 

    Output:
        1. left_CLIPPING: left CLIPPING object (None if no clipping found)
        2. right_CLIPPING: right CLIPPING object (None if no clipping found)
    '''    
    # Initialize as None
    left_CLIPPING, right_CLIPPING = (None, None)

    # Select first and last operation from cigar to search for clipping
    firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
    lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]

    ## Clipping >= X bp at the LEFT
    #  Note: soft (Operation=4) or hard clipped (Operation=5)     
    if ((firstOperation == 4) or (firstOperation == 5)) and (firstOperationLen >= minCLIPPINGlen):

        ## Create CLIPPING object if:
        # a) No interval specified OR 
        # b) Clipping within target interval 
        if (targetInterval == None) or (gRanges.overlap(alignmentObj.reference_start, alignmentObj.reference_start, targetInterval[0], targetInterval[1])[0]):
            
            # Collect piece of sequence flanking the clipping event
            flankingSeq = alignmentObj.query_sequence if overhang == None else events.pick_flanking_seq_CLIPPING(alignmentObj.query_sequence, alignmentObj.query_alignment_start, 'left', overhang)
            
            # Create CLIPPING object
            left_CLIPPING = events.CLIPPING(alignmentObj.reference_name, alignmentObj.reference_start, alignmentObj.reference_start, 'left', flankingSeq, alignmentObj, sample)

    ## Clipping > X bp at the RIGHT
    if ((lastOperation == 4) or (lastOperation == 5)) and (lastOperationLen >= minCLIPPINGlen):

        ## Create CLIPPING object if:
        # a) No interval specified OR 
        # b) Clipping within target interval 
        if (targetInterval == None) or (gRanges.overlap(alignmentObj.reference_end, alignmentObj.reference_end, targetInterval[0], targetInterval[1])[0]):

            # Collect piece of sequence flanking the clipping event
            flankingSeq = alignmentObj.query_sequence if overhang == None else events.pick_flanking_seq_CLIPPING(alignmentObj.query_sequence, alignmentObj.query_alignment_end, 'right', overhang)

            # Create CLIPPING object
            right_CLIPPING = events.CLIPPING(alignmentObj.reference_name, alignmentObj.reference_end, alignmentObj.reference_end, 'right', flankingSeq, alignmentObj, sample)

    return left_CLIPPING, right_CLIPPING


def collectINDELS(alignmentObj, targetSV, minINDELlen, targetInterval, overhang, sample):
    '''
    Collect insertions and deletions longer than a threshold that are completely spanned within an input read alignment

    Input:
        1. alignmentObj: pysam read alignment object instance
        2. targetSV: list with target SV (INS: insertion; DEL: deletion)
        3. minINDELlen: minimum INS and DEL lenght
        4. targetInterval: tuple containing begin and end position of the target genomic interval to extract events from. If 'None' all the events spanned by the read alignment will be reported
        5. overhang: Number of flanking base pairs around the SV event to be collected from the supporting read. If 'None' the complete read sequence will be collected)        
        6. sample: type of sample (TUMOUR, NORMAL or None). 

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
        if ('INS' in targetSV) and (operation == 1) and (length >= minINDELlen):

            ## Create INS if:
            # a) No interval specified OR 
            # b) Insertion within target interval 
            if (targetInterval == None) or (gRanges.overlap(posRef, posRef, targetInterval[0], targetInterval[1])[0]):

                # Collect piece of sequence flanking the INS event
                flankingSeq = alignmentObj.query_sequence if overhang == None else events.pick_flanking_seq_INS(alignmentObj.query_sequence, posQuery, length, overhang)

                # Create INS object
                INS = events.INS(alignmentObj.reference_name, posRef, posRef, length, flankingSeq, alignmentObj, sample)
                INS_events.append(INS)

        ## b) DELETION to the reference >= Xbp
        if ('DEL' in targetSV) and (operation == 2) and (length >= minINDELlen):
            
            ## Create DEL if:
            # a) No interval specified OR 
            # b) Deletion within target interval 
            if (targetInterval == None) or (gRanges.overlap(posRef, posRef + length, targetInterval[0], targetInterval[1])[0]):

                # Collect piece of sequence flanking the DEL event
                flankingSeq = alignmentObj.query_sequence if overhang == None else events.pick_flanking_seq_DEL(alignmentObj.query_sequence, posQuery, overhang)

                # Create DEL object
                DEL = events.DEL(alignmentObj.reference_name, posRef, posRef + length, length, flankingSeq, alignmentObj, sample)
                DEL_events.append(DEL)
                
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


