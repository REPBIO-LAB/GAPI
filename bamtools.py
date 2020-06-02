'''
Module 'bamtools' - Contains functions for extracting data from bam files
'''

## DEPENDENCIES ##
# External
import pysam
import subprocess
import sys
from cigar import Cigar
import numpy as np
import time
import os
import Bio.SeqUtils
from Bio.SeqUtils import lcc

# Internal
import log
import unix
import events
import gRanges
import formats
import sequences

## FUNCTIONS ##
def get_refs(bam):
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


def get_ref_lengths(bam):
    '''
    Make dictionary containing the length for each reference

	Input:
		1. bam: indexed BAM file
	
	Output:
		1. lengths: Dictionary containing reference ids as keys and as values the length for each reference
    '''
    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, 'rb')

    ## Make dictionary with the length for each reference (move this code into a function)
    refLengths = dict(list(zip(bamFile.references, bamFile.lengths)))

    ## Close bam file
    bamFile.close()

    return refLengths


def alignment_length_cigar(CIGAR):
    '''
    Compute alignment on the reference length from CIGAR string

    Input:
        1. CIGAR: CIGAR string

    Output:
        1. alignmentLen: alignment on the reference length
    '''
    ## 1. Read CIGAR string using proper module
    cigarTuples = Cigar(CIGAR)

    ## 2. Iterate over the operations and compute the alignment length
    alignmentLen = 0

    for cigarTuple in list(cigarTuples.items()):

        length = int(cigarTuple[0])
        operation = cigarTuple[1]

        ### Update reference alignment length
        ## a) Operations consuming query and reference
        # - Op M, tag 0, alignment match (can be a sequence match or mismatch)
        # - Op =, tag 7, sequence match
        # - Op X, tag 8, sequence mismatch
        if (operation == 'M') or (operation == '=') or (operation == 'X'):
            alignmentLen += length

        ## b) Operations only consuming reference
        # - Op D, tag 2, deletion from the reference
        # - Op N, tag 3, skipped region from the reference
        elif (operation == 'D') or (operation == 'N'):
            alignmentLen += length
            
    return alignmentLen


def alignment_interval_query(CIGAR, orientation):
    '''
    Compute alignment on the reference length from CIGAR string

    Input:
        1. CIGAR: CIGAR string
        2. orientation: alignment orientation (+ or -) 

    Output:
        1. beg: begin position in query
        2. end: end position in query
    '''
    ## 1. Read CIGAR string using proper module
    cigar = Cigar(CIGAR)

    ## 2. Iterate over the operations and compute query alignment length and start position in query
    alignmentLen = 0
    counter = 0 # Count operations

    for cigarTuple in list(cigar.items()):

        length = int(cigarTuple[0])
        operation = cigarTuple[1]

        ## Set start position in query based on first operation 
        if counter == 0:

            # a) Soft or Hard clipping
            if (operation == 'S') or (operation == 'H'):
                startPos = length

            # b) No clipping
            else:
                startPos = 0
            
        #### Update query alignment length
        # - Op M, alignment match (can be a sequence match or mismatch)
        # - Op =, sequence match
        # - Op X, sequence mismatch
        # - Op I, insertion to the reference
        if (operation == 'M') or (operation == '=') or (operation == 'X') or (operation == 'I'):
            alignmentLen += length

        ## Update operations counter
        counter += 1

    ## 3. Compute alignment interval in raw query
    ## Compute read length
    readLen = len(cigar)

    # a) Query aligned in +
    if orientation == '+':
        beg = startPos
        end = startPos + alignmentLen

    # b) Query aligned in - (reversed complemented to align)
    else:
        beg = readLen - startPos - alignmentLen
        end = readLen - startPos
        
    return beg, end


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

def BAM2BED(BAM, outDir):
    '''
    Convert BAM file into BED using bedtools

	Input:
		1. BAM: Path to BAM file 
        2. outDir: Output directory

	Output:
		1. BED: Path to BED file
    '''
    ## 0. Create logs directory 
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Convert BAM into BED
    BED_path = outDir + '/alignments.bed'
    err = open(logDir + '/BAM2BED.err', 'w') 
    command = 'bedtools bamtobed -split -i ' + BAM + ' > ' + BED_path
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'BAM2BED'
        msg = 'BAM to BED conversion failed' 
        log.step(step, msg)
    
    ## 2. Add header to BED file
    header = "#ref \t beg \t end \t name \t score \t strand \n"
    with open(BED_path, 'r') as original: data = original.read()
    with open(BED_path, 'w') as modified: modified.write(header + data)

    return BED_path

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

def alignments2PAF(alignments):
    '''
    Convert as set of pysam aligned segments into a PAF object

	Input:
		1. alignments: list of aligned segments

	Output:
		1. PAF: PAF object containing alignments
    '''
    
    ## 1. Initialize PAF object
    PAF = formats.PAF()

    ## 2. Convert each aligned segment into a PAF_line object and add to PAF
    for alignment in alignments:

        # Discard unmapped sequences
        if not alignment.is_unmapped:
            
            strand = '-' if alignment.is_reverse else '+'
            fields = [alignment.query_name, alignment.infer_read_length(), alignment.query_alignment_start, alignment.query_alignment_end, strand, alignment.reference_name, alignment.reference_length, alignment.reference_start, alignment.reference_end, 0, 0, alignment.mapping_quality]
            line = formats.PAF_line(fields)
            PAF.alignments.append(line)

    return PAF


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
        1. bins: List of bins. Each list item corresponds to a list [ref, beg, end]
    '''

    # A) Create bins de novo
    if targetBins == None:

        ## Split the reference genome into a set of genomic bins
        bins = makeGenomicBins(bam, binSize, targetRefs)

    # B) Read bins from bed file
    else:
        BED = formats.BED()
        BED.read(targetBins, 'List', None)        
        bins = [ [line.ref, line.beg, line.end] for line in BED.lines]
    
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
    ## Obtain the length for each reference
    refLengths = get_ref_lengths(bam)

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

    return bins


def collectSV_paired(ref, binBeg, binEnd, tumourBam, normalBam, confDict):
    '''
    Collect structural variant (SV) events in a genomic bin from tumour and matched normal bam files

    Input:
        1. ref: target reference
        2. binBeg: bin begin
        3. binEnd: bin end
        4. tumourBam: indexed tumour BAM file
        5. normalBam: indexed normal BAM file
        6. confDict:
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght

    Output:
        1. eventsDict: dictionary containing list of SV events grouped according to the SV type (only those types included in confDict[targetSV]):
            * INS -> list of INS objects
            * DEL -> list of DEL objects
            * LEFT-CLIPPING -> list of left CLIPPING objects
            * RIGHT-CLIPPING -> list of right CLIPPING objects
            * DISCORDANT -> list of DISCORDANT objects  
    '''
    ## Search for SV events in the tumour
    eventsDict_T = collectSV(ref, binBeg, binEnd, tumourBam, confDict, 'TUMOUR', True)

    ## Search for SV events in the normal
    eventsDict_N = collectSV(ref, binBeg, binEnd, normalBam, confDict, 'NORMAL', True)

    ## Join tumour and normal lists
    eventsDict = {}

    for SV_type in eventsDict_T:        
        eventsDict[SV_type] = eventsDict_T[SV_type] + eventsDict_N[SV_type]

    return eventsDict


def collectSV(ref, binBeg, binEnd, bam, confDict, sample, supplementary):
    '''
    Collect structural variant (SV) events in a genomic bin from a bam file

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. confDict:
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings, DISCORDANT: discordant)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)
        6. sample: type of sample (TUMOUR, NORMAL or None)
        # TODO: explain supplementary AND DO THE SAME FOR PAIRED MODE

    Output:
        1. eventsDict: dictionary containing list of SV events grouped according to the SV type (only those types included in confDict[targetSV]):
            * INS -> list of INS objects
            * DEL -> list of DEL objects
            * LEFT-CLIPPING -> list of left CLIPPING objects
            * RIGHT-CLIPPING -> list of right CLIPPING objects
            * DISCORDANT -> list of DISCORDANT objects  
    
    NOTE: * include secondary alignment filter???
    '''
    # Define target interval
    targetInterval = (binBeg, binEnd)

    ## Initialize dictionary to store SV events
    eventsDict = {}

    if 'INS' in confDict['targetSV']:
        eventsDict['INS'] = []

    if 'DEL' in confDict['targetSV']:
        eventsDict['DEL'] = []

    if 'CLIPPING' in confDict['targetSV']:
        eventsDict['LEFT-CLIPPING'] = []
        eventsDict['RIGHT-CLIPPING'] = []
    '''
    if 'DISCORDANT' in confDict['targetSV']:
        eventsDict['DISCORDANT'] = []
    '''

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
    # For each read alignment
    for alignmentObj in iterator:

        ### 1. Filter out alignments based on different criteria:
        MAPQ = int(alignmentObj.mapping_quality) # Mapping quality

        ## Unmapped reads   
        if alignmentObj.is_unmapped == True:
            continue

        ## No query sequence available
        if alignmentObj.query_sequence == None:
            continue

        ## Aligments with MAPQ < threshold
        if (MAPQ < confDict['minMAPQ']):
            continue

        ## Duplicates filtering enabled and duplicate alignment
        if (confDict['filterDuplicates'] == True) and (alignmentObj.is_duplicate == True):
            continue

        # Filter supplementary alignments if FALSE. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)
        if supplementary == False and alignmentObj.is_supplementary == True:
            continue
        
        ## 2. Collect CLIPPINGS
        if 'CLIPPING' in confDict['targetSV']:

            left_CLIPPING, right_CLIPPING = collectCLIPPING(alignmentObj, confDict['minCLIPPINGlen'], targetInterval, sample)

            # Left CLIPPING found
            if left_CLIPPING != None:
                eventsDict['LEFT-CLIPPING'].append(left_CLIPPING)
        
            # Right CLIPPING found
            if right_CLIPPING != None:
                    
                eventsDict['RIGHT-CLIPPING'].append(right_CLIPPING)
                    
        ## 3. Collect INDELS
        if ('INS' in confDict['targetSV']) or ('DEL' in confDict['targetSV']):

            INDEL_events = collectINDELS(alignmentObj, confDict['targetSV'], confDict['minINDELlen'], targetInterval, confDict['overhang'], sample)

            # Add events to the pre-existing lists                
            for INDEL_type, events in INDEL_events.items():
                eventsDict[INDEL_type] = eventsDict[INDEL_type] + events

        '''
        ## 4. Collect DISCORDANT
        if 'DISCORDANT' in confDict['targetSV']:

            DISCORDANTS = collectDISCORDANT(alignmentObj, sample)

            # Add discordant events
            for discordant in DISCORDANTS:
                eventsDict['DISCORDANT'].append(discordant)
        '''
        
    ## Close 
    bamFile.close()

    # return sv candidates
    return eventsDict


def collectCLIPPING(alignmentObj, minCLIPPINGlen, targetInterval, sample):
    '''
    For a read alignment check if the read is clipped on each side and return the corresponding clipping objects

    Input:
        1. alignmentObj: pysam read alignment object
        2. minCLIPPINGlen: minimum clipping lenght
        3. targetInterval: tuple containing begin and end position of the target genomic interval to extract events from. If 'None' all clippings will be reported
        4. sample: type of sample (TUMOUR, NORMAL or None). 

    Output:
        1. left_CLIPPING: left CLIPPING object (None if no clipping found)
        2. right_CLIPPING: right CLIPPING object (None if no clipping found)
    '''    
    # Initialize as None
    left_CLIPPING, right_CLIPPING = (None, None)

    ## 2. Determine if discordant is mate 1 or 2
    if alignmentObj.is_read1:
        pair = '1'

    else:
        pair = '2'

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
            
            # Create CLIPPING object
            left_CLIPPING = events.CLIPPING(alignmentObj.reference_name, alignmentObj.reference_start, alignmentObj.reference_start, firstOperationLen, 'left', alignmentObj.query_name, alignmentObj.query_sequence, alignmentObj.query_alignment_start, pair, alignmentObj, sample)

    ## Clipping > X bp at the RIGHT
    if ((lastOperation == 4) or (lastOperation == 5)) and (lastOperationLen >= minCLIPPINGlen):
 
        ## Create CLIPPING object if:
        # a) No interval specified OR 
        # b) Clipping within target interval 
        if (targetInterval == None) or (gRanges.overlap(alignmentObj.reference_end, alignmentObj.reference_end, targetInterval[0], targetInterval[1])[0]):

            # Create CLIPPING object
            right_CLIPPING = events.CLIPPING(alignmentObj.reference_name, alignmentObj.reference_end, alignmentObj.reference_end, lastOperationLen, 'right', alignmentObj.query_name, alignmentObj.query_sequence, alignmentObj.query_alignment_end, pair, alignmentObj, sample)         

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
        1. INDEL_events: dictionary containing list of SV events grouped according to the type of INDEL (only those types included in targetSV):
            * INS -> list of INS objects
            * DEL -> list of DEL objects
    '''
    ## Initialize dict
    INDEL_events = {}

    if ('INS' in targetSV):
        INDEL_events['INS'] = []

    if ('DEL' in targetSV):   
        INDEL_events['DEL'] = []

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
                flankingSeq, bkpPos = (alignmentObj.query_sequence, posQuery) if overhang == None else events.pick_flanking_seq_INS(alignmentObj.query_sequence, posQuery, length, overhang)
                
                # Create INS object
                INS = events.INS(alignmentObj.reference_name, posRef, posRef, length, alignmentObj.query_name, flankingSeq, bkpPos, alignmentObj, sample)
                INDEL_events['INS'].append(INS)

        ## b) DELETION to the reference >= Xbp
        if ('DEL' in targetSV) and (operation == 2) and (length >= minINDELlen):
            
            ## Create DEL if:
            # a) No interval specified OR 
            # b) Deletion within target interval 
            if (targetInterval == None) or (gRanges.overlap(posRef, posRef + length, targetInterval[0], targetInterval[1])[0]):

                # Collect piece of sequence flanking the DEL event
                flankingSeq, bkpPos = (alignmentObj.query_sequence, posQuery) if overhang == None else events.pick_flanking_seq_DEL(alignmentObj.query_sequence, posQuery, overhang)

                # Create DEL object
                DEL = events.DEL(alignmentObj.reference_name, posRef, posRef + length, length, alignmentObj.query_name, flankingSeq, bkpPos, alignmentObj, sample)
                INDEL_events['DEL'].append(DEL)
                
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

    return INDEL_events

def collectDISCORDANT_paired(ref, binBeg, binEnd, tumourBam, normalBam, confDict, supplementary):
    '''
    For the two bam files given (test and normal), for each read alignment check if the read is discordant (not proper in pair) and return the corresponding discordant objects

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. tumourBam: indexed test BAM file
        5. normalBam: indexed normal BAM file
        6. confDict:
            * minMAPQ -> minimum mapping quality
            * filterDuplicates -> If True filter out reads labeled as duplicates in bam file.
        7. sample: type of sample (TUMOUR, NORMAL or None)
        8. supplementary: Filter out supplementary alignments if False. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)

    Output:
        1. DISCORDANTS: list containing input discordant read pair events
    '''
    ## Search for SV events in the tumour
    DISCORDANTS_SAMPLE = collectDISCORDANT(ref, binBeg, binEnd, tumourBam, confDict, 'TUMOUR', supplementary)

    ## Search for SV events in the normal
    DISCORDANTS_NORMAL = collectDISCORDANT(ref, binBeg, binEnd, normalBam, confDict, 'NORMAL', supplementary)

    ## Join tumour and normal lists
    DISCORDANTS = DISCORDANTS_SAMPLE + DISCORDANTS_NORMAL

    return DISCORDANTS

def collectDISCORDANT(ref, binBeg, binEnd, bam, confDict, sample, supplementary):
    '''
    In a given indexed bam file for each  a read alignment check if the read is discordant (not proper in pair) and return the corresponding discordant objects

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. confDict:
            * minMAPQ -> minimum mapping quality
            * filterDuplicates -> If True filter out reads labeled as duplicates in bam file.
        6. sample: type of sample (TUMOUR, NORMAL or None)
        7. supplementary: Filter out supplementary alignments if False. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)

    Output:
        1. DISCORDANTS: list containing input discordant read pair events
    '''

    # Define target interval
    targetInterval = (binBeg, binEnd)

    # Initialize discordant events list
    DISCORDANTS = []

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
    # For each read alignment
    for alignmentObj in iterator:

        ### 1. Filter out alignments based on different criteria:
        MAPQ = int(alignmentObj.mapping_quality) # Mapping quality

        ## Unmapped reads   
        if alignmentObj.is_unmapped == True:
            continue

        ## No query sequence available
        if alignmentObj.query_sequence == None:
            continue

        ## Aligments with MAPQ < threshold
        if (MAPQ < confDict['minMAPQ']):
            continue

        ## Duplicates filtering enabled and duplicate alignment
        if (confDict['filterDuplicates'] == True) and (alignmentObj.is_duplicate == True):
            continue

        # Filter out supplementary alignments if False. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)
        if supplementary == False and alignmentObj.is_supplementary == True:
            continue

        # Filter SMS reads
        firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
        lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]
        if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation == 4) or (lastOperation == 5)):
            continue

        # If not proper pair (== discordant)
        if not alignmentObj.is_proper_pair:

            ## 1. Determine discordant orientation
            # a) Minus
            if alignmentObj.is_reverse:
                orientation = 'MINUS'

            # b) Plus
            else:
                orientation = 'PLUS'

            ## 2. Determine if discordant is mate 1 or 2
            if alignmentObj.is_read1:
                pair = '1'

            else:
                pair = '2'

            ## 3. Determine number of alignment blocks
            operations = [t[0] for t in alignmentObj.cigartuples]
            nbBlocks = operations.count(3) + 1 

            ## 4. Create discordant event
            # A) Read aligning in a single block (WG or RNA-seq read no spanning a splice junction)
            # TODO SR: collectDISCORDANT: Think and check if is neccessary to take into account the number of blocks
            if nbBlocks == 1:
                DISCORDANT = events.DISCORDANT(alignmentObj.reference_name, alignmentObj.reference_start, alignmentObj.reference_end, orientation, pair, alignmentObj.query_name, alignmentObj, sample, alignmentObj.is_duplicate)
                DISCORDANTS.append(DISCORDANT)

            # B) Read alignning in multiple blocks (RNA-seq read spanning one or multiple splice junctions) -> Create one discordant event per block
            else:

                blockBeg = alignmentObj.reference_start
                blockEnd = blockBeg

                # For each operation
                for cigarTuple in alignmentObj.cigartuples:

                    operation = int(cigarTuple[0])
                    length = int(cigarTuple[1])

                    # a) End of the block -> End current block by creating a discordant event and Initiate a new block
                    if operation == 3:

                        # Create discordant event for the block
                        DISCORDANT = events.DISCORDANT(alignmentObj.reference_name, blockBeg, blockEnd, orientation, pair, alignmentObj.query_name, alignmentObj, sample, alignmentObj.is_duplicate)
                        DISCORDANTS.append(DISCORDANT)

                        # Initialize new block
                        blockBeg = blockEnd + length
                        blockEnd = blockEnd + length

                    # b) Extend current block
                    else:
                        blockEnd = blockEnd + length   

                ## End last block by creating a discordant
                DISCORDANT = events.DISCORDANT(alignmentObj.reference_name, blockBeg, blockEnd, orientation, pair, alignmentObj.query_name, alignmentObj, sample, alignmentObj.is_duplicate)
                DISCORDANTS.append(DISCORDANT)

    return DISCORDANTS


# def collectMatesSeq(events, tumourBam, normalBam, checkUnmapped, maxMAPQ):
#     '''
#     From a list of events, get the mate sequence for each event.
    
#     Input:
#         1. events: list of events
#         2. tumourBam
#         3. normalBam
#         4. checkUnmapped: boolean. True -> collect only the sequence of those mates that are unmapped, False -> dont take into account if they are unmapped or not.
#         5. maxMAPQ: collect only the sequence of those mates witn MAPQ < maxMAPQ
#     Output:
#         1. It doesnt return anything, just add the mate sequence to event.mateSeq attribute.
#     '''

#     counter = 0

#     ## 2. Open BAM file for reading
#     start = time.time()
#     bamFile = pysam.AlignmentFile(tumourBam, "rb")
#     end = time.time()
#     print("TIEMPO DE bamFile = pysam.AlignmentFile" + str(end - start))

#     '''
#     listMates = []
#     bamIndex = tumourBam + '.bai'
#     listMates.append(tumourBam)
#     listMates.append(bamIndex)
#     for event in events:
#         mateEnd = event.mateStart + 1
#         cadena = str(event.mateRef) + ";" + str(event.mateStart) + ";" + str(mateEnd) + ";" + event.readName.split('/')[0]
#         listMates.append(cadena)
#     #print (listMates)
    

    
#     #start = timer()
#     # TODO: MAPQ ESTA COMO FIJO, PONER PARAMETRO!! (pero funcionar funciona, comprobado)
#     matesSeqsString = libyay.massadd(listMates)
#     #end = timer()
#     #print("TIEMPO DE LIBYAY " + str(end - start))
#     print (type(matesSeqsString))
#     print(matesSeqsString)
#     matesSeqsList = matesSeqsString.split(";")

#     i = 0
#     for event in events:
#         if event.readName == matesSeqsList[i]:
#             pruebaMateSeq = matesSeqsList[i+1]
#             print('pruebaMateSeq')
#             print (pruebaMateSeq)
#             i = i+2
#     '''
    
#     # TODO: First and second conditions can be together, since they have same outcome.
#     #msg = 'LEN EVENTTTS: ' + str(len(events))
#     #log.subHeader(msg)
#     start = time.time()
#     for event in events:
#         counter += 1
#         if event.sample == None:
#             collectMateSeq(event, bamFile, checkUnmapped, maxMAPQ)
#         elif event.sample == 'TUMOUR':
#             collectMateSeq(event, tumourBam, checkUnmapped, maxMAPQ)
#         elif event.sample == 'NORMAL':
#             collectMateSeq(event, normalBam, checkUnmapped, maxMAPQ)

#     end = time.time()
#     print("TIEMPO DE event in events mates" + str(end - start))

#     msg = '[COUNTER OF collectMatesSeq LOOP] '+ str(counter)
#     log.subHeader(msg)


# def collectMateSeq(event, bamFile, checkUnmapped, maxMAPQ):
#     '''
#     Get mate sequence for an event.
    
#     Input:
#         1. event: DISCORDANT event
#         2. bam file
#         3. checkUnmapped: boolean. True -> Pick only those sequences that are unmmapped or with mapping quality < maxMAPQ. False -> Pick only those sequences with mapping quality < maxMAPQ
#         4. maxMAPQ: int. Pick only those reads with MAPQ < maxMAPQ.
#     Output:
#         1. It doesnt return anything, just add the mate sequence to event.mateSeq attribute.
#     '''
#     #msg = '[Start collectMateSeq]'
#     #log.subHeader(msg)
    
#     ## 1. Define bin coordinates based on mate position
#     mateRef = event.mateRef
#     mateStart = event.mateStart

#     binBegMate = mateStart
#     binEndMate = mateStart + 1

#     readName = event.readName.split('/')[0]
#     #readName = 'ST-E00181:606:HMWM2CCXY:3:2224:3366:11048'

#     ## 2. Open BAM file for reading
#     # Extract alignments   
#     iteratorMate = bamFile.fetch(mateRef, binBegMate, binEndMate)
#     #iteratorMate = bamFile.fetch('hs37d5', 12994840, 12994841)
    
#     # 3. For each read alignment
#     for alignmentObjMate in iteratorMate:
        
#         # Check if the aligment has same query_name but different orientation (to ensure that its the mate and not the read itself)
#         if readName == alignmentObjMate.query_name:
            
#             #msg = '[collectMateSeq: if readName == alignmentObjMate.query_name]' + str(binBegMate)
#             #log.subHeader(msg)
            
#             matePair = '1' if alignmentObjMate.is_read1 else '2'
            
#             if matePair != event.pair:
#                 #msg = 'if matePair != event.pair' + str(binBegMate)
#                 #log.subHeader(msg)

#                 MAPQ = int(alignmentObjMate.mapping_quality)
                
#                 # Pick only those sequences that are unmmapped or with mapping quality < maxMAPQ
#                 if checkUnmapped == True:
#                     if (alignmentObjMate.is_unmapped == True) or (MAPQ < maxMAPQ):
#                         #msg = 'if (alignmentObjMate.is_unmapped == True) or (MAPQ < maxMAPQ)' + str(binBegMate)
#                         #log.subHeader(msg)

#                         #if len(alignmentObjMate.query_sequence) > 100:
#                         event.mateSeq = alignmentObjMate.query_sequence
#                         print (event.mateSeq)
#                         break
                    
#                 # Pick only those sequences with mapping quality < maxMAPQ
#                 else:
#                     if MAPQ < maxMAPQ:
#                         #msg = 'if MAPQ < maxMAPQ' + str(binBegMate)
#                         #log.subHeader(msg)
#                         event.mateSeq = alignmentObjMate.query_sequence
#                         print (event.mateSeq)
#                         break

##############################
#TODO REPE EN FAST!!
def collectMatesSeq(events, tumourBam, normalBam, checkUnmapped, maxMAPQ):
    '''
    From a list of events, get the mate sequence for each event.
    
    Input:
        1. events: list of events
        2. tumourBam
        3. normalBam
        4. checkUnmapped: boolean. True -> collect only the sequence of those mates that are unmapped, False -> dont take into account if they are unmapped or not.
        5. maxMAPQ: collect only the sequence of those mates witn MAPQ < maxMAPQ
    Output:
        1. It doesnt return anything, just add the mate sequence to event.mateSeq attribute.
    '''

    listMates = []
 
    bamIndex = tumourBam + '.bai'
    listMates.append(tumourBam)
    listMates.append(bamIndex)
    counter = 0

    for event in events:
        mateEnd = event.mateStart + 1
        cadena = str(event.mateRef) + ";" + str(event.mateStart) + ";" + str(mateEnd) + ";" + str(event.beg)
        #print ('event.readName.split')
        #print (event.readName.split('/')[0])
        listMates.append(cadena)   
        counter += 1 
    print ('counter')
    print (counter)

    
    
    #listMates = ['/mnt/netapp2/mobilegenomes/0/1_projects/1020_CELL-LINES-PE/2_alignment/Ca-ski/Ca-ski.sorted.dedup.bam', '/mnt/netapp2/mobilegenomes/0/1_projects/1020_CELL-LINES-PE/2_alignment/Ca-ski/Ca-ski.sorted.dedup.bam.bai', '11;88756760;88756761;ST-E00181:606:HMWM2CCXY:1:1120:32329:31511', '11;88756760;88756761;ST-E00181:606:HMWM2CCXY:1:1120:31923:31582', '1;26783406;26783407;ST-E00181:606:HMWM2CCXY:2:1116:22983:20032', '8;37092223;37092224;ST-E00181:606:HMWM2CCXY:3:1110:8988:30457', '1;26783116;26783117;ST-E00181:606:HMWM2CCXY:2:1116:22983:20032', '1;26783116;26783117;ST-E00181:606:HMWM2CCXY:2:1116:22983:20032', '2;115861748;115861749;ST-E00181:606:HMWM2CCXY:2:2221:24809:8886', '1;26784068;26784069;ST-E00181:606:HMWM2CCXY:3:2110:3173:32971', '1;26784093;26784094;ST-E00181:606:HMWM2CCXY:3:1123:12459:64105', '1;26785508;26785509;ST-E00181:606:HMWM2CCXY:1:1118:22353:39686', '1;26785865;26785866;ST-E00181:606:HMWM2CCXY:2:2210:25175:54700']
    #listMates = ['/mnt/netapp2/mobilegenomes/0/1_projects/1020_CELL-LINES-PE/2_alignment/Ca-ski/Ca-ski.sorted.dedup.bam', '/mnt/netapp2/mobilegenomes/0/1_projects/1020_CELL-LINES-PE/2_alignment/Ca-ski/Ca-ski.sorted.dedup.bam.bai', 'hs37d5;12994840;12994840;1']
    print (listMates)

    # TODO: MAPQ ESTA COMO FIJO, PONER PARAMETRO!! (pero funcionar funciona, comprobado)
    start = time.time()
    #matesSeqsString = libyay.massadd(listMates)
    #print (os.getpid())
    libyay.massadd(listMates)
    matesSeqsString = ""
    end = time.time()
    print("TIEMPO DE massadd" + str(end - start))


    matesSeqsList = matesSeqsString.split(";")
    print ('LENGTHEVA')
    #print (len(matesSeqsList))
    print (matesSeqsList)

    #for i in len(events):
    #for i in range(0,counter):
        #events[i].mateSeq = matesSeqsList[i]
        #print ('NEW')
        #print (events[i].readName)
        #print (events[i].mateSeq)

    # i = 0
    # for event in events:
    #     if event.readName == matesSeqsList[i]:
    #         event.mateSeq = matesSeqsList[i+1]
    #         print (event.readName)
    #         print (event.mateSeq)
    #         i = i+2
            # BREAK Y QUITARLO DE LA LISTA
###############################


def average_MAPQ_reads_interval(ref, beg, end, readIds, bam):
    '''
    Retrieve a set of target reads in a genomic interval. Then compute their average MAPQ

    Input:
        1. ref: target reference
        2. beg: interval begin position
        3. end: interval end position
        4. readIds: list of target read ids
        5. bam: pysam filehandler for bam file 

    Output:
        1. avMAPQ: average mapping quality
    '''
    ## 1. Collect alignments in the input interval
    iterator = bam.fetch(ref, beg, end)

    ## 2. Select only those alignments corresponding to the target reads 
    alignments = []

    for alignment in iterator:
        
        if alignment.query_name in readIds:
            alignments.append(alignment)
        
    ## 3. Compute average MAPQ for mates
    qualities = [alignment.mapping_quality for alignment in alignments]
    avMAPQ = np.mean(qualities)
    
    return avMAPQ


def collectDiscodantsLowMAPQSeq(ref, binBeg, binEnd, bam, discordantMatesMaxMAPQ, discordantMatesCheckUnmapped, discordantMatesSupplementary, discordantMatesMaxBasePerc, discordantMatesMinLcc, outDir):
    '''
    Collecting read name and sequence of discordant low quality reads from all bam refs

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. discordantMatesMaxMAPQ: Maximum mapping quality used for collecting dicordant read mates.
        6. discordantMatesCheckUnmapped: Boolean. If True, when a dicordant read mate is unmapped, collect it no matter its MAPQ
        7. discordantMatesSupplementary: Boolean. When False, avoid collecting dicordant read mates that are supplementary alignments.
        8. discordantMatesMaxBasePerc: Maximum base percentage of discordant read mates sequences.
        9. discordantMatesMinLcc: Minimum local complexity of discordant read mates sequences.
        10. outDir: Output directory
    
    Output:
        1. Doesn't return anything. It creates a fasta file with discordant low quality reads from all bam refs.
    '''
    # TODO SR: Think if filterDuplicates step is neccesary in collectSeq method and implement it if so.
    #filterDuplicates = True

    ## Initialize dictionary to store SV events
    eventsSeqDict = {}

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
    # For each read alignment
    for alignmentObj in iterator:

        ### Filter out alignments based on different criteria:
        MAPQ = int(alignmentObj.mapping_quality) # Mapping quality

        ## No query sequence available
        if alignmentObj.query_sequence == None:
            continue

        ## Aligments with MAPQ < threshold
        if (MAPQ > discordantMatesMaxMAPQ):
            continue

        # TODO SR: Think if filterDuplicates step is neccesary in collectSeq method and implement it if so.
        ## Duplicates filtering enabled and duplicate alignment
        #if (confDict['filterDuplicates'] == True) and (alignmentObj.is_duplicate == True):
            #continue

        # Filter supplementary alignments if TRUE. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)
        if discordantMatesSupplementary == False and alignmentObj.is_supplementary == True:
            continue
        
        ## Collect DISCORDANT low quality reads in dictionary -> eventsSeqDict[readName] = readSequence

        if not alignmentObj.is_proper_pair:

            # Pick sequences that are unmmapped or with mapping quality < discordantMatesMaxMAPQ
            if discordantMatesCheckUnmapped == True:
                if (alignmentObj.is_unmapped == True) or (MAPQ < discordantMatesMaxMAPQ):
                    # Calculate base percentage
                    basePercs = sequences.baseComposition(alignmentObj.query_sequence)[1]
                    # Delete total value of base percentage result
                    del basePercs['total']
                    # Only those sequences with base percentage lower than 85 are collected:
                    if all(perc < discordantMatesMaxBasePerc for perc in basePercs.values()):
                        # Check local complexity of sequences:
                        complexity = Bio.SeqUtils.lcc.lcc_simp(alignmentObj.query_sequence)
                        # Only those sequences with local complexity lower than 1.49 are collected:
                        if complexity > discordantMatesMinLcc:
                            eventsSeqDict[alignmentObj.query_name]=alignmentObj.query_sequence

            # Pick sequences with mapping quality < discordantMatesMaxMAPQ
            else:
                if MAPQ < discordantMatesMaxMAPQ:
                    # Calculate base percentage
                    basePercs = sequences.baseComposition(alignmentObj.query_sequence)[1]
                    # Delete total value of base percentage result
                    del basePercs['total']
                    # Only those sequences with base percentage lower than 85 are collected:
                    if all(perc < discordantMatesMaxBasePerc for perc in basePercs.values()):
                        # Check local complexity of sequences:
                        complexity = Bio.SeqUtils.lcc.lcc_simp(alignmentObj.query_sequence)
                        # Only those sequences with local complexity lower than 1.49 are collected:
                        if complexity > discordantMatesMinLcc:
                            eventsSeqDict[alignmentObj.query_name]=alignmentObj.query_sequence
        
    ## Close 
    bamFile.close()
    collectVirusDir = outDir + '/COLLECT_VIRUS'
    # Set output FASTA file name
    allFastas_all = collectVirusDir + "/allFastas_all.fasta"
    # Create FASTA object
    seqsFastaObj= formats.FASTA()
    # Create FASTA dictionary
    seqsFastaObj.seqDict = eventsSeqDict
    # Write output FASTA
    seqsFastaObj.write(allFastas_all, 'append', True)

    # Delete useless variables
    del eventsSeqDict
    del seqsFastaObj

    return

def samtools_index_bam(BAM, outDir):
    '''
    Index bam file using samtools

    Input:
        1. BAM: Input bam file complete path.
    
    Output:
        1. Doesn't return anything. Creates bam index files.
    '''
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    command = 'samtools index ' + BAM
    err = open(logDir + '/samtools_index_bam.err', 'w') 
    status = subprocess.call(command, stderr=err, shell=True)

    return

def BAM2FastaDict(BAM):
    '''
    Pick reads names and sequences from bam file and write tham in FASTA format.

    Input:
        1. BAM: Input bam file complete path.

    Output:
        1. fastaDict: fastaDict[readName] = readSequence
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}
    # For each read alignment
    for alignmentObj in iterator:

        if alignmentObj.query_name in fastaDict.keys():
            fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
        else:
            fastaDict[alignmentObj.query_name] = []
            fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)



    return fastaDict

def filterBAM2FastaDict(BAM, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus):
    '''
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}

    # For each read alignment
    for alignmentObj in iterator:
        alignmentPass = False
        #numMatches = 0
        queryCoord = 0

        if not alignmentObj.is_unmapped:
            ctuples = alignmentObj.cigartuples
            allMatches = [t[1] for t in ctuples if t[0] == 0]
            totalMatch = sum (allMatches)
            if totalMatch >= minTotalMatchVirus:
                c = Cigar(alignmentObj.cigarstring)
                for citem  in list(c.items()):
                    # If cigar is query consuming, update query coordinates:
                    if citem[1] != 'M' and citem[1] != 'D' and citem[1] != '=':
                        queryCoord = queryCoord + int(citem[0])
                    elif citem[1] == 'M' or citem[1] == '=':
                        if citem[0] >= minParcialMatchVirus:
                            if (citem[0] <= maxMatchCheckMAPQVirus and alignmentObj.mapping_quality > minMAPQVirus) or citem[0] > maxMatchCheckMAPQVirus:
                                sequence = alignmentObj.query_sequence[queryCoord:(queryCoord + int(citem[0]))]
                                # Calculate base percentage
                                basePercs = sequences.baseComposition(sequence)[1]
                                # Delete total value of base percentage result
                                del basePercs['total']
                                # Only those sequences with base percentage lower than 85 are collected:
                                if all(perc < maxBasePercVirus for perc in basePercs.values()):
                                    # Calculate complexity
                                    complexity = Bio.SeqUtils.lcc.lcc_simp(sequence)
                                    if complexity > minLccVirus:
                                        alignmentPass = True
                                        break
                                    else:
                                        queryCoord = queryCoord + int(citem[0])
                                else:
                                    queryCoord = queryCoord + int(citem[0])
                        else:
                            queryCoord = queryCoord + int(citem[0])
        '''
        print (numMatches)
        if numMatches > 90: # TODO SR: put this as an option
            alignmentPass = True
        elif alignmentObj.mapping_quality >= viralBamMAPQ and numMatches >= int(selectedPartialMatches):
            alignmentPass = True
        # New condition:
        # TODO SR: Put as options!!!
        elif alignmentObj.mapping_quality >= 40 and numMatches >= 60:
            alignmentPass = True
        '''
        if alignmentPass == True:
            # Add to fasta dict
            if alignmentObj.query_name in fastaDict.keys():
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
            else:
                fastaDict[alignmentObj.query_name] = []
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)

    return fastaDict
