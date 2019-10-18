'''
Module 'events' - Contains classes for dealing with structural variation events at single read level
'''

## DEPENDENCIES ##
# External

# Internal
import annotation
import virus
import structures

###############
## FUNCTIONS ##
###############

def separate(events):
    '''
    Separate events according to their type into a dictionary containing multiple lists
    '''
    eventTypes = {}

    for event in events:

        ## Divide clippings into left and right 
        if event.type == 'CLIPPING':
            eventType = 'LEFT-CLIPPING' if event.clippedSide == 'left' else 'RIGHT-CLIPPING'

        elif event.type == 'DISCORDANT':
            eventType = 'MINUS-DISCORDANT' if event.side == 'MINUS' else 'PLUS-DISCORDANT'

        else:
            eventType = event.type

        # Initialize event type list    
        if eventType not in eventTypes:
            eventTypes[eventType] = []

        # Add event to list
        eventTypes[eventType].append(event)

    return eventTypes


def pick_flanking_seq_INS(readSeq, readPos, length, overhang):
    '''
    Pick inserted sequence + insertion flanking sequence from INS supporting read

    Input:
        1. readSeq: INS supporting read
        2. readPos: read position where INS start
        3. length: INS length
        4. overhang: number of flanking base pairs around the INS event to be collected from the supporting read sequence 

    Output:
        1. seq: piece of supporting read sequence spanning the INS event
        2. seqPos: INS beginning breakpoint position at the output sequence
    '''
    
    ## A) Enough sequence for extracting overhang at the begin
    #              readPos
    # ----*-----------|-------------------*---
    #      <--------->|<--MEI--><--------->
    #        overhang |seqPos    overhang
    if (readPos - overhang) > 0:
        begPos = readPos - overhang
        endPos = readPos + length + overhang
        seq = readSeq[begPos:endPos]
        seqPos = overhang

    ## B) Not enough sequence, so we would get into - values. Set lower bound
    #              readPos
    #          *------|--------------------*---
    #      <--------->|<--MEI--><--------->
    #       overhang  |seqPos     overhang
    else:
        endPos = readPos + length + overhang
        seq = readSeq[:endPos]
        seqPos = readPos
        
    return seq, seqPos


def pick_flanking_seq_DEL(readSeq, readPos, overhang):
    '''
    Pick sequence flanking the deletion breakpoints from DEL supporting read

    Input:
        1. readSeq: DEL supporting read
        2. readPos: read position where DEL start
        3. overhang: number of flanking base pairs around the DEL event to be collected from the supporting read sequence 

    Output:
        1. seq: piece of supporting read sequence spanning the DEL event
        2. seqPos: DEL breakpoint position at the output sequence
    '''
    ## Pick deletion flanking read sequence
    begPos = readPos - overhang
    begPos = begPos if begPos >= 0 else 0 # set lower bound to 0
    endPos = readPos + overhang
    seq = readSeq[begPos:endPos]
    seqPos = overhang

    return seq, seqPos


def pick_flanking_seq_CLIPPING(readSeq, readPos, clippedSide, overhang):
    '''
    Pick sequence flanking the clipping breakpoint from CLIPPING supporting read

    Input:
        1. readSeq: CLIPPING supporting read
        2. readPos: CLIPPING breakpoint position at the read
        3. clippedSide: clipping orientation (left or right)
        4. overhang: number of flanking base pairs around the CLIPPING breakpoint position to be collected from the supporting read sequence 

    Output:
        1. seq: piece of supporting read sequence spanning the CLIPPING breakpoint event
        2. seqPos: CLIPPING breakpoint position at the output sequence
    '''

    ## Take into account clipping orientation to pick read sequence (TO DO) 
    # a) Left clipping (.........*###################) (*, readPos) 
    #                   ----------<-overhang->* (*, endPos)
    if clippedSide == 'left':
        endPos = readPos + overhang
        seq = readSeq[:endPos]
        seqPos = readPos

    # b) Rigth clipping (###################*.........) (*, readPos) 
    #                         * <-overhang->------------------ (*, begPos)
    else:
        begPos = readPos - overhang
        begPos = begPos if begPos >= 0 else 0 # set lower bound to 0
        seq = readSeq[begPos:]
        seqPos = overhang

    return seq, seqPos


def determine_clippingType(alignmentObj, clippedSide):
    '''
    Determine if soft or hard clipping

    Input:
        1. alignmentObj: pysam read alignment object instance
        2. clippedSide: Clipped side relative to the read (left or right)

    Output:
        1. clippingType: soft or hard clipping
    '''
    ### Extract operation
    ## a) Clipping at the begin
    if clippedSide == 'left':
        operation =  alignmentObj.cigartuples[0][0]

    ## b) Clipping at the end
    else:
        operation = alignmentObj.cigartuples[-1][0]

    ### Determine if soft or hard clipping
    clippingType = 'soft' if (operation == 4) else 'hard'

    return clippingType

def determine_discordant_identity(discordants, repeatsBinDb, transducedBinDb):
    '''
    Determine discortant read pair identity based on the mapping position of anchor´s mate

    Input:
        1. discordants: list containing input discordant read pair events
        2. repeatsBinDb: dictionary containing annotated retrotransposons organized per chromosome (keys) into genomic bins (values)
        3. transducedBinDb: dictionary containing source element transduced regions (keys) into genomic bins (values)

    Output:
        1. discordantsIdentity: dictionary containing lists of discordant read pairs organized taking into account their orientation and if the mate aligns in an annotated retrotransposon 
                               This info is encoded in the dictionary keys as follows. Keys composed by 3 elements separated by '_':
                                
                                    - Orientation: read orientation (PLUS or MINUS)
                                    - Event type: DISCORDANT   
                                    - Type: identity type. It can be retrotransposon family (L1, Alu, ...), source element (22q, 5p, ...), viral strain (HPV, ...)
    '''
    
    ## 1. Assess if discordant read pairs support transduction insertion if transduction database provided
    if transducedBinDb is not None:
        discordantsTd = annotation.intersect_mate_annotation(discordants, transducedBinDb, 'cytobandId')

        ## Separate discordants matching from those not matching source elements
        discordants = []

        if 'PLUS_DISCORDANT_None' in discordantsTd:
            discordants = discordants + discordantsTd['PLUS_DISCORDANT_None']
            discordantsTd.pop('PLUS_DISCORDANT_None', None)

        if 'MINUS_DISCORDANT_None' in discordantsTd:
            discordants = discordants + discordantsTd['MINUS_DISCORDANT_None']
            discordantsTd.pop('MINUS_DISCORDANT_None', None)
    else:

        discordantsTd = {}

    ## 2. Assess if discordant read pairs support retrotransposons insertion if repeats database provided
    if repeatsBinDb is not None:
        discordantsRt = annotation.intersect_mate_annotation(discordants, repeatsBinDb, 'family')

        if 'PLUS_DISCORDANT_None' in discordantsRt:
            discordantsRt.pop('PLUS_DISCORDANT_None', None)

        if 'MINUS_DISCORDANT_None' in discordantsRt:
            discordantsRt.pop('MINUS_DISCORDANT_None', None)

    else:
        discordantsRt = {}

    ## 3. Merge discordant read pairs supporting RT and transduction insertions if transduction database provided    
    discordantsIdentity = structures.merge_dictionaries([discordantsTd, discordantsRt])

    '''
    ## 2. Assess if discordant read pairs support viral insertion

    # Create a list containing all discordant events:
    discordantEvents = []
    for eventType in discordantDict.keys():
    discordantEvents.extend(discordantDict[eventType])

    # a) Single sample mode
    if self.mode == "SINGLE":
        discordantEventsIdent = virus.is_virusSR(discordantEvents, self.bam, None, binDir, self.viralDbIndex)

    # b) Paired sample mode (tumour & matched normal)
    else:
        discordantEventsIdent = virus.is_virusSR(discordantEvents, self.bam, self.normalBam, binDir, self.viralDbIndex)
    '''

    return discordantsIdentity

def discordants2mates(discordants):
    '''
    Generate discordant objects for the mates of a input list of discordant events

    Input:
        1. discordants: list of discordant event objects

    Output:
        1. mates: list of discordant event objects corresponding to input discordant mates
    '''
    mates = []

    for discordant in discordants:
        pair = '2' if discordant.pair == '1' else '1'
        mate = DISCORDANT(discordant.mateRef, discordant.mateStart, discordant.mateStart, None, pair, discordant.readName, None, discordant.sample)

        mates.append(mate)

    return mates
    
def merge_INS(INS_list):
    '''
    Merge a set of adjacent INS events supported by the same read into a single one

    Input:
        1. INS_list: list of INS events to be merged 

    Output:
        1. merged: INS event resulting from merging input events
    '''

    ## 1. Sort INS by begin position
    INS_list.sort(key=lambda event: event.beg)

    ## 2. Define merged INS length
    first = INS_list[0]
    last = INS_list[-1]
    
    ## Temporary (remove if condition with normal bams. Needed here since testing bams have duplicated reads due to an error during the processing)
    if first.beg == last.beg:
        length = first.length 

    # Keep once error fixed
    else:

        ## 2.1 Compute fragmented alignment length on the reference genome 
        # INS_1----ref_1----INS_2--ref_2--INS_3
        # dist == ref_1 + ref_2 == INS_3.end - INS_1.beg
        dist = last.end - first.beg 
        
        ## 2.2 Compute total INS fragments length
        # INS_1----ref_1----INS_2--ref_2--INS_3
        fragmentsLen = sum([INS.length for INS in INS_list])

        ## 2.3 Compute total insertion length
        length = dist + fragmentsLen

    ## 3. Create merged INS
    merged = INS(first.ref, first.beg, first.end, length, first.readName, first.readSeq, first.readBkp, None, first.sample)

    return merged

#############
## CLASSES ##
#############

class INS():
    '''
    Short insertion class. Insertion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, readName, readSeq, readBkp, alignmentObj, sample):
        '''
        '''
        INS.number += 1 # Update instances counter
        self.id = 'INS_' + str(INS.number)
        self.type = 'INS'
        self.ref = str(ref)
        self.beg = int(beg) # beg==end. 0-based INS breakpoint
        self.end = int(end)
        self.length = int(length)
        self.readName = readName
        self.readSeq = readSeq
        self.readBkp = readBkp        
        self.sample = sample
        self.clusterId = None

        # Supporting read alignment properties:
        if alignmentObj is None:
            self.reverse = None
            self.secondary = None
            self.supplementary = None
            self.mapQual = None
            self.supplAlignment = None
        else:
            self.reverse = alignmentObj.is_reverse
            self.secondary = alignmentObj.is_secondary
            self.supplementary = alignmentObj.is_supplementary
            self.mapQual = alignmentObj.mapping_quality
            self.supplAlignment = alignmentObj.get_tag('SA') if alignmentObj.has_tag('SA') else None
    
    def pick_insert(self):
        '''
        Pick and return the inserted sequence 
        '''
        begPos = self.readBkp
        endPos = self.readBkp + self.length 
        insert = self.readSeq[begPos:endPos]
        
        return insert


class DEL():
    '''
    Short deletion class. Deletion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, readName, readSeq, readBkp, alignmentObj, sample):
        '''
        '''
        DEL.number += 1 # Update instances counter
        self.id = 'DEL_' + str(DEL.number)
        self.type = 'DEL'
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.length = int(length)
        self.readName = readName
        self.readSeq = readSeq
        self.readBkp = readBkp        
        self.sample = sample
        self.clusterId = None
    
        # Supporting read alignment properties:
        if alignmentObj is None:
            self.reverse = None
            self.secondary = None
            self.supplementary = None
            self.mapQual = None
            self.supplAlignment = None
        else:
            self.reverse = alignmentObj.is_reverse
            self.secondary = alignmentObj.is_secondary
            self.supplementary = alignmentObj.is_supplementary
            self.mapQual = alignmentObj.mapping_quality
            self.supplAlignment = alignmentObj.get_tag('SA') if alignmentObj.has_tag('SA') else None


class CLIPPING():
    '''
    Clipping class
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, clippedSide, readName, readSeq, readBkp, alignmentObj, sample):
        '''
        '''
        CLIPPING.number += 1 # Update instances counter
        self.id = 'CLIPPING_' + str(CLIPPING.number)
        self.type = 'CLIPPING'
        self.ref = str(ref)
        self.beg = int(beg) # beg==end. 0-based CLIPPING breakpoint
        self.end = int(end)
        self.length = length
        self.clippedSide = clippedSide
        self.clippingType = determine_clippingType(alignmentObj, self.clippedSide)
        self.readName = readName
        self.readSeq = readSeq
        self.readBkp = readBkp        
        self.sample = sample
        self.clusterId = None

        # Supporting read alignment properties:
        if alignmentObj is None:
            self.reverse = None
            self.secondary = None
            self.supplementary = None
            self.mapQual = None
            self.supplAlignment = None
        else:
            self.reverse = alignmentObj.is_reverse
            self.secondary = alignmentObj.is_secondary
            self.supplementary = alignmentObj.is_supplementary
            self.mapQual = alignmentObj.mapping_quality
            self.supplAlignment = alignmentObj.get_tag('SA') if alignmentObj.has_tag('SA') else None
            self.refLen = alignmentObj.reference_length


class DISCORDANT():
    '''
    Discordant class
    '''
    number = 0 # Number of instances
    
    def __init__(self, ref, beg, end, orientation, pair, readName, alignmentObj, sample):
        DISCORDANT.number += 1 # Update instances counter
        self.id = 'DISCORDANT_' + str(DISCORDANT.number)
        self.type = 'DISCORDANT'
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.orientation = orientation
        self.pair = pair
        self.readName = readName 
        self.sample = sample
        self.clusterId = None
        self.identity = None
        
        ## Mate info
        self.mateSeq = None

        if alignmentObj is None:
            self.mateRef = None
            self.mateStart = None           
        else:
            self.mateRef = alignmentObj.next_reference_name
            self.mateStart = alignmentObj.next_reference_start

    

