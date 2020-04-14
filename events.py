'''
Module 'events' - Contains classes for dealing with structural variation events at single read level
'''

## DEPENDENCIES ##
# External

# Internal
import annotation
import virus
import structures
import bamtools
import formats
from cigar import Cigar
import gRanges


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
            eventType = 'MINUS-DISCORDANT' if event.orientation == 'MINUS' else 'PLUS-DISCORDANT'

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


def determine_discordant_identity(discordants, repeatsBinDb, transducedBinDb, bam, normalBam, binDir, targetINT2Search):
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
    # TODO SR: If both, RT and virus, are analysed check if it is a virus only if it is not a RT.
    discordantsIdentity1 = {}
    if 'ME' in targetINT2Search:
        ## 1. Assess if discordant read pairs support transduction insertion if transduction database provided
        if transducedBinDb is not None:
            discordantsTd = annotation.intersect_mate_annotation(discordants, transducedBinDb, 'cytobandId')

            ## Separate discordants matching from those not matching source elements
            discordants = []

            if 'PLUS-DISCORDANT-None' in discordantsTd:
                discordants = discordants + discordantsTd['PLUS-DISCORDANT-None']
                discordantsTd.pop('PLUS-DISCORDANT-None', None)

            if 'MINUS-DISCORDANT-None' in discordantsTd:
                discordants = discordants + discordantsTd['MINUS-DISCORDANT-None']
                discordantsTd.pop('MINUS-DISCORDANT-None', None)
        else:

            discordantsTd = {}

        ## 2. Assess if discordant read pairs support retrotransposons insertion if repeats database provided
        discordantsNone = {}
        discordantsNone['PLUS-DISCORDANT-None'] = []
        discordantsNone['MINUS-DISCORDANT-None'] = []
        if repeatsBinDb is not None:
            discordantsRt = annotation.intersect_mate_annotation(discordants, repeatsBinDb, 'family')

            if 'PLUS-DISCORDANT-None' in discordantsRt:
                discordantsNone['PLUS-DISCORDANT-None'].extend(discordantsRt['PLUS-DISCORDANT-None'])
                discordantsRt.pop('PLUS-DISCORDANT-None', None)

            if 'MINUS-DISCORDANT-None' in discordantsRt:
                discordantsNone['MINUS-DISCORDANT-None'].extend(discordantsRt['MINUS-DISCORDANT-None'])
                discordantsRt.pop('MINUS-DISCORDANT-None', None)

        else:
            discordantsRt = {}
        ## 3. Merge discordant read pairs supporting RT and transduction insertions if transduction database provided    
        discordantsIdentity1 = structures.merge_dictionaries([discordantsTd, discordantsRt, discordantsNone])

        # Remaining discordants for viruses:
        # NOTE SR: If viruses and MEs are searched at once, those discordants identified as MEs are NOT checked for viral identities.
        discordants = discordantsNone['PLUS-DISCORDANT-None'] + discordantsNone['MINUS-DISCORDANT-None']
    


    ## 2. Assess if discordant read pairs support viral insertion

    # Create a list containing all discordant events:
    #discordantEvents = []
    #for eventType in discordantDict.keys():
    #discordantEvents.extend(discordantDict[eventType])

    # TODO SR: Make viruses optional in determine_discordant_identity
    # TODO SR: Add RT analysis in determine_discordant_identity
    discordantEventsIdent = {}

    # Remove discordants events that were clasified as MEs (MAYBE CHANGE THIS WHEN ADDED NEW PARAMETER)

        
    if 'VIRUS' in targetINT2Search:
        discordantEventsIdent = virus.is_virusSR(discordants)
    #discordantsIdentity = structures.merge_dictionaries([discordantEventsIdent, discordantsIdentity1])
    discordantEventsIdent.update(discordantsIdentity1)

    #return discordantsIdentity
    return discordantEventsIdent

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

def events2Dict(events):
    '''
    Organize a set of input events into a dictionary based on their reference

    Input:
        1. events: list of input events. Events should contain the ref attribute

    Output:
        1. eventsDict: dictionary with the references as keys and the list of events from each reference as values
    '''    
    eventsDict = {}

    ## For each event
    for event in events:

        # Initialize list if needed
        if event.ref not in eventsDict:
            eventsDict[event.ref] = []
        
        # Add event to the list
        eventsDict[event.ref].append(event)

    return eventsDict

def events2nestedDict(events, eventType):
    '''
    Organize a set of input events into a nested dictionary

    Input:
        1. events: list of input events. Events should contain ref attribute
        2. eventType: type of events. Used a key in inner dictionary

    Output:
        1. eventsDict: nested dictionary with first level keys corresponding to references, second level keys to
                       eventType and the corresponding list of events as values
    '''    
    eventsDict = {}

    ## For each event
    for event in events:

        # Initialize list if needed
        if event.ref not in eventsDict:
            eventsDict[event.ref] = {}
            eventsDict[event.ref][eventType] = []
        
        # Add event to the list
        eventsDict[event.ref][eventType].append(event)

    return eventsDict

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
            self.CIGAR = None
            self.reverse = None
            self.secondary = None
            self.supplementary = None
            self.mapQual = None
            self.supplAlignment = None

        else:
            self.CIGAR = alignmentObj.cigarstring
            self.reverse = alignmentObj.is_reverse
            self.secondary = alignmentObj.is_secondary
            self.supplementary = alignmentObj.is_supplementary
            self.mapQual = alignmentObj.mapping_quality
            self.supplAlignment = alignmentObj.get_tag('SA') if alignmentObj.has_tag('SA') else None
            self.refLen = alignmentObj.reference_length
            self.is_duplicate = alignmentObj.is_duplicate

    def readCoordinates(self):
        '''
        Compute read level alignment coordinates
 
         Output:
            1. begQuery: query start alignment position
            2. endQuery: query end alignment position
        '''
        orientation = '-' if self.reverse else '+'
        begQuery, endQuery = bamtools.alignment_interval_query(self.CIGAR, orientation)

        return begQuery, endQuery

    def clipped_interval_coordinates(self):
        '''
        Compute read level coordinates for the clipped piece of the read sequence

                                   bkp
                    ######READ######|********CLIPPED********
                                   beg                     end

                                          bkp
                    ********CLIPPED********|######READ######
                    beg                   end
        '''
        ### 1. Compute alignment coordinates at read level 
        alignmentBeg, alignmentEnd = self.readCoordinates()

        ### 2. Compute read length
        cigar = Cigar(self.CIGAR)
        readLen = len(cigar)

        ### 3. Obtain clipped interval coordinates 
        # a) >>>>>>>>>>******CLIPPED******
        if (not self.reverse) and (self.clippedSide == 'right'):
            clippedBeg = alignmentEnd
            clippedEnd = readLen

        # b) ******CLIPPED******<<<<<<<<<<<<<<<
        elif (self.reverse) and (self.clippedSide == 'left'):
            clippedBeg = alignmentEnd
            clippedEnd = readLen

        # c) ******CLIPPED******>>>>>>>>>>
        elif (not self.reverse) and (self.clippedSide == 'left'):
            clippedBeg = 0
            clippedEnd = alignmentBeg

        # d) <<<<<<<<<<<<<<<******CLIPPED******
        else:
            clippedBeg = 0
            clippedEnd = alignmentBeg

        return clippedBeg, clippedEnd


    def parse_supplAlignments_field(self):
        '''
        Parse supplementary alignment optional field. Create supplementary alignment objects and return them organized into a dictionary 
        based on the reference where they do align
    
        SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;) Other canonical alignments in a chimeric alignment, formatted as a semicolon-delimited list. 
        Each element in the list represents a part of the chimeric alignment. Conventionally, at a supplementary line, the first element points to 
        the primary line. Strand is either ‘+’ or ‘-’, indicating forward/reverse strand, corresponding to FLAG bit 0x10. Pos is a 1-based coordinate.
        
        Output: 
            1. supplAlignmentsDict: dictionary containing lists of supplementary alignment objects organized by reference
        '''
        
        ## Return empty list if no supplementary alignment available
        if self.supplAlignment is None:
            return {}

        ## Go ahead if supplementary alignments available
        supplAlignmentsDict = {}
        
        # For each supplementary alignment
        for supplAlignment in self.supplAlignment.split(';')[:-1]:

            # Extract info
            ref, beg, strand, CIGAR, mapQ, NM = supplAlignment.split(',')
            
            # Compute suppl. alignment end from CIGAR string
            alignmentLen = bamtools.alignment_length_cigar(CIGAR)
            end = int(beg) + alignmentLen

            # Create suppl. alignment object
            supplObject = SUPPLEMENTARY(ref, beg, end, strand, CIGAR, mapQ, NM, self.readName)

            # Initialize ref if necessary
            if supplObject.ref not in supplAlignmentsDict:
                supplAlignmentsDict[supplObject.ref] = []
            
            # Add object to the dictionary
            supplAlignmentsDict[supplObject.ref].append(supplObject)
        
        return supplAlignmentsDict


    def search4candidateCompl(self, supplAlignments):
        '''
        Search for supplementary alignments candidate to be complementary to the clipping event. 

        A suppl. alignment is candidate if the aligned piece of the read corresponds to the 
        clipped read sequence interval

                                   bkp
                    ######READ######|********CLIPPED********
                                     ####SUPPL####           * Candidate
                                                ##SUPPL##    * Candidate
                     ###SUPPL###                             * No candidate
        ####SUPPL####                                        * No candidate

        Input:
            1. supplAlignments: list of supplementary alignments
        
        Output:
            2. filteredSupplAlignments: list of filtered supplementary alignments
        '''

        ## 1. Compute clipped read sequence interval coordinates
        clippedBeg, clippedEnd = self.clipped_interval_coordinates()

        ## 2. Select only those supplementary alignments corresponding to the clipped piece of the read
        filteredSupplAlignments = []

        for alignment in supplAlignments:

            ## Compute supplementary alignment read level coordinates
            supplAlignmentBeg, supplAlignmentEnd = alignment.readCoordinates()

            ## Assess overlap 
            overlap, overlapLen = gRanges.overlap(supplAlignmentBeg, supplAlignmentEnd, clippedBeg, clippedEnd)

            ## Pick alignment if overlap is found
            if overlap:
                filteredSupplAlignments.append(alignment)
        
        return filteredSupplAlignments


    def search4complAlignments(self, supplAlignments):
        '''
        Search for supplementary alignments complementary to the clipping event. 
        
        A complementary alignment span a piece of the clipped read sequence, so contributes to explain 
        the full conformation of the read alignment. 

        Input:
            1. supplAlignments: list of supplementary alignments
        
        Output:
            1. complementaryAlignments: list of complementary alignments
        '''
        ## 1. Select candidate complementary alignments
        filteredSupplAlignments = self.search4candidateCompl(supplAlignments)

        ## 2. Initialize PAF
        PAF = formats.PAF()

        ## 3. Add clipping alignment to PAF
        # Create PAF alignment
        clippingBeg, clippingEnd = self.readCoordinates()
        strand = '-' if self.reverse else '+'
        fields = ['clipping', 0, clippingBeg, clippingEnd, strand, self.ref, 0, self.beg, self.end, 0, 0, self.mapQual]
        line = formats.PAF_line(fields)

        # Add alignment
        PAF.alignments.append(line)

        ## 4. Add supplementary alignments to PAF 
        supplAlignmentsDict = {}

        # For each supplementary alignment:
        for index, supplAlignment in enumerate(filteredSupplAlignments):
        
            # Create PAF alignment
            seqName = 'supplementary_' + str(index)
            supplBeg, supplEnd = supplAlignment.readCoordinates()
            fields = [seqName, 0, supplBeg, supplEnd, supplAlignment.orientation, supplAlignment.ref, 0, supplAlignment.beg, supplAlignment.end, 0, 0, supplAlignment.mapQ]
            line = formats.PAF_line(fields)

            # Add PAF alignment to the PAF
            PAF.alignments.append(line)

            # Add suppl. alignment to the dictionary
            supplAlignmentsDict[seqName] = supplAlignment

        ## 5. Search for complementariety
        chain = PAF.chain(10000, 20) 

        ## 6. Collect sorted complementary entries from chain
        # a) Anchor to the left of the compl. alignments 
        # ----anchor---->|---complementary_0--->|---complementary_1--->
        if chain.alignments[0].qName == 'clipping':
            anchor = chain.alignments[0]
            anchorSide = 'left'
            complementaryEntriesPAF = chain.alignments[1:]
        
        # b) Anchor to the right of the compl. alignments 
        # <---complementary_1---|<---complementary_0---|<----anchor----
        elif chain.alignments[-1].qName == 'clipping':
            anchor = chain.alignments[-1]
            anchorSide = 'right'
            complementaryEntriesPAF = reversed(chain.alignments[:-1])

        # c) No clipping in the chain
        else:
            complementaryEntriesPAF = []

        ## 7. Generate list of complementary alignments
        ## For each complementary alignment, set index specifying the order with respect the clipped piece of sequence        
        # Index 0 for suppl. alignment adjacent to the clipped piece of read squence
        # Index 1 for the next one, etc... 
        complementaryAlignments = []

        for index, entry in enumerate(complementaryEntriesPAF):
            
            ## Add index and anchor relative position to complementary alignment object 
            complementary = supplAlignmentsDict[entry.qName]
            complementary.readIndex = index
            complementary.anchorSide = anchorSide

            ## For complementary alignment adjacent to anchor (index == 0), 
            # collect potential inserted sequence at the BND junction
            if complementary.readIndex == 0:
            
                ## Compute distance between both bkp at read level
                # a) ----anchor---->|---complementary_0--->
                #         anchor_end|entry_beg                 
                if anchorSide == 'left':
                    bkpDist = entry.qBeg - anchor.qEnd

                # b) <---complementary_0---|<----anchor----
                #                 entry_end|anchor_beg
                else:
                    bkpDist = anchor.qBeg - entry.qEnd

                ## Retrieve inserted sequence if dist >= threshold (50)
                if bkpDist >= 50:
                    
                    complementary.insertSize = bkpDist
                    
                    # a) Right clipping (###anchor###----clipped----)
                    if self.clippedSide == 'right':
                        insertBeg = self.readBkp
                        insertEnd = self.readBkp + bkpDist

                    # b) Left clipping (----clipped----###anchor###)
                    else:
                        insertBeg = self.readBkp - bkpDist
                        insertEnd = self.readBkp                         

                    complementary.insertSeq = self.readSeq[insertBeg : insertEnd]
                 
            ## Add complementary alignment to the list 
            complementaryAlignments.append(complementary)
                
        return complementaryAlignments

class SUPPLEMENTARY():
    '''
    Supplementary alignment class
    '''
    number = 0 # Number of instances
    
    def __init__(self, ref, beg, end, orientation, CIGAR, mapQ, NM, readName):
        SUPPLEMENTARY.number += 1 # Update instances counter
        self.id = 'SUPPLEMENTARY_' + str(SUPPLEMENTARY.number)
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.orientation = orientation
        self.CIGAR = CIGAR
        self.mapQ = mapQ
        self.NM = NM
        self.readName = readName
        self.readIndex = None
        self.anchorSide = None
        self.insertSize = None
        self.insertSeq = None

    def readCoordinates(self):
        '''
        Compute read level alignment coordinates

        Output:
            1. begQuery: query start alignment position
            2. endQuery: query end alignment position
        '''
        begQuery, endQuery = bamtools.alignment_interval_query(self.CIGAR, self.orientation)
        return begQuery, endQuery
        
class DISCORDANT():
    '''
    Discordant class
    '''
    number = 0 # Number of instances
    
    def __init__(self, ref, beg, end, orientation, pair, readName, alignmentObj, sample, identity, duplicate):
        DISCORDANT.number += 1 # Update instances counter
        self.id = 'DISCORDANT_' + str(DISCORDANT.number)
        self.type = 'DISCORDANT'
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.orientation = orientation
        self.pair = str(pair)
        self.readName = readName 
        self.sample = sample
        self.is_duplicate = duplicate
        self.clusterId = None
        self.mapQual = alignmentObj.mapq
        self.cigarTuples = alignmentObj.cigartuples
        if identity == None:
            self.identity = None
            self.specificIdentity = None
        else:
            # As viral identity is a list, that can have more than one element, check number of elements of identity
            if len(identity) == 1:
                identity = ''.join(identity)
                self.identity = identity.split("|")[0]
                try:
                    self.specificIdentity = identity.split("|")[1]
                except IndexError:
                    self.specificIdentity = None

            # If length of identity is higher than 1, make lists of identity and specific identity.
            elif len(identity) > 1:
                genIdent = []
                specificIdent = []
                for ident in identity:
                    genIdent.append(ident.split('|')[0])
                    try:
                        specificIdent.append(ident.split('|')[1])
                    except IndexError:
                        pass

                # Uniq lists and sort them.
                identity = list(set(genIdent))
                identity.sort()
                if len(identity) == 1:
                    identity = ''.join(identity)
                else:
                    identity = '_'.join(identity)
                
                if specificIdent != []:
                    specificIdentity = list(set(specificIdent))
                    specificIdentity.sort()
                    if len(specificIdentity) == 1:
                        specificIdentity = ''.join(specificIdentity)
                    else:
                        specificIdentity = '_'.join(specificIdentity)
                else:
                    specificIdentity = None
                
                self.identity = identity
                self.specificIdentity = specificIdentity
        
        ## Mate info
        self.mateSeq = None

        if alignmentObj is None:
            self.mateRef = None
            self.mateStart = None           
        else:
            self.mateRef = alignmentObj.next_reference_name
            self.mateStart = alignmentObj.next_reference_start

    def fullReadName(self):
        '''
        Return the supporting read name including mate info
        '''
        fullReadName = self.readName + '/' + self.pair

        return fullReadName
    
    def fullReadName_mate(self):
        '''
        Return the supporting read name including mate info
        '''

        matePair = '2' if self.pair == '1' else '1'
        fullReadName = self.readName + '/' + matePair

        return fullReadName    

