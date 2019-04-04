'''
Module 'events' - Contains classes for dealing with structural variation events at single read level
'''

## DEPENDENCIES ##
# External

# Internal

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
    '''
    ## Pick inserted fragment sequence + flanking read sequence
    begPos = readPos - overhang
    begPos = begPos if begPos >= 0 else 0 # set lower bound to 0
    endPos = readPos + length + overhang
    seq = readSeq[begPos:endPos]

    return seq

def pick_flanking_seq_DEL(readSeq, readPos, overhang):
    '''
    Pick sequence flanking the deletion breakpoints from DEL supporting read

    Input:
        1. readSeq: DEL supporting read
        2. readPos: read position where DEL start
        3. overhang: number of flanking base pairs around the DEL event to be collected from the supporting read sequence 

    Output:
        1. seq: piece of supporting read sequence spanning the DEL event
    '''
    ## Pick deletion flanking read sequence
    begPos = readPos - overhang
    begPos = begPos if begPos >= 0 else 0 # set lower bound to 0
    endPos = readPos + overhang
    seq = readSeq[begPos:endPos]

    return seq


def pick_flanking_seq_CLIPPING(readSeq, readPos, clippedSide, overhang):
    '''
    Pick sequence flanking the clipping breakpoint from CLIPPING supporting read

    Input:
        1. readSeq: CLIPPING supporting read
        2. readPos: read position where CLIPPING start
        3. clippedSide: clipping orientation (left or right)
        4. overhang: number of flanking base pairs around the CLIPPING breakpoint position to be collected from the supporting read sequence 

    Output:
        1. seq: piece of supporting read sequence spanning the CLIPPING breakpoint event
    '''

    ## Take into account clipping orientation to pick read sequence (TO DO) 
    # a) Left clipping (.........*readPos*###################) 
    #                   ------------------<-overhang-> *endPos
    if clippedSide == 'left':
        endPos = readPos + overhang
        seq = readSeq[:endPos]

    # b) Rigth clipping (###################*readPos*.........) 
    #                   begPos* <-overhang->------------------ 
    else:
        begPos = readPos - overhang
        begPos = begPos if begPos >= 0 else 0 # set lower bound to 0
        seq = readSeq[begPos:]

    return seq


#############
## CLASSES ##
#############

class INS():
    '''
    Short insertion class. Insertion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, readSeq, readName, sample):
        '''
        '''

        INS.number += 1 # Update instances counter
        self.id = INS.number
        self.type = 'INS'
        self.ref = str(ref)
        self.beg = int(beg) # beg==end. 0-based INS breakpoint
        self.end = int(end)
        self.length = int(length)
        self.readSeq = readSeq
        self.readName = readName
        self.sample = sample

class DEL():
    '''
    Short deletion class. Deletion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, readSeq, readName, sample):
        '''
        '''
        DEL.number += 1 # Update instances counter
        self.id = DEL.number
        self.type = 'DEL'
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.length = int(length)
        self.readSeq = readSeq
        self.readName = readName
        self.sample = sample
    
class CLIPPING():
    '''
    Clipping class
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, clippedSide, clippingType, readSeq, readName, sample):
        '''
        '''
        CLIPPING.number += 1 # Update instances counter
        self.id = CLIPPING.number
        self.type = 'CLIPPING'
        self.ref = str(ref)
        self.beg = int(beg) # beg==end. 0-based CLIPPING breakpoint
        self.end = int(end)
        self.clippedSide = clippedSide
        self.clippingType = clippingType
        self.readSeq = readSeq
        self.readName = readName
        self.sample = sample
