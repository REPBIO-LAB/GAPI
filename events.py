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

#############
## CLASSES ##
#############

class CLIPPING():
    '''
    Clipping class
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, clippedSide, readPos, alignmentObj, sample):
        '''
        '''
        CLIPPING.number += 1 # Update instances counter
        self.id = CLIPPING.number
        self.type = 'CLIPPING'
        self.ref = str(ref)
        self.beg = int(beg) # beg==end. 0-based CLIPPING breakpoint
        self.end = int(end)
        self.clippedSide = clippedSide
        self.readPos = int(readPos)
        self.alignmentObj = alignmentObj
        self.sample = sample

        ## Determine if soft or hard clipping
        self.clippingType = self.determine_clippingType(alignmentObj, clippedSide)

    def determine_clippingType(self, alignmentObj, clippedSide):
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

        ### Define if soft or had clipping
        clippingType = 'soft' if (operation == 4) else 'hard'

        return clippingType

class INS():
    '''
    Short insertion class. Insertion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, readPos, alignmentObj, sample):
        '''
        '''

        INS.number += 1 # Update instances counter
        self.id = INS.number
        self.type = 'INS'
        self.ref = str(ref)
        self.beg = int(beg) # beg==end. 0-based INS breakpoint
        self.end = int(end)
        self.length = int(length)
        self.readPos = int(readPos)
        self.alignmentObj = alignmentObj
        self.sample = sample

class DEL():
    '''
    Short deletion class. Deletion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, readPos, alignmentObj, sample):
        '''
        '''
        DEL.number += 1 # Update instances counter
        self.id = DEL.number
        self.type = 'DEL'
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.length = int(length)
        self.readPos = int(readPos)
        self.alignmentObj = alignmentObj
        self.sample = sample

