'''
Module 'events' - Contains classes for dealing with structural variation events at single read level
'''

## DEPENDENCIES ##
# External

# Internal

###############
## FUNCTIONS ##
###############

#############
## CLASSES ##
#############

class CLIPPING():
    '''
    Clipping class
    '''
    number = 0 # Number of instances

    def __init__(self, alignmentObj, clippedSide, sample):
        CLIPPING.number += 1 # Update instances counter
        self.id = CLIPPING.number
        self.type = 'CLIPPING'
        self.ref = str(alignmentObj.reference_name)
        self.beg, self.end = (alignmentObj.reference_start, alignmentObj.reference_start) if clippedSide == 'left' else (alignmentObj.reference_end, alignmentObj.reference_end)
        self.clippedSide = clippedSide
        self.readId = alignmentObj.query_name
        self.sample = sample
        self.clusterId = None

        ## call functions
        self.clippingType = ''
        self.determine_clippingType(alignmentObj, clippedSide)

    def determine_clippingType(self, alignmentObj, clippedSide):
        '''
        Determine if soft or hard clipping

        Input:
            1. alignmentObj: pysam read alignment object instance
            2. clippedSide: Clipped side relative to the read (left or right)

        Output:
            - Update 'clippingType' class attribute
        '''
        ### Extract operation
        ## a) Clipping at the begin
        if clippedSide == 'left':
            operation =  alignmentObj.cigartuples[0][0]

        ## b) Clipping at the end
        else:
            operation = alignmentObj.cigartuples[-1][0]

        ### Define if soft or had clipping
        self.clippingType = 'soft' if (operation == 4) else 'hard'


class INS():
    '''
    Short insertion class. Insertion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, seq, readId, sample):
        INS.number += 1 # Update instances counter
        self.id = INS.number
        self.type = 'INS'
        self.ref = str(ref)
        self.beg = int(beg) # beg==end. 0-based leftmost coordinate of the insertion point
        self.end = int(end)
        self.length = int(length)
        self.seq = seq
        self.readId = readId
        self.sample = sample
        self.clusterId = None


class DEL():
    '''
    Short deletion class. Deletion completely spanned by the read sequence
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, length, readId, sample):
        DEL.number += 1 # Update instances counter
        self.id = DEL.number
        self.type = 'DEL'
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.length = int(length)
        self.readId = readId
        self.sample = sample
        self.clusterId = None


