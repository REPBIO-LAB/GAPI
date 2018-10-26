'''
Module 'variants' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External

# Internal


#############
## CLASSES ##
#############

## PRIMARY CLASSES ##
# Primary classes are the most basic types of class. They constitute the blocks for building other classes
class CLIPPING():
    '''
    Clipping class
    '''
    def __init__(self, alignmentObj, clippedSide, sample):
        self.type = 'CLIPPING' 
        self.ref = str(alignmentObj.reference_name)
        self.pos = int(alignmentObj.reference_start if clippedSide == 'left' else alignmentObj.reference_end)
        self.length = None
        self.clippedSide = clippedSide
        self.sample = sample

        ## Set supporting read id
        mate = '/1' if alignmentObj.is_read1 else '/2'
        self.readId = alignmentObj.query_name + mate  

        ## call functions
        self.clippingType = ''
        self.determine_clippingType(alignmentObj, clippedSide)

    def determine_clippingType(self, alignmentObj, clippedSide):
        '''
        Determine if soft or hard clipping

        Input: 
            1. alignmentObj: pysam read alignment object instance
            2. clippedSide: clipped side relative to the read (left or right)

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
    def __init__(self, ref, pos, length, seq, readId, sample):
        self.type = 'INS' 
        self.ref = str(ref)
        self.pos = int(pos)
        self.length = int(length)
        self.seq = seq
        self.readId = readId   
        self.sample = sample 


class DEL():
    '''
    Short deletion class. Deletion completely spanned by the read sequence
    '''
    def __init__(self, ref, beg, end, length, readId, sample):
        self.type = 'DEL' 
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
        self.length = int(length)
        self.readId = readId   
        self.sample = sample 


## SECONDARY CLASSES ##
class cluster():
    '''
    Events cluster class
    '''
    def __init__(self, events):
        '''

        '''

        # Define cluster chromosome, begin and end position
        self.ref = events[0].ref
        self.pos = events[0].pos
        self.end = events[-1].pos
        self.events = events
    
    def add(self, new_events, side):
        '''
      
        '''
        if side == 'left':
            self.beg = new_events[0].pos
            self.events = new_events + self.events 
        else:
            self.end = new_events[-1].pos
            self.events = self.events + new_events 
        




