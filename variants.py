'''
Module 'variants' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External
import  sys

# Internal
import log

###############
## FUNCTIONS ##
###############

def makeCluster(events, clusterType):
    '''
    Function to create a cluster object instance

    Input:
        1. events: list of events that will compose the cluster
        2. clusterType: type of cluster (INS: insertion; DEL: deletion; CLIPPING: clipping)
 
    Output:
        1. cluster: cluster object instance
    '''
    cluster = ''
    
    ## a) Create INS cluster
    if (clusterType == 'INS'):
        
        ref = events[0].ref
        beg = events[0].pos
        end = events[-1].pos
        cluster = INS_cluster(ref, beg, end, events)

    ## b) Create DEL cluster
    elif (clusterType == 'DEL'):

        ref = events[0].ref
        beg = events[0].beg
        end = events[-1].end
        cluster = DEL_cluster(ref, beg, end, events)

    ## c) Create CLIPPING cluster
    elif (clusterType == 'CLIPPING'):

        ref = events[0].ref
        beg = events[0].pos
        end = events[-1].pos
        cluster = CLIPPING_cluster(ref, beg, end, events)

    ## d) Unexpected cluster type
    else:
        log.info('Error at \'makeCluster\'. Unexpected cluster type')
        sys.exit(1)

    return cluster

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
    def __init__(self, ref, beg, end, events):
        '''
        '''

        # Define cluster chromosome, begin and end position
        self.ref = ref
        self.beg = beg
        self.end = end
        self.events = events
    
    def add(self, newEvents, side):
        '''
        Incorporate events into the cluster. 

        Input: 
            1. newEvents: sorted list of events (from lower to upper position)
            2. side: add events to the 'left' or to the 'right' of the cluster

        Output:
            - Update 'clippingType' class attribute 
        '''
        if side == 'left':
            self.beg = newEvents[0].pos
            self.events = newEvents + self.events 
        else:
            self.end = newEvents[-1].pos
            self.events = self.events + newEvents 
        


class INS_cluster(cluster):
    '''
    Insertion (INS) cluster subclass
    '''
    def __init__(self, ref, beg, end, events):
        cluster.__init__(self, ref, beg, end, events)

class DEL_cluster(cluster):
    '''
    Deletion (DEL) cluster subclass
    '''
    def __init__(self, ref, beg, end, events):
        cluster.__init__(self, ref, beg, end, events)

class CLIPPING_cluster(cluster):
    '''
    Clipping cluster subclass
    '''
    def __init__(self, ref, beg, end, events):
        cluster.__init__(self, ref, beg, end, events)
