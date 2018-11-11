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

def createCluster(events, clusterType):
    '''
    Function to create a cluster object instance

    Input:
        1. events: list of events that will compose the cluster
        2. clusterType: type of cluster (INS: insertion; DEL: deletion; CLIPPING: clipping)
 
    Output:
        1. cluster: cluster object instance
    '''
    cluster = ''
    
    ref = events[0].ref
    beg =  min([event.beg for event in events])
    end = max([event.end for event in events]) 

    ## a) Create INS cluster
    if (clusterType == 'INS'):
        cluster = INS_cluster(ref, beg, end, events)

    ## b) Create DEL cluster
    elif (clusterType == 'DEL'):
        cluster = DEL_cluster(ref, beg, end, events)

    ## c) Create CLIPPING cluster
    elif (clusterType == 'CLIPPING') or (clusterType == 'LEFT-CLIPPING') or (clusterType == 'RIGHT-CLIPPING'):
        cluster = CLIPPING_cluster(ref, beg, end, events)

    ## d) Unexpected cluster type
    else:
        log.info('Error at \'createCluster\'. Unexpected cluster type')
        sys.exit(1)

    return cluster


def mergeClusters(clusters, clusterType):
    '''
    Merge a set of clusters into a single cluster instance

    Input:
        1. clusters: list of clusters that will be merged
        2. clusterType: type of cluster (INS: insertion; DEL: deletion; CLIPPING: clipping)
 
    Output:
        1. cluster: merged cluster instance
    '''
    events = []
    for cluster in clusters:
        events = events + cluster.events

    mergedCluster = createCluster(events, clusterType)

    return mergedCluster



#############
## CLASSES ##
#############

## PRIMARY CLASSES ##
# Primary classes are the most basic types of class. They constitute the blocks for building other classes
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
        self.sample = sample
        self.clusterId = None

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
    

## SECONDARY CLASSES ##
# Classes composed by primary classes
class cluster():
    '''
    Events cluster class. A cluster is composed by a set of events, all of them of the same type. 
    Each event is supported by a single read. One cluster can completely represent a single structural 
    variation event or partially if multiple clusters are required (see 'metaCluster' class)
    '''
    number = 0 # Number of instances 

    def __init__(self, ref, beg, end, events):
        '''
        '''
        cluster.number += 1 # Update instances counter
        self.id = cluster.number

        # Define cluster chromosome, begin and end position
        self.ref = ref
        self.beg = beg
        self.end = end
        self.events = self.setClusterId(events) # Asign the cluster id to the events
    
    def add(self, newEvents, side):
        '''
        Incorporate events into the cluster. 

        Input: 
            1. newEvents: List of events. List have to be sorted in increasingly order if side specified. 
            2. side: add events to the 'left', to the 'right' of the cluster or add events and resort by begin position (if 'None')
        '''
        # Asign the cluster id to the events
        newEvents = self.setClusterId(newEvents)

        # a) Add to the left side of the cluster
        if side == 'left':
            self.events = newEvents + self.events 
            self.beg = self.events[0].beg # Update begin

        # b) Add to the right side of the cluster
        elif side == 'right':
            self.events = self.events + newEvents 
            self.end = self.events[-1].end # Update end
        
        # c) Add events to the cluster and resort by begin position
        else:
            self.events = self.events + newEvents
            self.events.sort(key=lambda x: x.beg, reverse=True) # Resort
            self.beg = self.events[0].beg # Update begin
            self.end = max([event.end for event in self.events]) # Update end

    def setClusterId(self, events):
        '''
        Set cluster identifier for a set of events

        Input: 
            1. events: List of events.
        Output:
            1. events: List of events in the same sorting order with the cluster identifier asigned
        '''
        for event in events:
            event.clusterId = self.id
        
        return events
    
    def nbEvents(self):
        '''
        Return the number of events composing the cluster
        '''
        return len(self.events)


            
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
