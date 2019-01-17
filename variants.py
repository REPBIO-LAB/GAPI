'''
Module 'variants' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External
import sys
import numpy as np
import collections 

# Internal
import log
import formats
import unix
import consensus
import alignment
import bamtools 

###############
## FUNCTIONS ##
###############

def createCluster(events, eventsType):
    '''
    Function to create a cluster object instance

    Input:
        1. events: list of events that will compose the cluster
        2. eventsType: type of events to cluster (INS: insertion; DEL: deletion; CLIPPING: clipping)

    Output:
        1. cluster: cluster object instance
    '''
    cluster = ''

    ref = events[0].ref
    beg =  min([event.beg for event in events])
    end = max([event.end for event in events])

    ## a) Create INS cluster
    if (eventsType == 'INS'):
        cluster = INS_cluster(ref, beg, end, events, eventsType)

    ## b) Create DEL cluster
    elif (eventsType == 'DEL'):
        cluster = DEL_cluster(ref, beg, end, events, eventsType)

    ## c) Create CLIPPING cluster
    elif (eventsType == 'CLIPPING') or (eventsType == 'LEFT-CLIPPING') or (eventsType == 'RIGHT-CLIPPING'):
        cluster = CLIPPING_cluster(ref, beg, end, events, eventsType)

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

def polishClusters(clusters, clusterType):
    '''
    Function to polish a set of cluster objects. It does not produce any output just modify cluster objects through the polishing procedure

    Input:
        1. clusters: bin database containing a set of cluster objects
        2. clusterType: type of cluster (INS-CLUSTER: insertion; DEL-CLUSTER: deletion; LEFT-CLIPPING-CLUSTER: left clipping; RIGHT-CLIPPING-CLUSTER: right clipping)
    '''
    ## For each cluster
    for cluster in clusters.collect(clusterType):
        
        ## Polish
        cluster.polish()
        
def consensusClusters(clusters, clusterType, FASTQ, reference, confDict, rootDir):
    '''
    Function to create a consensus sequence for each cluster. 

    Input:
        1. clusters: bin database containing a set of cluster objects
        2. clusterType: type of cluster (INS-CLUSTER: insertion; DEL-CLUSTER: deletion; LEFT-CLIPPING-CLUSTER: left clipping; RIGHT-CLIPPING-CLUSTER: right clipping)
        3. FASTQ: FASTQ object containing all the reads to supporting SV clusters
        4. reference: path to reference genome in fasta format
        5. confDict: configuration dictionary (Complete...)
        6: rootDir: root directory to write files and directories
    '''

    outDir = rootDir + '/Consensus/'
    unix.mkdir(outDir)

    ## 1. Create consensus sequence for each cluster
    for cluster in clusters.collect(clusterType):
        
        ## Create consensus
        cluster.consensus(FASTQ, reference, confDict, outDir)

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


## SECONDARY CLASSES ##
# Classes composed by primary classes
class cluster():
    '''
    Events cluster class. A cluster is composed by a set of events, all of them of the same type.
    Each event is supported by a single read. One cluster can completely represent a single structural
    variation event or partially if multiple clusters are required (see 'metaCluster' class)
    '''
    number = 0 # Number of instances

    def __init__(self, ref, beg, end, events, clusterType):
        '''
        '''
        cluster.number += 1 # Update instances counter
        self.id = cluster.number

        # Define cluster chromosome, begin and end position
        self.ref = ref
        self.beg = beg
        self.end = end
        self.events = self.setClusterId(events) # Asign the cluster id to the events
        self.clusterType = clusterType

        # Cluster metrics
        self.nbOutliers = 0

        # Consensus sequences
        self.consensus_FASTA = None
        self.consensusIns = None
        self.consensusLen = None

    def add(self, newEvents, side):
        '''
        Incorporate events into the cluster.

        Input:
            1. newEvents: List of events. List have to be sorted in increasingly order if side specified.
            2. side: Add events to the 'left', to the 'right' of the cluster or add events and resort by begin position (if 'None')
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

    def supportingReads(self):
        '''
        Return list of cluster supporting reads
        '''
        readIds = [event.readId for event in self.events]

        return readIds

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
        nbTumour = 0
        nbNormal = 0

        for event in self.events:
            # a) Event identified in the TUMOUR sample
            if event.sample == "TUMOUR":
                nbTumour += 1
            
            # b) Event identified in the matched NORMAL sample
            elif event.sample == "NORMAL":
                nbNormal += 1
            
            # c) SINGLE sample mode
            else:
                nbTumour = None
                nbNormal = None
                break

        nbTotal = len(self.events)

        return nbTotal, nbTumour, nbNormal

    def meanLen(self):
        '''
        Compute the mean length for the events composing the cluster and the mean coefficient of variation
        '''
        lengths = []

        # Make a list with lengths of all events of the same cluster. 
        for event in self.events:

            # The event has the attribute length
            if hasattr(event, 'length'):
                lengths.append(event.length)

        # a) Length values available
        if lengths:

            ## Compute mean insertion length and standard deviation
            mean = np.mean(lengths)
            std = np.std(lengths)
            
            ## compute coefficient of variation (CV)
            # CV = std / mean * 100
            # std: standard deviation
            # mean: mean
            cv = std / mean * 100

        # b) Length not available
        else:
            mean = None
            std = None
            cv = None

        return mean, std, cv

    def polish(self):
        '''
        Apply successive rounds of polishing to the cluster by removing events whose length deviates from the cluster average
        '''

        # Check if length attribute is available for each event
        lengthsBool = [hasattr(event, 'length') for event in self.events]

        # A) Attemp polishing if length attribute available for all the events 
        if False not in lengthsBool:

            ## 1. Compute length metrics for the initial cluster 
            initialNbEvents = self.nbEvents()[0]
            mean, std, cv = self.meanLen()

            ## 2. Apply successive rounds of polishing while cv > threshold 
            while cv > 15: 
             
                ## 2.1 Set length cutoffs 
                cutOff = std * 1 
                lowerBound, upperBound = mean - cutOff, mean + cutOff

                ## 2.2 Generate list of events composing the polished cluster
                eventsAfterPolish = []

                # Evaluate for each event if it´s an ourlier or not
                for event in self.events:

                    # a) No outlier. Event length within boundaries -> include event into polished cluster   
                    if (event.length >= lowerBound) and (event.length <= upperBound):
                        eventsAfterPolish.append(event)

                    # b) Outlier. Event lenght outside boundaries

                ## 2.3 Recompute length metrics for polished cluster
                ## Prior polishing
                cvPrior = cv
                eventsPrior = self.events

                ## After polishing round
                self.events = eventsAfterPolish
                mean, std, cv = self.meanLen()

                ## 2.4 Stop polishing and use previous cluster state if current polishing round does not reduce the cv
                if cv >= cvPrior:
                    self.events = eventsPrior # Use previous cluster state
                    break

            ## 3. Compute the number of outliers
            finalNbEvents = self.nbEvents()[0]
            self.nbOutliers = initialNbEvents - finalNbEvents

            ## 4. Recompute cluster begin and end after polishing
            self.beg = self.events[0].beg
            self.end = self.events[-1].end

        # B) Don´t attemp polishing if length attribute not available for some of the events 
        
    def consensus(self, FASTQ, reference, confDict, outDir): 
        '''
        Use all the cluster supporting reads to (1) generate a consensus sequence for the cluster and (2) perform local realignment of the consensus 
        to obtain more accurate SV breakpoints and inserted sequence (only for INS)

        Update cluster with all this information
        
        Input:
            1. FASTQ: FASTQ object containing all the SV supporting reads  
            2. reference: Path to the reference genome in fasta format. An index of the reference generated with samtools faidx must be located in the same directory
            3. confDict: Configuration dictionary 
            4. outDir: Output directory
        '''
        ## 1. Create FASTQ object containing cluster supporting reads
        FASTQ_cluster = formats.FASTQ()
        readIds = self.supportingReads()

        # Collect from the FASTQ containing all the reads only those supporting the SV cluster 
        for readId in readIds:
            FASTQ_entry = FASTQ.fastqDict[readId]
            FASTQ_cluster.add(FASTQ_entry)

        ## 2. Generate consensus sequence for the cluster
        self.consensus_FASTA = consensus.racon(FASTQ_cluster, outDir)

        # Exit if consensus sequence could not be generated
        if self.consensus_FASTA == None:
            return

        ## 3. Realign consensus sequence into the SV event genomic region
        ## 3.1 Write consensus sequence into fasta file
        consensusFile = outDir + '/consensus.fasta'
        self.consensus_FASTA.write(consensusFile)

        ## 3.2 Define SV event surrounding region
        offset = 500
        targetBeg = self.beg - offset
        targetEnd = self.end + offset
        targetLen = targetEnd - targetBeg
        targetInterval = self.ref + ':' + str(targetBeg) + '-' + str(targetEnd)

        ## 3.3 Do realignment
        BAM = alignment.targeted_alignment_minimap2(consensusFile, targetInterval, reference, outDir)

        ## 4. Extract consensus SV event from consensus sequence realignment
        confDict["targetSV"] = self.clusterType
        INS_events, DEL_events, CLIPPING_left_events, CLIPPING_right_events, FASTQ = bamtools.collectSV(targetInterval, 0, targetLen, BAM, confDict, None)
               
        ## 5. Redefine cluster properties based on consensus SV event  
        # A) Insertion cluster
        if self.clusterType == 'INS':
            
            # a) Single event 
            if len(INS_events) == 1:
                INS = INS_events[0]
                self.beg = INS.beg + targetBeg # Map alignments into chromosomal coordinates
                self.end = INS.end + targetBeg
                self.consensusLen = INS.length
                self.consensusIns = INS.seq

            # b) Multiple possible events. Select most likely (TO DO) 
            elif len(INS_events) > 1:
                print("WARNING. Multiple possible INS events")

            # c) No identified events
             
        # B) Deletion cluster
        elif self.clusterType == 'DEL':

            # a) Single event 
            if len(DEL_events) == 1:
                DEL = DEL_events[0]
                self.beg = DEL.beg + targetBeg # Map alignments into chromosomal coordinates
                self.end = DEL.end + targetBeg
                self.consensusLen = self.end - self.beg 

            # b) Multiple possible events. Select most likely (TO DO)
            elif len(DEL_events) > 1:
                print("WARNING. Multiple possible DEL events")

            # c) No identified events

        # C) Clipped cluster 



class INS_cluster(cluster):
    '''
    Insertion (INS) cluster subclass
    '''
    def __init__(self, ref, beg, end, events, clusterType):

        cluster.__init__(self, ref, beg, end, events, clusterType)

class DEL_cluster(cluster):
    '''
    Deletion (DEL) cluster subclass
    '''
    def __init__(self, ref, beg, end, events, clusterType):

        cluster.__init__(self, ref, beg, end, events, clusterType)

class CLIPPING_cluster(cluster):
    '''
    Clipping cluster subclass
    '''
    def __init__(self, ref, beg, end, events, clusterType):

        cluster.__init__(self, ref, beg, end, events, clusterType)