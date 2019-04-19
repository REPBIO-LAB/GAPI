'''
Module 'variants' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External
import sys
import numpy as np
import collections 
# [SR CHANGE]
import pysam

# Internal
import log
import formats
import unix
import clustering
import events
import structures
import alignment
import bamtools 
import repeats
import retrotransposons
import filters
## [SR CHANGE]
import bkp

###############
## FUNCTIONS ##
###############

def create_cluster(events, clusterType):
    '''
    Function to create a cluster object instance

    Input:
        1. events: list of events that will compose the cluster
        2. clusterType: type of cluster (META: metacluster; INS: insertion; DEL: deletion; CLIPPING: clipping; DISCORDANT: discordant paired-end)

    Output:
        1. cluster: cluster object instance
    '''
    cluster = ''

    ## a) Create META cluster
    if (clusterType == 'META'):
        cluster = META_cluster(events)

    ## b) Create INS cluster
    elif (clusterType == 'INS'):
        cluster = INS_cluster(events)

    ## c) Create DEL cluster
    elif (clusterType == 'DEL'):
        cluster = DEL_cluster(events)

    ## d) Create CLIPPING cluster
    elif (clusterType == 'CLIPPING') or (clusterType == 'LEFT-CLIPPING') or (clusterType == 'RIGHT-CLIPPING'):
        cluster = CLIPPING_cluster(events)

    # [SR CHANGE]
    ## e) Create DISCORDANT cluster
    elif 'DISCORDANT' in clusterType:
        cluster = DISCORDANT_cluster(events)

    ## f) Unexpected cluster type
    else:
        log.info('Error at \'create_cluster\'. Unexpected cluster type')
        sys.exit(1)

    return cluster

def merge_clusters(clusters, clusterType):
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

    mergedCluster = create_cluster(events, clusterType)

    return mergedCluster

## [SR CHANGE]
def create_discordantClusters(eventsBinDb, confDict):
    '''
    '''
    discordantClustersDict = {}

    if 'DISCORDANT' in confDict['targetSV']:
        eventTypes = eventsBinDb.collectEventTypes()
        for eventType in eventTypes:
            discordantClustersDict[eventType] = clustering.reciprocal_clustering(eventsBinDb, 1, confDict['minClusterSize'], eventType, 0, eventType)

    return discordantClustersDict

# [SR CHANGE]
def create_metaclusters(eventsBinDb, confDict, bam, normalBam, mode):
    '''
    
    Input:
        1. eventsBinDb: Data structure containing a set of events organized in genomic bins
        2. confDict: 
            * targetSV            -> list with types of SV events included in the bin database provided as input 
            * maxEventDist          -> Maximum distance bewteen two adjacent breakpoints for INS and CLIPPING clustering
            * minClusterSize  -> minimum number of reads composing a root cluster

    Output:
        1. metaclustersBinDb: bin database structure containing metaclusters
    ''' 
    ## 1. Abort metaclustering if smallest bin longer than maxEventDist ##
    binSize = eventsBinDb.binSizes[0]
    
    if (binSize > confDict['maxEventDist']):
        step = 'ERROR'
        msg = 'Smallest bin can not be longer than maxEventDist'
        log.step(step, msg)
        sys.exit()

    allMetaclusters = []
    clusterType = 'META'

    ## 2. Do metaclustering ##
    # 2.1) Create INS + CLIPPING metaclusters
    if all(event_type in confDict['targetSV'] for event_type in ['INS', 'CLIPPING']):
        eventTypes = ['INS', 'LEFT-CLIPPING', 'RIGHT-CLIPPING']
        metaclusters = clustering.distance_clustering(eventsBinDb, binSize, eventTypes, clusterType, confDict['maxEventDist'], confDict['minClusterSize'])         
        allMetaclusters = allMetaclusters + metaclusters

    else:

        # 2.2) Create INS metaclusters
        if 'INS' in confDict['targetSV']:
            eventTypes = ['INS']
            metaclusters = clustering.distance_clustering(eventsBinDb, binSize, eventTypes, clusterType, confDict['maxEventDist'], confDict['minClusterSize'])         
            allMetaclusters = allMetaclusters + metaclusters

        # 2.3) Create CLIPPING metaclusters
        if 'CLIPPING' in confDict['targetSV']:
            eventTypes = ['LEFT-CLIPPING', 'RIGHT-CLIPPING']
            metaclusters = clustering.distance_clustering(eventsBinDb, binSize, eventTypes, clusterType, confDict['maxEventDist'], confDict['minClusterSize'])         
            allMetaclusters = allMetaclusters + metaclusters

    # 2.4) Create DEL metaclusters 
    if 'DEL' in confDict['targetSV']:
        eventTypes = ['DEL']
        # TO DO  
        # metaclusters = clustering.reciprocal_clustering()   
        # allMetaclusters = allMetaclusters + metaclusters

    # 2.5) Create DISCORDANT metaclusters
    if 'DISCORDANT' in confDict['targetSV']:
        eventTypes = eventsBinDb.collectEventTypes()
        for eventType in eventTypes:
        # TO DO  
            metaclusters = clustering.reciprocal_clustering(eventsBinDb, 1, confDict['minClusterSize'], eventType, 0, 'META')
            allMetaclusters = allMetaclusters + metaclusters
        # metaclusters = clustering.reciprocal_clustering()
        # allMetaclusters = allMetaclusters + metaclusters


    for metacluster in allMetaclusters:
        # [SR CHANGE]:
        #print('METACLUSTER: ', metacluster, len(metacluster.events), [(clusterType, len(subcluster.events)) for clusterType, subcluster in metacluster.subclusters.items()])
        # [SR CHANGE]:
        print('METACLUSTER: ', str(metacluster) +' '+ str(len(metacluster.events)) +' '+ str(metacluster.ref) +' '+ str(metacluster.beg) +' '+ str(metacluster.end) +' '+ str(metacluster.intOrigin))
        # [SR CHANGE]:
        for event in metacluster.events:
                #print (str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' ' + str(event.identity) + ' ' + str(event.side))
                print (str(metacluster) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type))
        #CLIPPING_cluster = metacluster.supportingCLIPPING(1, confDict, bam, normalBam, mode)
        #bkp.clippingBkp(CLIPPING_cluster)

        if mode == 'PAIRED':
            metacluster.setIntOrigin()

        print('METACLUSTER ADDED: ', str(metacluster) +' '+ str(len(metacluster.events)) +' '+ str(metacluster.ref) +' '+ str(metacluster.beg) +' '+ str(metacluster.end) +' '+ str(metacluster.intOrigin))
        # [SR CHANGE]:
        for event in metacluster.events:
                if event.type == 'DISCORDANT':
                    print (str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' ' + str(event.identity) + ' ' + str(event.side) + ' ' + str(event.sample))
                else:
                    print (str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' None ' + str(event.clippedSide) + ' ' + str(event.sample))


    ## 3. Organize metaclusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    metaclustersDict = {}
    metaclustersDict['METACLUSTERS'] = allMetaclusters 

    metaclustersBinDb = structures.create_bin_database(eventsBinDb.ref, eventsBinDb.beg, eventsBinDb.end, metaclustersDict, binSizes)
    return metaclustersBinDb


def find_chimeric_alignments(clusterA, clusterB):
    '''
    Search for chimeric read alignments connecting two clipping clusters. Select one as representative if multiple are identified. 

    Input:
        1. clusterA: Clipping cluster object
        2. clusterB: Clipping cluster object

    Output:
        1. primary: clipping event for representative primary alignment
        2. supplementary: clipping event for representative supplementary alignment
        4. chimericSorted: list of tuples. Each tuple is composed by two clipping events corresponding to primary and supplementary alignments, respectively. 
    ''' 
    ### 1. Identify chimeric read alignments connecting both clipping clusters
    chimeric = []

    # For each clipping event composing cluster A        
    for clippingA in clusterA.events:

        # For each clipping event composing cluster B
        for clippingB in clusterB.events:

            # Clipping events from the same read
            if (clippingA.readName == clippingB.readName):
                        
                # a) Clipping A primary while clipping B supplementary
                if (clippingA.supplementary == False) and (clippingB.supplementary == True):
                    primary = clippingA
                    supplementary = clippingB
                    chimeric.append((primary, supplementary))

                # b) Left clipping supplementary while right clipping primary
                elif (clippingA.supplementary == True) and (clippingB.supplementary == False):
                    primary = clippingB
                    supplementary = clippingA
                    chimeric.append((primary, supplementary))

                # c) Both clippings supplementary (Discard! These cases are not informative as supplementary alignments are hardclipped so don´t allow to pick the inserted fragment)
                
    # Exit if not chimeric alignments found
    if not chimeric:
        return False, None, None, None

    ### 2. Select the clipping events with the longest supplementary alignment as representative
    chimericSorted = sorted(chimeric, key=lambda alignments: alignments[1].refLen, reverse=True)
    primary = chimericSorted[0][0]
    supplementary = chimericSorted[0][1]

    return True, primary, supplementary, chimericSorted

def find_insertion_at_clipping_bkp(primary, supplementary):
    '''
    Search for inserted DNA between the two clipping event breakpoints. 

    Input:
        1. primary: Clipping event object corresponding to the primary alignment
        2. supplementary: Clipping event object corresponding to the supplementary alignment

    Output:

        1. insert: DNA fragment inserted between clipping breakpoints
    ''' 
    ### Check if there is an unaligned piece of read sequence between the primary and supplementary alignment breakpoints
    ## A) Primary right clipped 
    # -: aligned; _: clipped            readBkpA
    # primary -----------------------------*______________________
    #                                       <<<<insert>>>>*------- supplementary
    #                                                  readBkpB == total read length -  length of the supplementary piece of sequence aligned
    if (primary.clippedSide == 'right'):
        readBkpA = primary.readBkp
        readBkpB = len(primary.readSeq) - len(supplementary.readSeq)

    ## B) Primary left clipped
    # -: aligned; _: clipped        readBkpB
    # primary        ______________________*-----------------------------
    # supplementary  -------*<<<<insert>>>>                            
    #               readBkpA == length of the supplementary piece of sequence aligned
    else:
        readBkpA = len(supplementary.readSeq)
        readBkpB = primary.readBkp             

    # Extract inserted sequence between clipping breakpoints
    insert = primary.readSeq[readBkpA:readBkpB]

    return insert


#############
## CLASSES ##
#############

class cluster():
    '''
    Events cluster class. A cluster is composed by a set of events.
    Each event is supported by a single read. One cluster can completely represent a single structural
    variation event or partially if multiple clusters are required (see 'metaCluster' sub class)
    '''
    number = 0 # Number of instances

    def __init__(self, events, clusterType):
        '''
        '''
        cluster.number += 1 # Update instances counter
        self.id = cluster.number

        # Define list of events composing the cluster and cluster type
        self.events = events
        self.clusterType = clusterType

        # Set cluster's reference, begin and end position
        self.ref, self.beg, self.end = self.coordinates() 

        # Cluster filtering
        self.filters = None
        self.nbOutliers = 0

    def sort(self):
        '''
        Sort events in increasing coordinates order
        '''
        self.events.sort(key=lambda event: event.beg)

    def coordinates(self):
        '''
        Compute cluster ref, beg and end coordinates. 
        
        Begin and end will correspond to the left and rightmost positions, respectively
        '''
        # Sort events from lower to upper beg coordinates
        self.sort()  

        # Define cluster coordinates 
        ref = self.events[0].ref
        beg = self.events[0].beg
        end = max([event.end for event in self.events])
        
        return ref, beg, end

    def add(self, events):
        '''
        Incorporate events into the cluster and redefine cluster beg and end
        positions accordingly

        Input:
            1. events: List of events to be added to the cluster
        '''
        # Add events to the cluster  
        self.events = self.events + events

        # Resort and redefine cluster begin and end coordinates
        self.ref, self.beg, self.end = self.coordinates() 
            
    def collect_reads(self):
        '''
        Create FASTA object containing cluster supporting reads.
        
        Output:
            1. FASTA: FASTA object containing cluster supporting reads
        '''
        ## Initiate FASTA object
        FASTA = formats.FASTA()

        ## Add reads supporting the events to the FASTA
        for event in self.events:

            FASTA.seqDict[event.readName] = event.readSeq

        return FASTA

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


class INS_cluster(cluster):
    '''
    Insertion (INS) cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'INS')

        # Insertion features
        self.status = None
        self.insType = None
        self.family = None 
        self.srcId = None
        self.percResolved = None
        self.strand = None
        self.hits = None
    
        # Inserted seq
        self.consensus_FASTA = None
        self.consensusLen = None
        self.isConsensus = None
        self.insertSeq = None


class DEL_cluster(cluster):
    '''
    Deletion (DEL) cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'DEL')


class CLIPPING_cluster(cluster):
    '''
    Clipping cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'CLIPPING')

## [SR CHANGE]
class DISCORDANT_cluster(cluster):
    '''
    Discordant cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'DISCORDANT')

class META_cluster():
    '''
    Meta cluster class
    '''
    number = 0 # Number of instances

    def __init__(self, events):
        '''
        '''
        META_cluster.number += 1 # Update instances counter
        self.id = META_cluster.number

        # Define list of events composing the cluster 
        self.events = events

        # Set cluster's reference, begin and end position
        self.ref, self.beg, self.end = self.coordinates() 

        # Cluster filtering
        self.filters = None

        # Organize events into subclusters
        self.subclusters = self.create_subclusters()

        # Tag germline or somatic
        self.intOrigin = None

    def sort(self):
        '''
        Sort events in increasing coordinates order
        '''
        self.events.sort(key=lambda event: event.beg)

    def coordinates(self):
        '''
        Compute cluster ref, beg and end coordinates. 
        
        Begin and end will correspond to the left and rightmost positions, respectively
        '''
        # Sort events from lower to upper beg coordinates
        self.sort()  

        # Define cluster coordinates 
        ref = self.events[0].ref
        beg = self.events[0].beg
        end = max([event.end for event in self.events])
        
        return ref, beg, end

    def create_subclusters(self):
        '''
        '''
        ## 1. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(self.events)

        ## 2. Create subclusters ##
        subclusters = {}

        for eventType, eventList in eventTypes.items():

            ## Create subcluster
            subcluster = create_cluster(eventList, eventType) 

            ## Add subcluster to the dict
            subclusters[eventType] = subcluster 

        return subclusters

    def add(self, eventsList):
        '''

        Input:
            1. events: List of events to be added to the metacluster
        '''
        ## 1. Add events to the cluster ##
        previous = self.events
        self.events = self.events + eventsList

        ## 2. Resort and redefine cluster begin and end coordinates ##
        self.ref, self.beg, self.end = self.coordinates()

        ## 3. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(eventsList)

        ## 4. Add events to the subclusters ##
        for eventType, eventList in eventTypes.items():
            
            # a) Create subcluster if not pre-existing one
            if eventType not in self.subclusters:
         
                ## Create subcluster
                subcluster = create_cluster(eventList, eventType) 
            
                ## Add subcluster to the dict
                self.subclusters[eventType] = subcluster 

            # b) Add events to pre-existing subcluster
            else:
                self.subclusters[eventType].add(eventList)
                
    def nbEvents(self):
        '''
        Return the number of events composing the metacluster. 
        
        Extend to return a as well a dictionary containing the number of events composing each subcluster
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

    def supportingCLIPPING(self, buffer, confDict, bam, normalBam, mode):
        '''
        Look for clipping reas within the region of the existing discordant clusters. 
        It doesn't return anything
        Take into account that discordant and clipped reads could be duplicated!!
        '''

        clippingEventsDict = {}
        clippingEventsDict['RIGHT-CLIPPING'] = []
        clippingEventsDict['LEFT-CLIPPING'] = []
        sample = None

        ## Define region
        if self.beg > buffer:
            binBeg = self.beg - buffer

        else:
            binBeg = self.beg
        
        # TODO check as above
        binEnd = self.end

        ref = self.ref

        if normalBam != None:
            sample = 'NORMAL'
            ## If there is normal bam:
            normalBamFile = pysam.AlignmentFile(normalBam, "rb")

            ## Extract alignments
            iterator = normalBamFile.fetch(ref, binBeg, binEnd)

            for alignmentObj in iterator:
                    try:
                        ## Collect clipping
                        # TODO hacer para paired tb
                        targetInterval = None
                        CLIPPING_left_alignmentObj, CLIPPING_right_alignmentObj = bamtools.collectCLIPPING(alignmentObj, confDict['minCLIPPINGlen'], targetInterval, confDict['overhang'], sample)
                        if CLIPPING_left_alignmentObj != None:
                            clippingEventsDict['LEFT-CLIPPING'].append(CLIPPING_left_alignmentObj)
                        if CLIPPING_right_alignmentObj != None:
                            clippingEventsDict['RIGHT-CLIPPING'].append(CLIPPING_right_alignmentObj)
                    except TypeError:
                        continue
            sample = 'TUMOUR'

        ## Open BAM file for reading
        bamFile = pysam.AlignmentFile(bam, "rb")

        ## Extract alignments
        iterator = bamFile.fetch(ref, binBeg, binEnd)

        for alignmentObj in iterator:
                try:
                    ## Collect clipping
                    # TODO hacer para paired tb
                    targetInterval = None
                    CLIPPING_left_alignmentObj, CLIPPING_right_alignmentObj = bamtools.collectCLIPPING(alignmentObj, confDict['minCLIPPINGlen'], targetInterval, confDict['overhang'], sample)
                    if CLIPPING_left_alignmentObj != None:
                        clippingEventsDict['LEFT-CLIPPING'].append(CLIPPING_left_alignmentObj)
                    if CLIPPING_right_alignmentObj != None:
                        clippingEventsDict['RIGHT-CLIPPING'].append(CLIPPING_right_alignmentObj)
                except TypeError:
                    continue

        ## When the discordant cluster is RIGHT, add the biggest right clipping cluster if any:
        if all (event.side == 'PLUS' for event in self.events):
            ## Get clipping clusters:
            clippingRightEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == 'RIGHT-CLIPPING')
            CLIPPING_cluster = self.add_clippingEvents(ref, binBeg, binEnd, clippingRightEventsDict, ['RIGHT-CLIPPING'], confDict)

        ## When the discordant cluster is LEFT, add the biggest left clipping cluster if any:
        elif all (event.side == 'MINUS' for event in self.events):
            ## Get clipping clusters:
            clippingLeftEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == 'LEFT-CLIPPING')
            CLIPPING_cluster = self.add_clippingEvents(ref, binBeg, binEnd, clippingLeftEventsDict, ['LEFT-CLIPPING'], confDict)

        # TODO si es reciproco:
        else:
            CLIPPING_cluster = self.add_clippingEvents(ref, binBeg, binEnd, clippingEventsDict, ['RIGHT-CLIPPING', 'LEFT-CLIPPING'], confDict)

        return CLIPPING_cluster

    def add_clippingEvents(self, ref, binBeg, binEnd, clippingEventsDict, eventTypes, confDict):
        binSizes = [100, 1000]
        clippingBinDb = structures.create_bin_database(ref, binBeg, binEnd, clippingEventsDict, binSizes)
        binSize = clippingBinDb.binSizes[0]
        CLIPPING_clusters = clustering.distance_clustering(clippingBinDb, binSize, eventTypes, 'CLIPPING', confDict['maxEventDist'], confDict['minClusterSize']) 

        # Si hay algun clipping cluster:
        if len (CLIPPING_clusters) > 0:
            ## Choose the clipping cluster with the highest number of events:
            # Coger el cluster de la lista de clusters si si length es igual a la maxima length de todos los clusters de la lista. (como devuelve una lista de un solo elemento, cojo el primer elemento de la lista.)
            # TODO: si hay dos con el mismo numero de eventos.
            CLIPPING_cluster = [cluster for cluster in CLIPPING_clusters if len(cluster.events) == max([len(cluster.events) for cluster in CLIPPING_clusters])][0]

            ## Add cluster's reads to the discordant metacluster:
            self.add(CLIPPING_cluster.events)

            return CLIPPING_cluster

            ## Remove events from discordant cluster that are higher (more to the right) than the clippingEnd
            # discordantCluster.removeDiscordant(clippingEnd, 'right')

    def setIntOrigin(self):
        '''
        Set germline or somatic
        '''
        nbNormalEvents = 0

        for event in self.events:
            if event.sample == "NORMAL":
                nbNormalEvents += 1
        # HACER ESTO ARGUMENTOO!!!
        if nbNormalEvents > 3:
            self.intOrigin = 'germline'
        else:
            self.intOrigin = 'somatic'