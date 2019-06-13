'''
Module 'clusters' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External
import sys
import numpy as np
import collections 
import pysam
import math

# Internal
import log
import formats
import unix
import clustering
import events
import structures
import alignment
import bamtools 
import assembly
import repeats
import retrotransposons
import filters
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


def create_discordantClusters(discordantBinDb, minClusterSize, buffer):
    '''
    Group discordant read pairs according to cluster type and reciprocal overlap into clusters 

    Input:
        1. discordantBinDb: data structure containing a set of discordant read pairs organized in genomic bins  
        2. minClusterSize: minimum cluster size
        3. buffer: number of base pairs to extend begin and end coordinates for each discordant prior clustering


    Output:
        1. discordantClustersDict: dictionary containing for each possible discordant cluster type (keys) a list of clusters (values)
    '''
    discordantClustersDict = {}

    # For each discordant read pair cluster type
    for clusterType in discordantBinDb.eventTypes:
        
        # Do clustering based on reciprocal overlap
        discordantClustersDict[clusterType] = clustering.reciprocal_overlap_clustering(discordantBinDb, 1, minClusterSize, clusterType, buffer, clusterType)

    return discordantClustersDict


def create_metaclusters(eventsBinDb, confDict):
    '''
    Group different types of SV events into metaclusters 

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
        # metaclusters = clustering.reciprocal_overlap_clustering()   
        # allMetaclusters = allMetaclusters + metaclusters

    # 2.5) Create DISCORDANT metaclusters
    if 'DISCORDANT' in confDict['targetSV']:
        eventTypes = eventsBinDb.eventTypes
        
        for eventType in eventTypes:
            
            metaclusters = clustering.reciprocal_overlap_clustering(eventsBinDb, 1, confDict['minClusterSize'], eventType, 0, 'META')
            allMetaclusters = allMetaclusters + metaclusters

    if 'DISCORDANT' in confDict['targetSV']:
        eventTypes = eventsBinDb.eventTypes
        for eventType in eventTypes:
            discordantClustersDict[eventType] = clustering.reciprocal_overlap_clustering(eventsBinDb, 1, confDict['minClusterSize'], eventType, 0, eventType)

    ## 3. Organize metaclusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    metaclustersDict = {}
    metaclustersDict['METACLUSTERS'] = allMetaclusters 

    metaclustersBinDb = structures.create_bin_database_interval(eventsBinDb.ref, eventsBinDb.beg, eventsBinDb.end, metaclustersDict, binSizes)
    return metaclustersBinDb


def make_consensus(clustersBinDb, confDict, reference, clusterType, rootOutDir):
    '''
    Make consensus supporting sequence and event for a set of cluster objects

    Input:
        1. clustersBinDb: Data structure containing a set of clusters organized in genomic bins
        2. confDict: 
            * technology     -> sequencing technology (NANOPORE, PACBIO or ILLUMINA)
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)

        3. reference: path to reference genome in fasta format    
        4. clusterType: target cluster type to generate consensus
        5. rootOutDir: root output directory

    Output:
        1. consensusBinDb: bin database structure containing metaclusters for which consensus sequences and events have been generated
    ''' 
    consensusDict = {}

    ## For each cluster in the database
    for cluster in clustersBinDb.collect([clusterType]):

        outDir = rootOutDir + '/Consensus/' + str(cluster.id)

        cluster.SV_type, cluster.consensus = cluster.make_consensus(confDict, reference, outDir)

        # Discard clusters without known SV type
        if cluster.SV_type is not None:

            # Initialize SV type at the dictionary
            if cluster.SV_type not in consensusDict:
                consensusDict[cluster.SV_type] = []

            # Add cluster 
            consensusDict[cluster.SV_type].append(cluster)
            
        else:
            print('UNKNOWN_TYPE: ', cluster)

    ## 3. Organize metaclusters into bins according to their SV type    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    consensusBinDb = structures.create_bin_database_interval(clustersBinDb.ref, clustersBinDb.beg, clustersBinDb.end, consensusDict, binSizes)

    return consensusBinDb


def find_chimeric_alignments(clusterA, clusterB):
    '''
    Search for chimeric read alignments connecting two clipping clusters. Select one as representative if multiple are identified. 

    Input:
        1. clusterA: Clipping cluster object
        2. clusterB: Clipping cluster object

    Output:
        1. primary: clipping event for representative primary alignment
        2. supplementary: clipping event for representative supplementary alignment
        3. chimericSorted: list of tuples. Each tuple is composed by two clipping events corresponding to primary and supplementary alignments, respectively. 
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
        return None, None, None

    ### 2. Select the clipping events with the longest supplementary alignment as representative
    chimericSorted = sorted(chimeric, key=lambda alignments: alignments[1].refLen, reverse=True)
    primary = chimericSorted[0][0]
    supplementary = chimericSorted[0][1]

    return primary, supplementary, chimericSorted


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


def determine_INS_type(metaclusters, index, confDict, rootOutDir):
    '''
    Function to determine what has been inserted for each cluster. 

    Input:
        1. metaclusters: list containing metaclusters supporting INS events
        2. index: Minimap2 index for fasta file containing retrotransposon related sequences 
        3. confDict: ...
        4: rootOutDir: root directory to write files and directories
    '''

    ## For each metacluster in the list
    for metacluster in metaclusters:

        ## Determine the insertion type
        outDir = rootOutDir + '/INS_type/' + str(metacluster.id)
        metacluster.determine_INS_type(index, confDict, outDir)


def merge_fragmented_INDELS(metaclusters):
    '''
    Identify and merge fragmented alignments over INDELs

    Input:
        1. metaclusters: list of metaclusters  

    Output:
        1. Modify metacluster instances 
    '''

    ## For each metacluster in the list
    for metacluster in metaclusters:
        metacluster.merge_fragmented_INDELS()


#############
## CLASSES ##
#############

class cluster():
    '''
    Events cluster class. A cluster is composed by a set of events.
    Each event is supported by a single read. One cluster can completely represent a single structural
    variation event or partially if multiple clusters are required (see 'metaCluster' class)
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

    def sort_by_length(self):
        '''
        Sort events in increasing length ordering
        '''
        return sorted(self.events, key=lambda event: event.length)
        
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

    def pick_median_length(self):
        '''
        Return event whose length is at the median amongst all cluster supporting events
        '''
        ## Sort events by their length
        sortedEvents = self.sort_by_length()

        ## Compute the index for the event with the median length
        median = (len(sortedEvents) - 1)/2  # minus 1 because the first element is index 0
        medianIndex = int(math.ceil(median))

        ## Pick event located at the median 
        return sortedEvents[medianIndex]

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

    def correct_fragmentation(self):
        '''
        Identify and merge fragmented alignments over INS
        
        Before merging:
        ############<<<INS>>>##<<INS>>###<<<INS>>>##############
        
        After merging:
        ############<<<<<<<<INS>>>>>>>>##############
        
        Output:
            1. Modify INS cluster instance 
            2. merged_list: list of merged INS objects
            3. fragmented_list: list of fragmented INS events that has been merged
        '''
        ## 1. Organize INS events into a dictionary according to their supporting read
        eventsByReads =  {}

        for INS in self.events:
        
            # a) First event supported by that read
            if INS.readName not in eventsByReads:
                eventsByReads[INS.readName] = [INS]

            # b) There are previous events supported by that read
            else:
                eventsByReads[INS.readName].append(INS)
        
        ## 2. Merge INS events supported by the same read
        merged_list = []
        fragmented_list = []

        # For each read
        for readId, INS_list in eventsByReads.items():
            
            ## Read supporting multiple INS events
            if len(INS_list) > 1:

                ## 2.1 Do merging of fragmented INS
                merged = events.merge_INS(INS_list)

                ## 2.2 Delete INS that have been merged 
                self.events = [INS for INS in self.events if INS not in INS_list]

                ## 2.3 Add merged INS
                merged_list.append(merged)
                self.events.append(merged)             

                ## 2.4 Update fragmented alignments list
                fragmented_list = fragmented_list + INS_list
    
        return merged_list, fragmented_list

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

        # Organize events into subclusters
        self.subclusters = self.create_subclusters()

        # Tag germline or somatic
        self.intOrigin = None

        # Set some metacluster properties as None
        self.filters = None
        self.consensus = None
        self.SV_type = None

        # Metacluster SV features
        self.INS_features = {}
        self.DEL_features = {}
        self.BKP_features = {}

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
                
    def collect_reads(self):
        '''
        Create FASTA object containing metacluster supporting reads.
        
        Output:
            1. FASTA: FASTA object containing metacluster supporting reads
        '''
        ## Initiate FASTA object
        FASTA = formats.FASTA()

        ## For each event composing the metacluster         
        for event in self.events: 

            ## A) Read sequence supporting the event not included in the FASTA yet
            if event.readName not in FASTA.seqDict:
                FASTA.seqDict[event.readName] = event.readSeq 

            ## B) Read sequence supporting the event already included in the FASTA
            else:

                # Current read sequence aligned as primary -> replace previously included sequence (therefore supplementary)
                if not event.supplementary:
                    FASTA.seqDict[event.readName] = event.readSeq 
     
        return FASTA

    def merge_fragmented_INDELS(self):
        '''
        Identify and merge fragmented alignments over INDELs
        
        **** DEL fragmentation ****
        Before merging:
        ############---DEL---##--DEL--###---DEL---##############
        
        After merging:
        ############----------DEL----------##############
        '''
        ## Metacluster contains an INS cluster
        if 'INS' in self.subclusters:

            ## 1) Correct INS fragmentation at cluster level
            merged, fragmented = self.subclusters['INS'].correct_fragmentation()

            ## 2) Remove fragmented INS events from metacluster
            self.events = [event for event in self.events if event not in fragmented]

            ## 3) Add merged INS events to the metacluster
            self.events = self.events + merged

        ## MEtacluster contains a DEL cluster
        #if 'DEL' in self.subclusters:
            #self.subclusters['DEL'].correct_fragmentation()

        # Sort events from lower to upper beg coordinates
        self.sort()  

    def nbEvents(self):
        '''
        Return the number of events composing the metacluster. 
        '''
        ## Initialize counters
        nbTumour = 0
        nbNormal = 0

        nbINS = 0
        nbDEL = 0
        nbCLIPPING = 0   

        # For each event composing the metacluster
        for event in self.events:

            ## Tumour and matched normal counts
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
            
            ## Event type counts
            # a) INS event
            if event.type == 'INS':
                nbINS += 1

            # b) DEL event
            elif event.type == 'DEL':
                nbDEL += 1

            # c) CLIPPING event
            else:
                nbCLIPPING += 1            

        nbTotal = len(self.events)

        return nbTotal, nbTumour, nbNormal, nbINS, nbDEL, nbCLIPPING
        
    
    def supportingCLIPPING(self, buffer, confDict, bam, normalBam, mode):
        # Note: This function works but you have to allow duplicates in the clipping 

        # Make custom conf. dict for only selecting duplicates
        clippingConfDict = dict(confDict)
        clippingConfDict['targetSV'] = ['CLIPPING']
        clippingConfDict['minMAPQ'] = 0

        clippingEventsDict = {}

        ## Define region
        binBeg = self.beg - buffer if self.beg > buffer else 0
        
        # TODO check as above
        binEnd = self.end

        ref = self.ref

        if mode == "SINGLE":
            clippingEventsDict = bamtools.collectSV(ref, binBeg, binEnd, bam, clippingConfDict, None)
        elif mode == "PAIRED":
            clippingEventsDict = bamtools.collectSV_paired(ref, binBeg, binEnd, bam, normalBam, clippingConfDict)

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

        # TODO if it is reciprocal
        else:
            CLIPPING_cluster = self.add_clippingEvents(ref, binBeg, binEnd, clippingEventsDict, ['RIGHT-CLIPPING', 'LEFT-CLIPPING'], confDict)

        return CLIPPING_cluster
        
    def add_clippingEvents(self, ref, binBeg, binEnd, clippingEventsDict, eventTypes, confDict):
        binSizes = [100, 1000]
        clippingBinDb = structures.create_bin_database_interval(ref, binBeg, binEnd, clippingEventsDict, binSizes)
        binSize = clippingBinDb.binSizes[0]
        CLIPPING_clusters = clustering.distance_clustering(clippingBinDb, binSize, eventTypes, 'CLIPPING', confDict['maxEventDist'], confDict['minClusterSize']) 

        # If there is a clipping cluster
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

        return nbTotal, nbTumour, nbNormal, nbINS, nbDEL, nbCLIPPING

           
    def select_template(self, technology, outDir):
        '''
        Select one metacluster supporting read as template to be used in the metacluster consensus generation step

        Input: 
            1. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
            2. outDir: output directory

        Output:
            1. templateFile: Path to FASTA file containing the template or 'None' if no template was found
        '''

        # A) Metacluster contains an INS cluster
        if ('INS' in self.subclusters):

            ## Select INS event with median length as template
            templateEvent = self.subclusters['INS'].pick_median_length()

            ## Write template into output file 
            templateFasta = formats.FASTA()
            templateFasta.seqDict[templateEvent.readName] = templateEvent.readSeq
            templateFile = outDir + '/template.fa'
            templateFasta.write(templateFile)   

        # B) Metacluster contains a DEL cluster
        elif ('DEL' in self.subclusters):

            ## Select DEL event with median length as template
            templateEvent = self.subclusters['DEL'].pick_median_length()

            ## Write template into output file 
            templateFasta = formats.FASTA()
            templateFasta.seqDict[templateEvent.readName] = templateEvent.readSeq
            templateFile = outDir + '/template.fa'
            templateFasta.write(templateFile)   

        # C) Metacluster composed by only two CLIPPING clusters (left and right) 
        elif all (clusterType in self.subclusters for clusterType in ['LEFT-CLIPPING', 'RIGHT-CLIPPING']):
                    
            ## Search for chimeric alignment spanning the SV event
            templateEvent, supplementary, chimeric = find_chimeric_alignments(self.subclusters['RIGHT-CLIPPING'], self.subclusters['LEFT-CLIPPING'])

            # a) Chimeric alignment found -> Write template into output file 
            if (templateEvent != None): 

                templateFasta = formats.FASTA()
                templateFasta.seqDict[templateEvent.readName] = templateEvent.readSeq
                templateFile = outDir + '/template.fa'
                templateFasta.write(templateFile)   

            # b) Chimeric alignment NOT found -> search for complementary clippings
            else:

                ## Generate fasta files containing clusters supporting reads:
                readsA = self.subclusters['RIGHT-CLIPPING'].collect_reads() 
                readsB = self.subclusters['LEFT-CLIPPING'].collect_reads()

                readsA_file = outDir + '/seqA.fa'
                readsB_file = outDir + '/seqB.fa'

                readsA.write(readsA_file)
                readsB.write(readsB_file)

                ## Assemble clippings based on overlap 
                templateFile = assembly.assemble_overlap(readsA_file, readsB_file, technology, outDir)

        # D) Metacluster composed by a single CLIPPING cluster. 
        # For now set the template as None, but at one point I should set a criteria for picking a template
        # in these cases. I think the template will be the read with the longest piece of clipped sequence
        else:
            templateFile = None
        
        return templateFile


    def make_consensus(self, confDict, reference, outDir):
        '''
        Make consensus supporting sequence and event for the metacluster

        Input: 
            1. confDict: 
                * technology     -> sequencing technology (NANOPORE, PACBIO or ILLUMINA)
                * rounds         -> number of polishing rounds to be attempled. 0 means no polishing
                * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
                * minMAPQ        -> minimum mapping quality
                * minCLIPPINGlen -> minimum clipping lenght
                * minINDELlen    -> minimum INS and DEL lenght
                * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)            2. outDir: output directory
            
            2. reference: path to reference genome in fasta format    
            3. outDir: Output directory
        
        Output:
            1. SV_type: Structural variation type (INS, DEL, BKP or None)
            2. consensus: Consensus SV event object
        '''

        ## 0. Create directory 
        unix.mkdir(outDir)

        ## 1. Define template
        templateFile = self.select_template(confDict['technology'], outDir)

        if templateFile is None:
            return None, None
        
        ## 2. Collect metacluster supporting reads 
        supportingReads = self.collect_reads()

        ## Remove template from FASTA
        template = formats.FASTA()
        template.read(templateFile)
        templateReadName = list(template.seqDict.keys())[0]

        if templateReadName in supportingReads.seqDict:
            del supportingReads.seqDict[templateReadName]

        ## Write supporting reads FASTA 
        supportingReadsFile = outDir + '/supportingReads.fa'
        supportingReads.write(supportingReadsFile)
            
        ## 3. Template polishing to generate consensus
        consensusFile = assembly.polish_racon(templateFile, supportingReadsFile, confDict['technology'], confDict['rounds'], outDir)

        ## If no polished sequence obtained use directly the template as consensus
        if consensusFile == None:
            consensusFile = templateFile

        ## 4. Realign consensus sequence into the SV event genomic region
        ##  Define SV cluster surrounding region
        # ------------------<***SV_cluster***>-----------------
        # <--nbFlankingBp-->                  <--nbFlankingBp-->
        nbFlankingBp = 10000
        intervalBeg = self.beg - nbFlankingBp
        intervalBeg = intervalBeg if intervalBeg >= 0 else 0 ## Set lower bound
        intervalEnd = self.end + nbFlankingBp
        intervalCoord = self.ref + ':' + str(intervalBeg) + '-' + str(intervalEnd)
            
        ## Do realignment
        # ------------------<***SV_cluster***>-----------------
        #         -------------consensus_seq-------------
        BAM = alignment.targeted_alignment_minimap2(consensusFile, intervalCoord, reference, outDir)
            
        ## Consensus realignment failed
        if BAM is None:
            return None, None

        ## Extract events from consensus sequence realignment
        # ------------------<***SV_cluster***>-----------------
        #        -------------consensus_seq-------------
        #               <--->----------------<---> overhang (100 bp)
        #                   event_search_space         
        overhang = 100
        clusterIntervalLen = self.end - self.beg
        targetBeg = nbFlankingBp - overhang
        targetEnd = nbFlankingBp + clusterIntervalLen + overhang            
        eventsDict = bamtools.collectSV(intervalCoord, targetBeg, targetEnd, BAM, confDict, None)

        ## 5. Define metacluster type based on events collected from consensus sequence realignment
        # A) Single INS event
        if len(eventsDict['INS']) == 1:

            ## Set metacluster type
            SV_type = 'INS'
                
            ## Convert coordinates
            consensus = alignment.targetered2genomic_coord(eventsDict['INS'][0], self.ref, intervalBeg)

        # B) Multiple INS events 
        # Raw alignment    -------------[INS]-[INS]---[INS]-------------
        # Consensus        -------------[       INS       ]-------------
        # This is consequence of fragmented alignments. So Merge events into a single consensus INS
        elif len(eventsDict['INS']) > 1:

            ## Set metacluster type
            SV_type = 'INS'

            ## Do merging
            merged = events.merge_INS(eventsDict['INS'])

            ## Convert coordinates
            consensus = alignment.targetered2genomic_coord(merged, self.ref, intervalBeg)

        # C) Single DEL event
        elif (len(eventsDict['DEL']) == 1) and ('DEL' in self.subclusters):

            ## Set metacluster type
            SV_type = 'DEL'
                
            ## Convert coordinates
            consensus = alignment.targetered2genomic_coord(eventsDict['DEL'][0], self.ref, intervalBeg)

        # D) Multiple DEL events 
        elif (len(eventsDict['DEL']) > 1) and ('DEL' in self.subclusters):

            ## Set metacluster type
            SV_type = None

            ## Do merging

            ## Convert coordinates
            consensus = None

        # E) One left and one right CLIPPING 
        elif (len(eventsDict['LEFT-CLIPPING']) == 1) and (len(eventsDict['RIGHT-CLIPPING']) == 1):

            ## Set metacluster type
            SV_type = 'INS'

            ## Create INS event
            rightClipping = eventsDict['RIGHT-CLIPPING'][0]
            leftClipping = eventsDict['LEFT-CLIPPING'][0]
            length = leftClipping.readBkp - rightClipping.readBkp

            # a) Misalignment leading to aberrant clipping 
            # (rarely happens, at one point investigate further)
            if length <= 0:
                SV_type = None
                consensus = None

            # b) Correct clipping
            else:
                event = events.INS(rightClipping.ref, rightClipping.beg, rightClipping.end, length, rightClipping.readName, rightClipping.readSeq, rightClipping.readBkp, None, None)
        
                ## Convert coordinates
                consensus = alignment.targetered2genomic_coord(event, self.ref, intervalBeg)

        # F) Another possibility 
        else:

            # a) Metacluster contains an INS cluster
            if ('INS' in self.subclusters):

                ## Set metacluster type
                SV_type = 'INS'

                ## Select INS event with median length as consensus
                consensus = self.subclusters['INS'].pick_median_length()

            # b) No INS cluster composing the metacluster
            else:
                SV_type = None
                consensus = None
    
        ## Cleanup
        unix.rm([outDir])

        return SV_type, consensus

    def determine_INS_type(self, index, confDict, outDir): 
        '''
        Determine the type of insertion (retrotransposon, simple repeat, virus, ...) and collect insertion information.

        Input: 
            1. index: Minimap2 index for fasta file containing retrotransposon related sequences 
            2. confDict: Configuration dictionary 
            3. outDir: Output directory
        '''

        ## 0. Create output directory 
        unix.mkdir(outDir)

        ## 1. Pick inserted sequence
        insert = self.consensus.pick_insert()

        ## 2. Write target seq into file 
        FASTA = formats.FASTA()
        FASTA.seqDict[str(self.id)] = insert

        FASTA_file = outDir + '/insert.fa'
        FASTA.write(FASTA_file)

        ## 3. Determine to what corresponds the insertion
        # A) Tandem repeat/simple repeat expansion? 
        minPercSimple = 70
        self.INS_features['insType'], self.INS_features['status'], self.INS_features['percResolved'] = repeats.is_simple_repeat(FASTA_file, minPercSimple, outDir)

        ## Stop if insertion classified as simple repeat
        if self.INS_features['status'] == 'resolved':
            return

        # B) Retrotransposon insertion? 
        self.INS_features['insType'], self.INS_features['family'], self.INS_features['srcId'], self.INS_features['status'], self.INS_features['percResolved'], self.INS_features['strand'], self.INS_features['hits'] = retrotransposons.is_retrotransposition(FASTA_file, index, outDir)
        
        ## Stop if insertion classified as retrotransposon
        if (self.INS_features['status'] == 'resolved') or (self.INS_features['status'] == 'partially_resolved'):
            return
       
        # C) Viral insertion? 
        #viruses.is_virus()

        # D) Telomemric insertion? 

        # E) Mitochondrial or rearranged DNA insert?
        #¿¿¿rearrangements???.is_chromosomal_dna()
        
        ## Cleanup
        unix.rm([outDir])
