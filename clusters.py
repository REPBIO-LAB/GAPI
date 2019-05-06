'''
Module 'variants' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External
import sys
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
    #elif (clusterType == 'DISCORDANT'):
    #    cluster = DISCORDANT_cluster(events)

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

def create_metaclusters(eventsBinDb, confDict):
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
        eventTypes = ['DISCORDANT']
        # TO DO  
        # metaclusters = clustering.reciprocal_clustering()
        # allMetaclusters = allMetaclusters + metaclusters

    for metacluster in allMetaclusters:
        print('METACLUSTER: ', metacluster, len(metacluster.events), [(clusterType, len(subcluster.events)) for clusterType, subcluster in metacluster.subclusters.items()])


    ## 3. Organize metaclusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    metaclustersDict = {}
    metaclustersDict['METACLUSTERS'] = allMetaclusters 

    metaclustersBinDb = structures.create_bin_database(eventsBinDb.ref, eventsBinDb.beg, eventsBinDb.end, metaclustersDict, binSizes)
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
        1. : 
    ''' 
    for cluster in clustersBinDb.collect([clusterType]):

        print('CREATE_CONSENSUS_FOR: ', cluster, cluster.ref, cluster.beg, cluster.end)

        clusterId = '_'.join([str(cluster.ref), str(cluster.beg), str(cluster.end)])
        outDir = rootOutDir + '/' + clusterId
        SV_type, consensus = cluster.make_consensus(confDict, reference, outDir)

        print('CONSENSUS_RESULT: ', SV_type, consensus)
        print('-------------------------------')

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
    
        # Set some metacluster properties as None
        self.filters = None
        self.SV_type = None
        self.consensus = None

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

            print('A) Metacluster contains an insertion cluster')

            ## Select INS event with median length as template
            templateEvent = self.subclusters['INS'].pick_median_length()

            ## Write template into output file 
            templateFasta = formats.FASTA()
            templateFasta.seqDict[templateEvent.readName] = templateEvent.readSeq
            templateFile = outDir + '/template.fa'
            templateFasta.write(templateFile)   

        # B) Metacluster contains a DEL cluster
        elif ('DEL' in self.subclusters):

            print('B) Metacluster contains a deletion cluster')

            ## Select DEL event with median length as template
            templateEvent = self.subclusters['DEL'].pick_median_length()

            ## Write template into output file 
            templateFasta = formats.FASTA()
            templateFasta.seqDict[templateEvent.readName] = templateEvent.readSeq
            templateFile = outDir + '/template.fa'
            templateFasta.write(templateFile)   

        # C) Metacluster composed by only two CLIPPING clusters (left and right) 
        elif all (clusterType in self.subclusters for clusterType in ['LEFT-CLIPPING', 'RIGHT-CLIPPING']):
            
            print('C) Metacluster composed by only two clipping clusters (left and right)')
        
            ## Search for chimeric alignment spanning the SV event
            templateEvent, supplementary, chimeric = find_chimeric_alignments(self.subclusters['RIGHT-CLIPPING'], self.subclusters['LEFT-CLIPPING'])

            # a) Chimeric alignment found -> Write template into output file 
            if (templateEvent != None): 

                print('C.a) Chimeric alignment found')

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

                if (templateFile != None): 
                    print('C.b) complementary clippings found')

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

        if templateFile == None:
            print('TEMPLATE EXTRACTION FAILED. SKIP CONSENSUS GENERATION')
            return None, None
            
        template = formats.FASTA()
        template.read(templateFile)
        templateReadName = list(template.seqDict.keys())[0]

        ## 2. Collect metacluster supporting reads 
        supportingReads = self.collect_reads()

        ## Remove template from FASTA
        if templateReadName in supportingReads.seqDict:
            del supportingReads.seqDict[templateReadName]

        ## Write supporting reads FASTA 
        supportingReadsFile = outDir + '/supportingReads.fa'
        supportingReads.write(supportingReadsFile)
            
        ## 3. Template polishing to generate consensus
        consensusFile = assembly.polish_racon(templateFile, supportingReadsFile, confDict['technology'], 2, outDir)

        ## If polishing failed use directly the template as consensus
        if consensusFile == None:
            print('POLISHING FAILED. USE TEMPLATE AS CONSENSUS')
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
        if BAM == None:
            print('CONSENSUS REALIGNMENT FAILED')

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

        ## 5. Define metacluster type and properties based on events collected from consensus sequence realignment
        # A) Single INS event
        if len(eventsDict['INS']) == 1:
            print('CONSENSUS!! A) Single INS event')

            ## Set metacluster type
            SV_type = 'INS'
                
            ## Convert coordinates
            event = eventsDict['INS'][0]
            event = alignment.targetered2genomic_coord(event, self.ref, intervalBeg)

            ## Incorporate INS in the metacluster as consensus  
            consensus = event
            
        # B) Multiple INS events 
        # Raw alignment    -------------[INS]-[INS]---[INS]-------------
        # Consensus        -------------[       INS       ]-------------
        # This is consequence of fragmented alignments. So Merge events into a single consensus INS
        elif len(eventsDict['INS']) > 1:
            print('CONSENSUS!! B) Multiple INS events')

            ## Set metacluster type
            SV_type = None

            ## Do merging

            ## Convert coordinates

            ## Incorporate merged INS in the metacluster as consensus  
            consensus = None

        # C) Single DEL event
        elif len(eventsDict['DEL']) == 1:
            print('CONSENSUS!! C) Single DEL event')

            ## Set metacluster type
            SV_type = 'DEL'
                
            ## Convert coordinates
            event = eventsDict['DEL'][0]
            event = alignment.targetered2genomic_coord(event, self.ref, intervalBeg)

            ## Incorporate DEL in the metacluster as consensus  
            consensus = event

        # D) Multiple DEL events (same scenario as B) 
        elif len(eventsDict['DEL']) > 1:
            print('CONSENSUS!! D) Multiple DEL events (same scenario as B)')

            ## Set metacluster type
            SV_type = None

            ## Do merging

            ## Convert coordinates

            ## Incorporate merged DEL in the metacluster as consensus 
            consensus = None

        # E) One left and one right CLIPPING (NEXT TO DO)
        elif (len(eventsDict['LEFT-CLIPPING']) == 1) and (len(eventsDict['RIGHT-CLIPPING']) == 1):
            print('CONSENSUS!! E) One left and one right CLIPPING')

            ## Set metacluster type
            SV_type = 'INS'

            ## Create INS event
            rightClipping = eventsDict['RIGHT-CLIPPING'][0]
            leftClipping = eventsDict['LEFT-CLIPPING'][0]
            length = leftClipping.readBkp - rightClipping.readBkp
            event = events.INS(rightClipping.ref, rightClipping.beg, rightClipping.end, length, rightClipping.readName, rightClipping.readSeq, rightClipping.readBkp, None, None)
        
            ## Convert coordinates
            event = alignment.targetered2genomic_coord(event, self.ref, intervalBeg)

            ## Incorporate INS in the metacluster as consensus  
            consensus = event

        # F) Single left CLIPPING
        elif (len(eventsDict['LEFT-CLIPPING']) == 1): 
            print('CONSENSUS!! F) Single left CLIPPING')
            SV_type = None
            consensus = None

        # G) Single right CLIPPING
        elif (len(eventsDict['RIGHT-CLIPPING']) == 1): 
            print('CONSENSUS!! G) Single right CLIPPING')
            SV_type = None
            consensus = None

        # H) Another possibility
        else:
            print('CONSENSUS!! H) Another possibility')
            SV_type = None
            consensus = None
    
        return SV_type, consensus