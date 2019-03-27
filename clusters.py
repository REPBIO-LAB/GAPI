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
import repeats
import retrotransposons

###############
## FUNCTIONS ##
###############

def createCluster(events, clusterType):
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
        cluster = META_cluster(events, clusterType)

    ## b) Create INS cluster
    elif (clusterType == 'INS'):
        cluster = INS_cluster(events, clusterType)

    ## c) Create DEL cluster
    elif (clusterType == 'DEL'):
        cluster = DEL_cluster(events, clusterType)

    ## d) Create CLIPPING cluster
    elif (clusterType == 'CLIPPING') or (clusterType == 'LEFT-CLIPPING') or (clusterType == 'RIGHT-CLIPPING'):
        cluster = CLIPPING_cluster(events, clusterType)

    ## e) Create DISCORDANT cluster
    #elif (clusterType == 'DISCORDANT'):
    #    cluster = DISCORDANT_cluster(events, clusterType)

    ## f) Unexpected cluster type
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
        
def consensusClusters(clusters, clusterType, supportingReads, reference, confDict, rootDir):
    '''
    Function to create a consensus sequence for each cluster. 

    Input:
        1. clusters: bin database containing a set of cluster objects
        2. clusterType: type of cluster (INS-CLUSTER: insertion; DEL-CLUSTER: deletion; LEFT-CLIPPING-CLUSTER: left clipping; RIGHT-CLIPPING-CLUSTER: right clipping)
        3. supportingReads: FASTQ (if qualities available) or FASTA (if qualities NOT available) object containing all the reads supporting SV clusters 
        4. reference: path to reference genome in fasta format
        5. confDict: configuration dictionary (Complete...)
        6: rootDir: root directory to write files and directories
    '''

    outDir = rootDir + '/Consensus/'
    unix.mkdir(outDir)

    ## 1. Create consensus sequence for each cluster
    for cluster in clusters.collect(clusterType):        

        ## Attempt to create a consensus if: 
        # a) No filter has been applied to the clusters OR
        # b) MAX-NBREADS filter has not been applied to the clusters OR
        # c) MAX-NBREADS filter has been applied and the cluster passes the filter 
        # Done to skip consensus generation for clusters with an abnormally high number of supporting reads. They are artefacts and takes a very long time to process them
        if (cluster.filters == None) or ('MAX-NBREADS' not in cluster.filters) or (cluster.filters['MAX-NBREADS']):
            
            ## Create consensus
            cluster.create_consensus(supportingReads, reference, confDict, outDir)

def insTypeClusters(clusters, index, confDict, rootDir):
    '''
    Function to determine what has been inserted for each cluster. 

    Input:
        1. clusters: bin database containing a set of cluster objects
        2. index: Minimap2 index for fasta file containing retrotransposon related sequences 
        3. confDict: type of cluster (INS-CLUSTER: insertion; DEL-CLUSTER: deletion; LEFT-CLIPPING-CLUSTER: left clipping; RIGHT-CLIPPING-CLUSTER: right clipping)
        4: rootDir: root directory to write files and directories
    '''
    outDir = rootDir + '/InsType/'
    unix.mkdir(outDir)

    ## 1. Determine insertion type for each cluster
    for cluster in clusters.collect('INS-CLUSTER'):

        ## Determine insertion type
        cluster.determine_insertion_type(index, confDict, outDir)


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

        # TEMPORARY (REMOVE ONCE OUTPUT VCF IS IMPLEMENTED!!!!!)
        # Inserted seq
        self.consensus_FASTA = None
        self.consensusLen = None
        self.isConsensus = None
        self.insertSeq = None

        # Cluster features
        self.status = None
        self.insType = None
        self.family = None 
        self.srcId = None
        self.percResolved = None
        self.strand = None
        self.hits = None

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
            
    def supportingReadIds(self):
        '''
        Return list of ids for cluster supporting reads
        '''
        readIds = [event.readId for event in self.events]

        return readIds

    def collect_supportingReads(self, reads, quality):
        '''
        Collect cluster supporting reads from a FASTQ/FASTA object 
        
        Input:
            1. reads: FASTQ or FASTA object containing a set of reads
            2. quality: True (sequence qualities available) or False (not available). 
        Output:
            1. supportingReads: FASTQ (if qualities available) or FASTA (if qualities NOT available) object containing the cluster supporting reads
        '''
        readIds = self.supportingReadIds()

        # a) Quality available  
        if quality:
            clusterSupportingReads = formats.FASTQ()

            # Collect from the FASTQ containing all the reads only those supporting the SV cluster 
            for readId in readIds:
                entry = reads.seqDict[readId]
                clusterSupportingReads.add(entry)

        # b) Quality not available
        else:
            clusterSupportingReads = formats.FASTA()

            # Collect from the FASTA containing all the reads only those supporting the SV cluster 
            for readId in readIds:
                seq = reads.seqDict[readId]
                clusterSupportingReads.seqDict[readId] = seq

        return clusterSupportingReads

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


class META_cluster(cluster):
    '''
    Meta cluster subclass
    '''
    def __init__(self, events, clusterType):

        cluster.__init__(self, events, clusterType)

class INS_cluster(cluster):
    '''
    Insertion (INS) cluster subclass
    '''
    def __init__(self, events, clusterType):

        cluster.__init__(self, events, clusterType)

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

    def polish(self):
        '''
        Apply successive rounds of polishing to the INS cluster by removing events whose length deviates from the cluster average
        '''
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

    def create_consensus(self, supportingReads, reference, confDict, rootDir): 
        '''
        Use all the cluster supporting reads to (1) generate a consensus sequence for the cluster and (2) perform local realignment of the consensus 
        to obtain more accurate SV breakpoints and inserted sequence (only for INS)

        Update cluster with all this information
        
        Input:
            1. supportingReads: FASTQ (if qualities available) or FASTA (if qualities NOT available) object containing all the reads supporting SV clusters 
            2. reference: Path to the reference genome in fasta format. An index of the reference generated with samtools faidx must be located in the same directory
            3. confDict: Configuration dictionary 
            4. rootDir: Root directory where temporary folders and files will be created
        '''
        ## 0. Create output directory
        outDir = rootDir + '/' + str(self.id)
        unix.mkdir(outDir)

        ## 1. Collect cluster supporting reads from the FASTA/FASTQ containing all the reads supporting SV clusters
        clusterSupportingReads = self.collect_supportingReads(supportingReads, confDict['quality'])

        ## 2. Generate consensus sequence for the cluster
        self.consensus_FASTA = consensus.racon(clusterSupportingReads, confDict['technology'], confDict['quality'], outDir)

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
        INS_events, DEL_events, CLIPPING_left_events, CLIPPING_right_events, supportingReads = bamtools.collectSV(targetInterval, 0, targetLen, BAM, confDict, None)
               
        ## 5. Redefine cluster properties based on consensus SV event  
        # a) Single event 
        if len(INS_events) == 1:
            INS = INS_events[0]
            self.isConsensus = True
            self.beg = INS.beg + targetBeg # Map alignments into chromosomal coordinates
            self.end = INS.end + targetBeg
            self.consensusLen = INS.length
            self.insertSeq = INS.seq

        # b) Multiple possible events. Select most likely (TO DO) 
        elif len(INS_events) > 1:
            print("WARNING. Multiple possible INS events")

        # c) No identified events
             
        ## Do Cleanup 
        unix.rm([outDir])

    def determine_insertion_type(self, index, confDict, rootDir): 
        '''
        Determine the type of insertion (retrotransposon, simple repeat, virus, ...) and collect insertion information.

        Input: 
            1. index: Minimap2 index for fasta file containing retrotransposon related sequences 
            2. confDict: Configuration dictionary 
            3. rootDir: Root directory where temporary folders and files will be created
        '''

        ## 0. Create output directory ##
        outDir = rootDir + '/' + str(self.id)
        unix.mkdir(outDir)

        ## 1. Define target INS sequence ##
        # a) INS seq already available
        if self.insertSeq != None:
            seq = self.insertSeq

        # b) INS seq not available -> pick arbitrary INS sequence from one read 
        else:
            self.isConsensus = False
            self.insertSeq = self.events[0].seq
            seq = self.insertSeq
        
        ## 2. Write target seq into file ##
        FASTA = formats.FASTA()
        FASTA.seqDict[str(self.id)] = seq

        FASTA_file = outDir + '/insert.fa'
        FASTA.write(FASTA_file)

        ## 3. Is an insertion of a simple/tandem repeat? ##
        minPercSimple = 70
        self.insType, self.status, self.percResolved = repeats.is_simple_repeat(FASTA_file, minPercSimple, outDir)

        ## Stop if insertion classified as simple repeat
        if self.status == 'resolved':
            return

        ## 4. Is a retrotransposition insertion? ##
        self.insType, self.family, self.srcId, self.status, self.percResolved, self.strand, self.hits = retrotransposons.is_retrotransposition(FASTA_file, index, outDir)
        
        ## Stop if insertion classified as retrotransposon
        if (self.status == 'resolved') or (self.status == 'partially_resolved'):
            return
       
        ## 5. Is a virus? ##
        #viruses.is_virus()

        ## 6. Is a chromosomal insertion? (mitochondrial, templated...) ##
        #¿¿¿rearrangements???.is_chromosomal_dna()

class DEL_cluster(cluster):
    '''
    Deletion (DEL) cluster subclass
    '''
    def __init__(self, events, clusterType):

        cluster.__init__(self, events, clusterType)

    def create_consensus(self, supportingReads, reference, confDict, rootDir): 
        '''
        Use all the cluster supporting reads to (1) generate a consensus sequence for the cluster and (2) perform local realignment of the consensus 
        to obtain more accurate SV breakpoints and inserted sequence (only for INS)

        Update cluster with all this information
        
        Input:
            1. supportingReads: FASTQ (if qualities available) or FASTA (if qualities NOT available) object containing all the reads supporting SV clusters 
            2. reference: Path to the reference genome in fasta format. An index of the reference generated with samtools faidx must be located in the same directory
            3. confDict: Configuration dictionary 
            4. rootDir: Root directory where temporary folders and files will be created
        '''
        ## 0. Create output directory
        outDir = rootDir + '/' + str(self.id)
        unix.mkdir(outDir)

        ## 1. Collect cluster supporting reads from the FASTA/FASTQ containing all the reads supporting SV clusters
        clusterSupportingReads = self.collect_supportingReads(supportingReads, confDict['quality'])

        ## 2. Generate consensus sequence for the cluster
        self.consensus_FASTA = consensus.racon(clusterSupportingReads, confDict['technology'], confDict['quality'], outDir)

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
        INS_events, DEL_events, CLIPPING_left_events, CLIPPING_right_events, supportingReads = bamtools.collectSV(targetInterval, 0, targetLen, BAM, confDict, None)
               
        ## 5. Redefine cluster properties based on consensus SV event  
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

        ## Do Cleanup 
        unix.rm([outDir])

class CLIPPING_cluster(cluster):
    '''
    Clipping cluster subclass
    '''
    def __init__(self, events, clusterType):

        cluster.__init__(self, events, clusterType)