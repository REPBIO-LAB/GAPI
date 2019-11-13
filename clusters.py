'''
Module 'clusters' - Contains classes for dealing with genomic variation
'''

## DEPENDENCIES ##
# External
import sys
import multiprocessing as mp
import numpy as np
import pysam
import math
import itertools
import os
from operator import itemgetter

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
import annotation
import sequences
import gRanges

###############
## FUNCTIONS ##
###############

def create_cluster(events, clusterType):
    '''
    Function to create a cluster object instance

    Input:
        1. events: list of events/clusters that will compose the cluster
        2. clusterType: type of cluster (META: metacluster; INS: insertion; DEL: deletion; CLIPPING: clipping; DISCORDANT: discordant paired-end)

    Output:
        1. cluster: cluster object instance
    '''
    cluster = ''

    ## a) Create META cluster
    # Note: events should correspond to a list of clusters 
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
    Merge a set of clusters/metaclusters into a single cluster/metacluster instance

    Input:
        1. clusters: list of clusters/metaclusters that will be merged
        2. clusterType: type of cluster to be merged (INS: insertion; DEL: deletion; CLIPPING: clipping; META: metacluster)

    Output:
        1. cluster: merged cluster/metacluster instance
    '''
    # A) Merge metaclusters
    if clusterType == 'META':
        subclusters = []
         
        for metacluster in clusters:
            subclusters = subclusters + list(metacluster.subclusters.values())

        mergedCluster = create_cluster(subclusters, clusterType)

    # B) Merge standard clusters
    else:

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
    
    NOTE: this function should be removed and integrated into the more generic one: 'create_clusters'
    '''
    discordantClustersDict = {}

    # For each discordant read pair cluster type
    for clusterType in discordantBinDb.eventTypes:
        
        # Do clustering based on reciprocal overlap
        discordantClustersDict[clusterType] = clustering.reciprocal_overlap_clustering(discordantBinDb, 1, minClusterSize, [clusterType], buffer, clusterType)

    return discordantClustersDict


def extra_clustering_by_matePos(discordantClusters, refLengths, minClusterSize):
    '''
    Apply an extra clustering step to discordant read pair clusters based on mate position

    Input:
        1. discordantClusters: dictionary containing for each possible discordant cluster type (keys) a list of clusters (values)
        2. refLengths: dictionary containing references as keys and their lengths as values
        3. minClusterSize: minimum cluster size

    Output:

    '''
    print('INPUT_extra_clustering_by_matePos: ', discordantClusters, refLengths, minClusterSize)

    ## For each discordant cluster type
    for clusterType, clusters in discordantClusters.items():
        print('CLUSTER_TYPE: ', clusterType, clusters)
        
        ## For each cluster
        for cluster in clusters:

            print('CLUSTER: ', cluster, len(cluster.events), cluster.events)
            cluster_discordants_by_matePos(cluster.events, refLengths, minClusterSize)

            print('--------------')

def cluster_discordants_by_matePos(discordants, refLengths, minClusterSize):
    '''
    Cluster discordant read pairs based on their mate alignment position

    Input:
        1. discordants: list of discordant events
        2. refLengths: dictionary containing references as keys and their lengths as values
        3. minClusterSize: minimum cluster size

    Output:
        1. discordantClusters: list of discordant clusters

    '''   
    print('INPUT_cluster_discordants_by_matePos: ', discordants, refLengths, minClusterSize)

    ## 1. Organize discordant into a dictionary according to supporting read id
    discordantsDict = {}

    for discordant in discordants:

        ## Note: create method to return readName + mateId (\1 and \2). Use this id as dictionary key. 
        # Othewise 
        discordantsDict[discordant.fullReadName()] = discordant

    print('discordantsDict: ', discordantsDict)
    
    ## 2. Produce discordant objects for mates:
    mates = events.discordants2mates(discordants)

    ## 3. Organize mates into a bin database prior clustering
    ## Create dictionary containing mates
    matesDict = events.events2dict(mates, 'DISCORDANT_MATE')

    ## Create bin database 
    matesBinDb = structures.create_bin_database(refLengths, matesDict)

    ## 4. Cluster mates according to their alignment positions
    mateClusters = []

    # For each reference
    for ref in matesBinDb:
        binDb = matesBinDb[ref]
        binLevel = binDb.binSizes[0]
        clusters = clustering.distance_clustering(binDb, binLevel, ['DISCORDANT_MATE'], 'DISCORDANT', binLevel, 1)
        mateClusters = mateClusters + clusters
    
    ## 5. Make clusters of discordants based on mate clusters
    discordantClusters = []

    # For each mate cluster in the reference
    for mateCluster in mateClusters:

        # Create cluster composed by the discordant read pairs corresponding to the mates
        discordants = [discordantsDict[mate.fullReadName_mate()] for mate in mateCluster.events]

        # Add group to the list
        #discordantClusters = discordantClusters + .append(discordants)

    return discordantClusters
    

def create_clusters(eventsBinDb, confDict):
    '''
    Group SV events into distinct clusters according to their SV type
    
    Input:
        1. eventsBinDb: Data structure containing a set of events organized in genomic bins
        2. confDict: 
            * maxInsDist: Maximum distance bewteen two adjacent INS to be clustered together
            * maxBkpDist: Maximum distance bewteen two adjacent breakpoints for CLIPPING clustering    
            * minClusterSize: minimum number of reads composing a root cluster
            * minPercOverlap: minimum percentage of reciprocal overlap for DEL clustering 

    Output:
        1. clustersBinDb: bin database structure containing SV clusters
    ''' 

    ## 1. Create clusters ##
    clustersDict = {}

    ## For each event type, group events into clusters 
    for SV_type in eventsBinDb.eventTypes:

        ## A) Perfom clustering for INS events
        if SV_type == 'INS':

            binLevel = eventsBinDb.binSizes[0]
            clustersDict[SV_type] = clustering.distance_clustering(eventsBinDb, binLevel, [SV_type], SV_type, confDict['maxInsDist'], confDict['minClusterSize'])      

        ## B) Perform clustering for CLIPPING events
        elif SV_type in ['CLIPPING', 'RIGHT-CLIPPING', 'LEFT-CLIPPING']:

            binLevel = eventsBinDb.binSizes[0]
            clustersDict[SV_type] = clustering.distance_clustering(eventsBinDb, binLevel, [SV_type], SV_type, confDict['maxBkpDist'], confDict['minClusterSize'])      

        ## C) Perform clustering based on reciprocal overlap for DEL 
        elif SV_type == 'DEL':

            clustersDict[SV_type] = clustering.reciprocal_overlap_clustering(eventsBinDb, 20, confDict['minClusterSize'], [SV_type], 0, SV_type)

        ## D) Perform clustering based on reciprocal overlap for DISCORDANT
        elif SV_type == 'DISCORDANT':

            clustersDict[SV_type] = clustering.reciprocal_overlap_clustering(eventsBinDb, 1, confDict['minClusterSize'], [SV_type], 100, SV_type)    

    ## 2. Organize clusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    clustersBinDb = structures.create_bin_database_interval(eventsBinDb.ref, eventsBinDb.beg, eventsBinDb.end, clustersDict, binSizes)

    return clustersBinDb


def polish_clusters(clustersBinDb, minClusterSize):
    '''
    Apply set of steps for refining clusters of SV events

    Input:
        1. clustersBinDb: Data structure containing a set of clusters organized in genomic bins  
        2. minClusterSize: minimum cluster size

    Output: 
        It does not return variables as output. Just modify the bin database
    '''
    ## 1. Polish INS clusters
    if 'INS' in clustersBinDb.eventTypes:
        
        # Collect all clusters
        INS_clusters = clustersBinDb.collect(['INS'])

        # Correct INS fragmentation
        for cluster in INS_clusters:    

            cluster.correct_fragmentation()

        # Search for subclusters and remove outliers
        for cluster in INS_clusters:

            subclusters = cluster.identify_subclusters(minClusterSize)

            # Replace cluster by newly created subclusters in the bin database
            if subclusters:

                clustersBinDb.remove([cluster], 'INS')
                clustersBinDb.add(subclusters, 'INS')

    ## 2. Polish DEL clusters (TO DO) 

    ## 3. Polish CLIPPING clusters (TO DO)


def create_metaclusters(clustersBinDb):
    '''    
    Group SV cluster events into metaclusters

    Input:
        1. clustersBinDb: Data structure containing a set of clusters organized in genomic bins  

    Output:
        1. metaclusters: list containing newly created metaclusters
    '''
    metaclusters = clustering.reciprocal_overlap_clustering(clustersBinDb, 1, 1, clustersBinDb.eventTypes, 50, 'META')

    return metaclusters


def SV_type_metaclusters(metaclusters, minINDELlen, technology, rootOutDir):
    '''
    Infer the SV type supported by each metacluster

    Input:
        1. metaclusters: list of metaclusters
        2. minINDELlen: minimum INS and DEL lenght
        3. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        4. rootOutDir: root output directory

    Output:
        1. metaclustersSVType: dictionary containing one key per SV type and the list of metaclusters identified as value
    '''
    metaclustersSVType = {}

    for metacluster in metaclusters:

        metaInterval = '_'.join([str(metacluster.ref), str(metacluster.beg), str(metacluster.end)])
        outDir = rootOutDir + '/' + metaInterval
        metacluster.determine_SV_type(minINDELlen, technology, outDir)

        # A) Initialize list containing metaclusters of a given SV type
        if metacluster.SV_type not in metaclustersSVType:
            metaclustersSVType[metacluster.SV_type] = [metacluster]

        # B) Add metacluster to the list        
        else:
            metaclustersSVType[metacluster.SV_type].append(metacluster)    

    return metaclustersSVType


def create_consensus(metaclusters, confDict, reference, targetSV, rootOutDir):
    '''
    Generate consensus events for a set of metacluster objects

    Input:
        1. metaclusters: Dictionary containing one key per SV type and the list of metaclusters identified as value
        2. confDict: 
            * technology     -> sequencing technology (NANOPORE, PACBIO or ILLUMINA)
            * rounds         -> number of polishing rounds to be attempled. 0 means no polishing
            * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
            * minMAPQ        -> minimum mapping quality
            * minCLIPPINGlen -> minimum clipping lenght
            * minINDELlen    -> minimum INS and DEL lenght
            * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)

        3. reference: Path to reference genome in fasta format    
        4. targetSV: Target SV types to generate consensus
        5. rootOutDir: Root output directory
    ''' 

    ## For each type of SV 
    for SV in targetSV:

        ## Abort if no metacluster from this SV type has been identified
        if SV not in metaclusters:
            continue

        ## For each metacluster
        for metacluster in metaclusters[SV]:
            metaInterval = '_'.join([str(metacluster.ref), str(metacluster.beg), str(metacluster.end)])
            outDir = rootOutDir + '/' + metaInterval

            ## 1. Polish metacluster´s consensus sequence
            metacluster.polish(confDict, reference, outDir)

            ## 2. Obtain consensus metacluster´s event
            metacluster.consensus_event(confDict, reference, 10000, outDir)

            ## Cleanup
            unix.rm([outDir])
    

def lighten_up_metaclusters(metaclusters):
    '''
    Make metacluster objects lighter by removing events and subcluster objects.

    Collect relevant info for downstream analysis as metacluster attributes

    Input:
        1. metaclusters: Dictionary containing one key per SV type and the list of metaclusters identified as value
    '''
    # For each SV type
    for SV_type in metaclusters:

        # For each metacluster
        for metacluster in metaclusters[SV_type]:

            ## Set some object attributes before lightening up
            metacluster.nbTotal, metacluster.nbTumour, metacluster.nbNormal, metacluster.nbINS, metacluster.nbDEL, metacluster.nbCLIPPING = metacluster.nbEvents()

            if 'INS' in metacluster.subclusters:
                metacluster.cv = metacluster.subclusters['INS'].cv_len()[1]

            # Remove events and subclusters from metacluster instance 
            metacluster.events = None
            metacluster.subclusters = None


def double_clipping_supports_INS(clusterA, clusterB, minINDELlen, technology, outDir):
    '''
    Assess if two clipping clusters support an insertion event.

    Input:
        1. clusterA: clipping cluster A
        2. clusterB: clipping cluster B
        3. minINDELlen: minimum INS and DEL lenght
        4. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        5. outDir: output file
    
    Output:
        1. boolean: clippings support an INS event (True) or not (False) 
        2. consensusFasta: fasta containing consensus read sequence spanning the insertion event. None if no insertion supporting evidences found 
    '''
    ## Search for chimeric alignments completely spanning the INS fragment
    primary, supplementary, chimeric, percSameStrand = find_chimeric_alignments(clusterA, clusterB)

    ## A) Chimeric alignment found with both pieces aligning in the same orientation
    if (primary is not None) and (percSameStrand >= 75):

        ## Search for inserted sequence at clipped clusters breakpoints
        insert = find_insertion_at_clipping_bkp(primary, supplementary)

        ## a) Inserted sequence longer than threshold
        if len(insert) >= minINDELlen:
            boolean = True            
            consensusFasta = formats.FASTA()
            consensusFasta.seqDict[primary.readName] = primary.readSeq
            
        ## b) No inserted sequence or shorter than threshold
        else:
            boolean = False
            consensusFasta = None

    ## B) Chimeric alignment NOT found -> search for complementary clippings
    else:

        ## Generate fasta files containing clusters supporting reads:
        readsA = clusterA.collect_reads() 
        readsB = clusterB.collect_reads()

        readsA_file = outDir + '/seqA.fa'
        readsB_file = outDir + '/seqB.fa'

        readsA.write(readsA_file)
        readsB.write(readsB_file)

        ## Assemble clippings based on overlap 
        contigFile = assembly.assemble_overlap(readsA_file, readsB_file, technology, outDir)
        
        ## a) Reciprocal overlap found 
        if contigFile is not None:

            boolean = True
            
            ## Read fasta with contig
            consensusFasta = formats.FASTA()
            consensusFasta.read(contigFile)             

        ## b) Reciprocal overlap not found
        else:
            boolean = False
            consensusFasta = None

    return boolean, consensusFasta

def find_chimeric_alignments(clusterA, clusterB):
    '''
    Search for chimeric read alignments connecting two clipping clusters. Select one as consensus if multiple are identified. 

    Input:
        1. clusterA: Clipping cluster object
        2. clusterB: Clipping cluster object

    Output:
        1. primary: clipping event for representative primary alignment
        2. supplementary: clipping event for representative supplementary alignment
        3. chimericSorted: list of tuples. Each tuple is composed by two clipping events corresponding to primary and supplementary alignments, respectively. 
        4. percSameStrand: percentage of chimeric alignments with both fragments having the same orientation
    ''' 
    ### 1. Identify chimeric read alignments connecting both clipping clusters
    chimeric = []
    sameStrand = 0
    oppositeStrand = 0

    # For each clipping event composing cluster A        
    for clippingA in clusterA.events:

        # For each clipping event composing cluster B
        for clippingB in clusterB.events:

            # Clipping events supported by the same read (chimeric alignments)
            if (clippingA.readName == clippingB.readName):
                
                ## Determine if fragments in chimeric alignments have the same or opposite orientations
                # A) Same
                if (clippingA.reverse == clippingB.reverse):
                    sameStrand += 1

                # B) Opposite
                else:
                    oppositeStrand += 1

                ## Determine which clipping is primary and supplementary
                # A) Clipping A primary; clipping B supplementary
                if (clippingA.supplementary == False) and (clippingB.supplementary == True):
                    primary = clippingA
                    supplementary = clippingB
                    chimeric.append((primary, supplementary))

                # B) Clipping A supplementary; clipping B primary
                elif (clippingA.supplementary == True) and (clippingB.supplementary == False):
                    primary = clippingB
                    supplementary = clippingA
                    chimeric.append((primary, supplementary))

                # C) Both clippings supplementary (Discard! These cases are not informative as supplementary alignments are hardclipped so don´t allow to pick the inserted fragment)
                
    # Exit if not chimeric alignments found
    if not chimeric:
        return None, None, None, None

    ### 2. Select the clipping events with the longest supplementary alignment as representative
    chimericSorted = sorted(chimeric, key=lambda alignments: alignments[1].refLen, reverse=True)
    primary = chimericSorted[0][0]
    supplementary = chimericSorted[0][1]

    ### 3. Compute fraction of chimeric alignments with same orientation
    percSameStrand = float(sameStrand) / (sameStrand + oppositeStrand) * 100

    return primary, supplementary, chimericSorted, percSameStrand


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


def INS_type_metaclusters(metaclusters, reference, refLengths, refDir, transductionSearch, processes, rootOutDir):
    '''
    For each metacluster provided as input determine the type of insertion

    Input:
        1. metaclusters: list of metaclusters supporting insertion events
        2. reference: Path to the reference genome in fasta format (bwa mem and minimap2 indexes must be located in the same folder)
        3. refLengths: Dictionary containing reference ids as keys and as values the length for each reference. 
        4. refDir: Directory containing reference databases. 
        5. transductionSearch: boolean specifying if transduction search is enabled (True) or not (False)
        6. processes: number of processes
        7. rootOutDir: Root output directory
    '''      
    ### 1. Load repeats, transduced regions and exons database 
    msg = '1. Load repeats, transduced regions and exons database'
    log.subHeader(msg)        
    annotDir = rootOutDir + '/ANNOT/'
    unix.mkdir(annotDir)
    annotations2load = ['REPEATS']

    if transductionSearch:    
        annotations2load.append('TRANSDUCTIONS')

    if True: # at one point include flag for pseudogene search
        annotations2load.append('EXONS')

    annotations = annotation.load_annotations(annotations2load, refLengths, refDir, processes, annotDir)

    ## Cleanup
    unix.rm([annotDir])

    ## 2. Create fasta containing all consensus inserted sequences 
    msg = '2. Create fasta containing all consensus inserted sequences'
    log.info(msg)   
    fastaPath = insertedSeq2fasta(metaclusters, rootOutDir)

    ## 3. Align consensus inserted sequences
    msg = '3. Align consensus inserted sequences'
    log.info(msg) 

    ## 3.1 Align consensus inserted sequences into the reference genome
    msg = '3.1 Align consensus inserted sequences into the reference genome'
    log.info(msg)    
    SAM_genome = alignment.alignment_bwa(fastaPath, reference, 'alignments_genome', processes, rootOutDir)

    ## Convert SAM to PAF
    PAF_genome = alignment.sam2paf(SAM_genome, 'alignments_genome', rootOutDir)

    ## Organize hits according to their corresponding metacluster
    allHits_genome = alignment.organize_hits_paf(PAF_genome) 

    ## 3.2 Align consensus inserted sequences into the reference genome (splicing-aware)
    msg = '3.2 Align consensus inserted sequences into the reference genome (splicing-aware)'
    log.info(msg)    

    ## Minimap index for the reference
    index = os.path.splitext(reference)[0] + '.mmi'
    SAM_splicing = alignment.alignment_minimap2_spliced(fastaPath, index, 'alignments_spliced', processes, rootOutDir)

    ## Convert SAM to BAM
    BAM_splicing = bamtools.SAM2BAM(SAM_splicing, rootOutDir)

    ## Convert BAM to BED
    BED_path = bamtools.BAM2BED(BAM_splicing, rootOutDir)

    ## Organize hits according to their corresponding metacluster
    allHits_splicing = formats.BED()
    allHits_splicing.read(BED_path, 'List', None)
    groupedEntries = allHits_splicing.group_entries_by_name()

    ## 3.3 Align consensus inserted sequences into the viral database
    msg = '3.3 Align consensus inserted sequences into the viral database'
    log.info(msg)    

    '''
    SAM_viral = alignment.alignment_bwa(fastaPath, reference, processes, rootOutDir)

    ## Convert SAM to PAF
    PAF_viral = alignment.sam2paf(SAM_viral, 'alignments_viral', rootOutDir)

    ## Organize hits according to their corresponding metacluster
    allHits_viral = organize_hits_paf(PAF_viral) 
    '''
    allHits_viral = {}

    ## 4. For each metacluster determine the insertion type
    msg = '4. For each metacluster determine the insertion type'
    log.info(msg)   
    
    # For each metacluster
    for metacluster in metaclusters:

        ## 4.1 Collect consensus inserted sequence hits
        metaId = str(metacluster.ref) + ':' + str(metacluster.beg) + '-' + str(metacluster.end)

        ## Hits in the reference genome
        if metaId in allHits_genome:
            hits_genome = allHits_genome[metaId]
        
        else:
            hits_genome = formats.PAF()

        ## Hits in the reference genome (splice-aware alignment)
        if metaId in groupedEntries:
            hits_splicing = groupedEntries[metaId]

        else:
            hits_splicing = []

        ## Hits in the viral database
        if metaId in allHits_viral:
            hits_viral = allHits_viral[metaId]

        else:
            hits_viral = formats.PAF()

        ## 4.2 Insertion type inference
        metacluster.determine_INS_type(hits_genome, hits_splicing, hits_viral, annotations['REPEATS'], annotations['TRANSDUCTIONS'], annotations['EXONS'])


def structure_inference_parallel(metaclusters, consensusPath, transducedPath, transductionSearch, processes, rootDir):
    '''
    Infer structure for a list of INS metacluster objects. Parallelize by distributing metaclusters by processes. 

    Input:
        1. metaclusters: list of metacluster objects
        2. consensusPath: path to fasta file containing retrotransposon consensus sequences
        3. transducedPath: path to fasta containing transduced sequences downstream of source elements
        4. transductionSearch: boolean specifying if transduction search is enabled (True) or not (False)
        5. processes: number of processes
        6. rootDir: Root output directory
    
    Output:
        1. metaclustersOut: list of metacluster objects with structure information stored at 'SV_features' dict attribute
    '''
    ## 1. Create tuple list for multiprocessing
    msg = '1. Create tuple list for multiprocessing'
    log.subHeader(msg)      
    metaclustersOut = []
    tupleList = []

    for metacluster in metaclusters:
        
        ## Skip structure inference if insertion type not available or not solo, partnered or orphan transduction
        # Note: investigate why INS_TYPE is not defined in some metaclusters
        if ('INS_TYPE' not in metacluster.SV_features) or (metacluster.SV_features['INS_TYPE'] not in ['solo', 'partnered', 'orphan']):
            metaclustersOut.append(metacluster) 
            continue

        ## Add to the list
        fields = (metacluster, consensusPath, transducedPath, transductionSearch, rootDir)
        tupleList.append(fields)

    ## 2. Infer structure
    msg = '2. Infer structure'
    log.subHeader(msg)       

    pool = mp.Pool(processes=processes)
    metaclustersStructure = pool.starmap(structure_inference, tupleList)
    pool.close()
    pool.join()
      
    ## 3. Create final list with metaclusters to be reported as output
    msg = '3. Create final list with metaclusters to be reported as output'
    metaclustersOut = metaclustersOut + metaclustersStructure 

    return metaclustersOut

def structure_inference(metacluster, consensusPath, transducedPath, transductionSearch, rootDir):
    '''
    Wrapper to call 'determine_INS_structure' method for a given INS metacluster provided as input

    Input:
        1. metacluster: INS metacluster 
        2. consensusPath: path to fasta file containing retrotransposon consensus sequences
        3. transducedPath: path to fasta containing transduced sequences downstream of source elements
        4. transductionSearch: boolean specifying if transduction search is enabled (True) or not (False)
        5. rootDir: Root output directory
    
    Output:
        1. metacluster: INS metacluster containing structural properties  
    '''
    # Create output directory
    metaInterval = '_'.join([str(metacluster.ref), str(metacluster.beg), str(metacluster.end)])
    outDir = rootDir + '/' + metaInterval
    unix.mkdir(outDir)

    # Infer structure
    structure = metacluster.determine_INS_structure(consensusPath, transducedPath, transductionSearch, outDir)

    # Add structure info to the metacluster
    metacluster.SV_features.update(structure) 
    
    # Remove output directory
    unix.rm([outDir])

    return metacluster

def insertedSeq2fasta(metaclusters, outDir):
    '''
    Collect all the consensus inserted sequences from a list of metaclusters supporting INS and 
    generate fasta file containing them 

    Input:
        1. metaclusters: list of metaclusters
        2. outDir: output directory

    Output:
        1. fastaPath: fasta file containing all the consensus inserted sequences
    '''
    ## Initiate FASTA object
    FASTA = formats.FASTA()

    # For each metacluster
    for metacluster in metaclusters:

        ## Skip insertion type inference if consensus event not available
        if metacluster.consensusEvent is None:
            continue               

        ## Retrieve inserted sequence and add to the FASTA
        metaclusterId = metacluster.ref + ':' + str(metacluster.beg) + '-' + str(metacluster.end)
        insert = metacluster.consensusEvent.pick_insert()

        FASTA.seqDict[metaclusterId] = insert
        
    ## Write fasta         
    fastaPath = outDir + '/inserted_sequences.fa'
    FASTA.write(fastaPath)    

    return fastaPath


def assignAligments2metaclusters_paf(metaclusters, PAF_path):
    '''
    Map alignments to their corresponding metacluster. 

    Input:
        1. metaclusters: list of metaclusters
        2. PAF_path: Path to path file containing alignments to asign

    Output:
        1. tupleList: List of tuples. Each tuple contain two elements: 
            1) metacluster object  
            2) PAF object containing the corresponding alignments for the metacluster consensus inserted sequence
    '''
    ## 1. Create a dictionary to organize the data
    hits = {}

    for metacluster in metaclusters:
        metaclusterId = metacluster.ref + ':' + str(metacluster.beg) + '-' + str(metacluster.end)
        hits[metaclusterId] = (metacluster, formats.PAF())

    ## 2. Read PAF file and add hits to the metaclusters
    ## Read PAF 
    PAF = formats.PAF()
    PAF.read(PAF_path)

    # For each read alignment 
    for alignment in PAF.alignments:
        hits[alignment.qName][1].alignments.append(alignment)
    
    ## 3. Generate list of tuples
    tupleList = list(hits.values())

    return tupleList


def assignAligments2metaclusters_sam(metaclusters, SAM_path):
    '''
    Map alignments to their corresponding metacluster. 

    Input:
        1. metaclusters: list of metaclusters
        2. SAM_path: Path to SAM file containing alignments to asign

    Output:
        1. metaclustersHits: Nested list. Each list element is composed by a list with two elements: 
            1) metacluster object  
            2) List of aligned segment objects
    '''
    ## 1. Create a dictionary to organize the data
    hits = {}

    for metacluster in metaclusters:
        metaclusterId = metacluster.ref + ':' + str(metacluster.beg) + '-' + str(metacluster.end)
        hits[metaclusterId] = [metacluster, []]

    ## 2. Read SAM file and add hits to the metaclusters
    ## Read SAM 
    SAM = pysam.AlignmentFile(SAM_path, "r")

    # For each read alignment 
    for alignment in SAM:
        hits[alignment.query_name][1].append(alignment)
    
    ## 3. Generate nested list
    metaclustersHits = list(hits.values())

    return metaclustersHits

def search4bridges_metaclusters(metaclusters, maxBridgeLen, minMatchPerc, refLengths, refDir, processes, outDir):
    '''
    Search for transduction or repeat bridges at BND junctions for a list of metacluster objects

    Input:
        1. metaclusters: list of input metacluster objects supporting BND
        2. maxBridgeLen: maximum supplementary cluster length to search for a bridge
        3. minMatchPerc: minimum percentage of the supplementary cluster interval to match in a transduction or repeats database to make a bridge call
        4. refLengths: dictionary containing reference ids as keys and as values the length for each reference. 
        5. refDir: directory containing reference databases. 
        6. processes: number of processes
        7. outDir: output directory

    Output: For each metacluster update... attribute
    '''
    
    ## 1. Load repeats annnotation and transduced regions beds
    annot2load = ['REPEATS', 'TRANSDUCTIONS']
    annotations = annotation.load_annotations(annot2load, refLengths, refDir, processes, outDir)

    ## 2. Search for bridges
    # For each metacluster
    for metacluster in metaclusters:

        metacluster.search4bridges(maxBridgeLen, minMatchPerc, annotations)


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
        self.id = 'CLUSTER_' + str(cluster.number)
        self.clusterId = None

        # Define list of events composing the cluster and cluster type
        self.events = events
        self.clusterType = clusterType

        # Set cluster's reference, begin and end position
        self.ref, self.beg, self.end = self.coordinates() 

        # Cluster filtering
        self.filters = None
        self.nbOutliers = 0

        # Update event's clusterId attribute
        for event in events:
            event.clusterId = self.id        

    def length(self):
        '''
        Compute cluster interval length
        '''
        return self.end - self.beg
        
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

    def add(self, events2add):
        '''
        Incorporate events into the cluster and redefine cluster beg and end
        positions accordingly

        Input:
            1. events2add: List of events to be added to the cluster
        '''
        # Add events to the cluster  
        self.events = self.events + events2add

        # Resort and redefine cluster begin and end coordinates
        self.ref, self.beg, self.end = self.coordinates() 

        # Update event's clusterId attribute
        for event in events2add:
            event.clusterId = self.id

    def remove(self, events2remove):
        '''
        Remove list of events from the cluster and redefine cluster beg and end
        positions accordingly

        Input:
            1. events2remove: List of events to be removed 
        '''
        ## 1. Remove events from the metacluster ##
        self.events = [event for event in self.events if event not in events2remove]

        ## 2. Resort and redefine metacluster begin and end coordinates ##
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

    def cv_len(self):
        '''
        Compute mean length and coefficient of variation (cv)

        Output:
            1. meanLen: Mean length
            2. cv: Coefficient of variation 
        '''
        lengths = [ event.length for event in self.events]
        meanLen = np.mean(lengths)
        std = np.std(lengths)
        cv = std / meanLen * 100 

        return meanLen, cv

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

    def correct_fragmentation(self):
        '''
        Correct fragmented alignments over INS
        
        Before merging:
        ############<<<INS>>>##<<INS>>###<<<INS>>>##############
        
        After merging:
        ############<<<<<<<<INS>>>>>>>>##############        
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

                ## 2.2 Add merged INS
                merged_list.append(merged)

                ## 2.3 Update fragmented alignments list
                fragmented_list = fragmented_list + INS_list
    
        ## 3. Update cluster
        self.add(merged_list)
        self.remove(fragmented_list)

    def identify_subclusters(self, minClusterSize):
        '''
        Identify subclusters of INS with similar lengths

        Input:
            1. minClusterSize: minimum size to create a subcluster

        Output:
            1. subclusters: list of subclusters
        '''

        subclusters = []

        ## Compute metrics based on events length
        meanLen, cv = self.cv_len()

        ## A) Cluster with heterogeneous lengths  
        if cv > 15:

            ## 1. Cluster events according to their length
            clustered = clustering.KMeans_clustering(self.events, None, 'length')

            ## 2. Filter out subclusters supported by less than X events
            clusteredFiltered = []

            for cluster in clustered.values():

                if len(cluster) >= minClusterSize:
                    clusteredFiltered.append(cluster)

            ## 3. Create subclusters
            # NOTE: Don´t create subclusters if:
            #   a) Cluster no fragmented into any subcluster passing support filter OR
            #   b) Cluster fragmented in more than 2 subclusters passing support filter (indicative of noisy region)
            if len(clusteredFiltered) in [1, 2]:

                for events in clusteredFiltered:

                    subcluster = INS_cluster(events)
                    subclusters.append(subcluster)

        return subclusters     

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

class SUPPLEMENTARY_cluster(cluster):
    '''
    Supplementary alignment cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'SUPPLEMENTARY')
        self.bkpSide = None
        self.annot = None


class DISCORDANT_cluster(cluster):
    '''
    Discordant cluster subclass
    '''
    def __init__(self, events):

        cluster.__init__(self, events, 'DISCORDANT')
        self.matesCluster = None


class META_cluster():
    '''
    Meta cluster class
    '''
    number = 0 # Number of instances

    def __init__(self, clusters):
        '''
        '''
        META_cluster.number += 1 # Update instances counter
        self.id = 'META_' + str(META_cluster.number)

        # Define list of events composing the cluster 
        self.events = list(itertools.chain(*[cluster.events for cluster in clusters]))

        # Set cluster's reference, begin and end position
        self.ref, self.beg, self.end = self.coordinates() 

        # Organize events into subclusters
        self.subclusters = self.create_subclusters()

        # Set some metacluster properties as None
        self.mutOrigin = None
        self.failedFilters = None
        self.consensusEvent = None                
        self.insertHits = None
        self.nbTotal, self.nbTumour, self.nbNormal, self.nbINS, self.nbDEL, self.nbCLIPPING = [None, None, None, None, None, None] 
        self.cv = None
        self.bridgeClusters = {}
        self.bridgeType = None

        # Update input cluster's clusterId attribute
        for cluster in clusters:
            cluster.clusterId = self.id

        # Initialize dictionaries 
        self.SV_features = {}
        self.supplAlignClusters = {}

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

            ## Set subcluster metacluster id attribute:
            subcluster.clusterId = self.id

            ## Add subcluster to the dict
            subclusters[eventType] = subcluster 

        return subclusters

    def add(self, clusters2add):
        '''
        Add a list of clusters to the metacluster 

        Input:
            1. clusters2add: List of clusters to be added 
        '''
        ## 0. Update metacluster's id attribute
        for cluster in clusters2add:
            cluster.clusterId = self.id

        ## 1. Add events within the input cluster to the metacluster ##
        events2add = list(itertools.chain(*[cluster.events for cluster in clusters2add]))
        self.events = self.events + events2add

        ## 2. Resort and redefine metacluster begin and end coordinates ##
        self.ref, self.beg, self.end = self.coordinates()

        ## 3. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(events2add)

        ## 4. Add events to the subclusters ##
        for eventType, eventList in eventTypes.items():
            
            # a) Create subcluster if not pre-existing one
            if eventType not in self.subclusters:
         
                ## Create subcluster
                subcluster = create_cluster(eventList, eventType) 
            
                ## Set subcluster metacluster id attribute:
                subcluster.clusterId = self.id

                ## Add subcluster to the dict
                self.subclusters[eventType] = subcluster 

            # b) Add events to pre-existing subcluster
            else:
                self.subclusters[eventType].add(eventList)

        # Update input cluster's clusterId attribute
        for cluster in clusters2add:
            cluster.clusterId = self.id

    def remove(self, events2remove):
        '''
        Remove a list of events from the metacluster and corresponding subclusters

        Input:
            1. events2remove: List of events to be removed 
        '''
        ## 1. Remove events from the metacluster ##
        self.events = [event for event in self.events if event not in events2remove]

        ## 2. Resort and redefine metacluster begin and end coordinates ##
        self.ref, self.beg, self.end = self.coordinates()

        ## 3. Separate events according to their type into multiple lists ##
        eventTypes = events.separate(events2remove)

        ## 4. Remove events from the subclusters ##
        for eventType, eventList in eventTypes.items():

            if eventType in self.subclusters:
                self.subclusters[eventType].remove(eventList)
                
            else:
                log.info('WARNING at remove method from META_cluster class. Event with unkown type')

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
 
    def polish(self, confDict, reference, outDir):
        '''
        Polish metacluster consensus sequence

        Input: 
            1. confDict: 
                * technology     -> sequencing technology (NANOPORE, PACBIO or ILLUMINA)
                * rounds         -> number of polishing rounds to be attempled. 0 means no polishing
            
            2. reference: path to reference genome in fasta format    
            3. outDir: Output directory
        
        Output: Update 'consensusFasta' attribute with the polished sequence
        '''
        ## 0. Create directory 
        unix.mkdir(outDir)        

        ## 1. Use cluster supporting reads to polish metacluster consensus sequence ##
        ## 1.1 Write raw consensus sequence
        unpolishedFasta = outDir + '/raw_consensus.fa'
        self.consensusFasta.write(unpolishedFasta)

        ## 1.2 Collect reads that will be used to polish the consensus sequence        
        supportingReads = self.collect_reads()

        ## Remove consensus from FASTA
        consensusReadName = list(self.consensusFasta.seqDict.keys())[0]

        if consensusReadName in supportingReads.seqDict:
            del supportingReads.seqDict[consensusReadName]

        ## Write supporting reads FASTA 
        supportingReadsFasta = outDir + '/supportingReads.fa'
        supportingReads.write(supportingReadsFasta)
            
        ## 2. Consensus polishing 
        polishedFasta = assembly.polish_racon(unpolishedFasta, supportingReadsFasta, confDict['technology'], confDict['rounds'], outDir)

        ## If polishing is successful replace consensus by new polished sequence
        if polishedFasta is not None:
            polished = formats.FASTA()
            polished.read(polishedFasta)     
            self.consensusFasta = polished

        return 

    def consensus_event(self, confDict, reference, offset, outDir):
        '''
        Define metacluster´s consensus event based on consensus sequence realignment

        Input: 
            1. confDict: 
                * targetSV       -> list with target SV (INS: insertion; DEL: deletion; CLIPPING: left and right clippings)
                * minMAPQ        -> minimum mapping quality
                * minCLIPPINGlen -> minimum clipping lenght
                * minINDELlen    -> minimum INS and DEL lenght
                * overhang       -> Number of flanking base pairs around the INDEL events to be collected from the supporting read. If 'None' the complete read sequence will be collected)            2. outDir: output directory
            
            2. reference: path to reference genome in fasta format    
            3. offset: number of base pairs to extend cluster begin and end coordinates when defining target region for consensus sequence realignment
            4. outDir: Output directory
        
        Output: Update 'consensusEvent' attribute in the metacluster class
        '''

        ## 1. Local realignment of the consensus sequence into the SV genomic interval ##
        ## 1.1 Write consensus sequence into a fasta file 
        consensusFile = outDir + '/consensus.fa'
        self.consensusFasta.write(consensusFile)

        ## 1.2 Define SV cluster genomic interval
        # ------------------<***SV_cluster***>-----------------
        #       <--offset-->                  <--offset-->
        intervalBeg = self.beg - offset
        intervalBeg = intervalBeg if intervalBeg >= 0 else 0 ## Set lower bound
        intervalEnd = self.end + offset
        intervalCoord = self.ref + ':' + str(intervalBeg) + '-' + str(intervalEnd)
            
        ## 1.3 Do realignment
        # ------------------<***SV_cluster***>-----------------
        #         -------------consensus_seq-------------
        BAM = alignment.targeted_alignment_minimap2(consensusFile, intervalCoord, reference, outDir)
 
        ## Continue if realignment is succesfull 
        if BAM is not None:
            
            ## 2. Search for metacluster´s consensus event  ##
            ## 2.1 Extract events from consensus sequence realignment
            # ------------------<***SV_cluster***>-----------------
            #        -------------consensus_seq-------------
            #               <--->----------------<---> overhang (100 bp)
            #                   event_search_space         
            overhang = 100
            clusterIntervalLen = self.end - self.beg
            targetBeg = offset - overhang
            targetEnd = offset + clusterIntervalLen + overhang            
            eventsDict = bamtools.collectSV(intervalCoord, targetBeg, targetEnd, BAM, confDict, None)

            ## 2.2 Define consensus event based on the events resulting from consensus sequence realignment
            ## A) Metacluster supports an INS and realignment leads to one INS event 
            if (self.SV_type == 'INS') and ('INS' in eventsDict) and (len(eventsDict['INS']) == 1):

                ## Convert coordinates
                self.consensusEvent = alignment.targetered2genomic_coord(eventsDict['INS'][0], self.ref, intervalBeg)

            ## B) Metacluster supports an INS and realignment leads to multiple INS events
            elif (self.SV_type == 'INS') and ('INS' in eventsDict) and (len(eventsDict['INS']) > 1):

                ## Do merging
                merged = events.merge_INS(eventsDict['INS'])

                ## Convert coordinates
                self.consensusEvent = alignment.targetered2genomic_coord(merged, self.ref, intervalBeg)

            ## C) Metacluster supports a DEL and realignment leads to one DEL event 
            elif (self.SV_type == 'DEL') and ('DEL' in eventsDict) and (len(eventsDict['DEL']) == 1):

                ## Convert coordinates
                self.consensusEvent = alignment.targetered2genomic_coord(eventsDict['DEL'][0], self.ref, intervalBeg)

            ## D) Metacluster supports a DEL and realignment leads to multiple DEL events (TO DO)
            #elif (self.SV_type == 'DEL') and (len(eventsDict['DEL']) > 1):
                    ## Do merging
                    ## Convert coordinates
                        
            ## E) Metacluster supports an INS and realignment leads to one left and one right CLIPPING 
            elif (self.SV_type == 'INS') and ('LEFT-CLIPPING' in eventsDict) and ('RIGHT-CLIPPING' in eventsDict) and (len(eventsDict['RIGHT-CLIPPING']) == 1) and (len(eventsDict['LEFT-CLIPPING']) == 1):

                ## Compute the inserted sequence length as the difference between both clipping breakpoints in the long sequence 
                rightClipping = eventsDict['RIGHT-CLIPPING'][0]
                leftClipping = eventsDict['LEFT-CLIPPING'][0]
                length = leftClipping.readBkp - rightClipping.readBkp

                ## Create consensus INS from clipping alignments if inserted sequence found
                if length >= confDict['minINDELlen']:

                    event = events.INS(rightClipping.ref, rightClipping.beg, rightClipping.end, length, rightClipping.readName, rightClipping.readSeq, rightClipping.readBkp, None, None)
        
                    ## Convert coordinates
                    self.consensusEvent = alignment.targetered2genomic_coord(event, self.ref, intervalBeg)
    
                else:
                    print('INS_NOT_FOUND!')

            ## F) Another possibility (Don´t do anything, leave previous. Later we may need to include new conditions)

    def determine_SV_type(self, minINDELlen, technology, outDir): 
        '''
        Determine the type of structural variant (SV) supported by the metacluster and select a consensus metacluster supporting sequence and event

        SV types:
            INSERTION: 
                        >>>>>>>>>>>>>/////INS/////>>>>>>>>>>>    * Completely spanned insertion
                        >>>>>>>>>>>>>/////INS///                 * Insertion partially spanned
                                      ////INS/////>>>>>>>>>>>
            DELETION:   >>>>>>>>>>>>>-----DEL----->>>>>>>>>>>    

            DUPLICATION (TO DO) 
            INVERSION (TO DO)
            BREAK END (TO DO)      

        Input:
            1. minINDELlen: minimum INS and DEL lenght
            2. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
            3. outDir: Output directory

        Output, set the following object attributes:
            
            1. SV_type: structural variant type supported by the metacluster. None if sv type 
            2. consensusEvent: consensus event attribute 
            3. consensusFasta: fasta file containing consensus metacluster supporting sequence 
        '''
        ## Create output directory 
        unix.mkdir(outDir)
        subClusterTypes = list(self.subclusters.keys())

        ## A) Metacluster supports an insertion:
        if ('INS' in subClusterTypes) and ('DEL' not in subClusterTypes):
            self.SV_type = 'INS'

            ## Select consensus INS event and sequence
            self.consensusEvent = self.subclusters['INS'].pick_median_length()
            
            self.consensusFasta = formats.FASTA()
            self.consensusFasta.seqDict[self.consensusEvent.readName] = self.consensusEvent.readSeq

        ## B) Metacluster supports a deletion:
        elif ('DEL' in subClusterTypes) and ('INS' not in subClusterTypes):
            self.SV_type = 'DEL'

            ## Select consensus DEL event and sequence
            self.consensusEvent = self.subclusters['DEL'].pick_median_length()

            self.consensusFasta = formats.FASTA()
            self.consensusFasta.seqDict[self.consensusEvent.readName] = self.consensusEvent.readSeq
                        
        ## C) Metacluster only composed by double clipping -> Long insertion candidate
        elif (len(subClusterTypes) == 2) and ('RIGHT-CLIPPING' in subClusterTypes) and ('LEFT-CLIPPING' in subClusterTypes):

            self.consensusEvent = None                

            ## Assess if clipping clusters support an insertion
            is_INS, self.consensusFasta = double_clipping_supports_INS(self.subclusters['RIGHT-CLIPPING'], self.subclusters['LEFT-CLIPPING'], minINDELlen, technology, outDir)

            ## a) Double clippings support an INS
            if is_INS:
                self.SV_type = 'INS'

            ## b) Double clippings not support an INS
            else: 
                self.SV_type = None

        ## D) Metacluster only composed by one clipping side (left or right)
        elif (len(subClusterTypes) == 1) and (('RIGHT-CLIPPING' in subClusterTypes) or ('LEFT-CLIPPING' in subClusterTypes)):

            self.SV_type = 'BND'

            ## Cluster supplementary alignment positions
            self.cluster_suppl_positions()

        ## E) Other combination -> Unknown SV type (Temporal, extend later)
        else:
            self.SV_type = None
            self.consensusEvent = None                
            self.consensusFasta = None


    def cluster_suppl_positions(self):
        '''
        Cluster supplementary alignments based on their begin/end alignment positions

        Output:
            Set 'supplAlignClusters' attribute. Dictionary containing references as keys and nested lists of supplementary 
            alignment clusters as values.
        '''

        ## 1. Collect for each clipping event its suppl. alignments
        supplAlignmentsDictList = []
                
        for clipping in self.events:
            supplAlignmentsDict = clipping.parse_supplAlignments_field()
            supplAlignmentsDictList.append(supplAlignmentsDict)
        
        ## 2. Merge suppl. alignments into a single dictionary
        allSupplAlignmentsDict = structures.merge_dictionaries(supplAlignmentsDictList)

        ## 3. Cluster suppl. alignments based on breakpoint position
        clustersDict = {}

        # For each reference
        for ref in allSupplAlignmentsDict:

            # Collect suppl. alignments list
            supplAlignments = allSupplAlignmentsDict[ref]

            ## Cluster suppl. alignments based on their beg and end alignment positions
            clustersBeg = clustering.distance_clustering_targetPos(supplAlignments, 100, 'beg')
            clustersEnd = clustering.distance_clustering_targetPos(supplAlignments, 100, 'end')

            ## Determine bkp side based on the biggest cluster
            biggestLenBeg = max([len(cluster.events) for cluster in clustersBeg])
            biggestLenEnd = max([len(cluster.events) for cluster in clustersEnd])

            # a) Bkp at the beg of supplementary alignment interval
            if biggestLenBeg >= biggestLenEnd:
                clusters = clustersBeg
            
            # b) Bkp at the end of supplementary alignment interval
            else:
                clusters = clustersEnd
                
            ## Add clusters to the dictionary
            self.supplAlignClusters[ref] = clusters

    def search4bridges(self, maxBridgeLen, minMatchPerc, annotations):
        '''        
        Search for transduction or repeat bridges at metacluster BND junction

        Input:
            1. maxBridgeLen: maximum supplementary cluster length to search for a bridge
            2. minMatchPerc: minimum percentage of the supplementary cluster interval to match in a transduction or repeats database to make a bridge call
            3. annotations: dictionary containing one key per type of annotation loaded and bin databases containing annotated features as values (None for those annotations not loaded)

        Output: Update 'bridgeClusters' and 'bridgeType' metacluster attributes
        '''
        ## For each reference        
        for ref in self.supplAlignClusters:

            ## For each suppl. cluster in the reference
            for supplCluster in self.supplAlignClusters[ref]:

                ## Cluster interval length within limit
                if supplCluster.length() <= maxBridgeLen:

                    ### 1. Assess if supplementary cluster spans a transduced region
                    if 'TRANSDUCTIONS' in annotations:
                        
                        ## Select transduction bin database for the corresponding ref            
                        tdBinDb = annotations['TRANSDUCTIONS'][ref]

                        ## Intersect supplementary cluster interval with transducted regions  
                        tdMatches = tdBinDb.collect_interval(supplCluster.beg, supplCluster.end, 'ALL')

                        ## Sort matches in decreasing coordinates order
                        sortedTdMatches = sorted(tdMatches, key=itemgetter(2), reverse=True)

                        ## Check if cluster matches a transduced area
                        if sortedTdMatches and sortedTdMatches[0][2] > minMatchPerc:
                            tdMatch = sortedTdMatches[0]
                        else:
                            tdMatch = None

                        ## Stop if suppl. cluster matches transduced area
                        if tdMatch is not None:
                            supplCluster.annot = tdMatch

                            # a) First suppl. cluster supporting a transduction bridge -> Initialize list
                            if 'TRANSDUCTION' not in self.bridgeClusters:
                                self.bridgeClusters['TRANSDUCTION'] = [supplCluster]

                            # b) Add suppl. cluster to pre-existing list 
                            else:
                                self.bridgeClusters['TRANSDUCTION'].append(supplCluster)                             

                            continue

                    ### 2. Assess if supplementary cluster spans an annotated repeat
                    if 'REPEATS' in annotations:

                        ## Select repeats bin database for the corresponding ref
                        repeatsBinDb = annotations['REPEATS'][ref]

                        ## Intersect supplementary cluster interval with annotated repeats 
                        repeatMatches = repeatsBinDb.collect_interval(supplCluster.beg, supplCluster.end, 'ALL')

                        ## Sort matches in decreasing coordinates order
                        sortedRepeatMatches = sorted(repeatMatches, key=itemgetter(2), reverse=True)

                        ## Check if cluster matches a transduced area
                        if sortedRepeatMatches and sortedRepeatMatches[0][2] > minMatchPerc:
                            repeatMatch = sortedRepeatMatches[0]
                        else:
                            repeatMatch = None

                        ## Cluster matches annotated repeat area
                        if repeatMatch is not None:
                            supplCluster.annot = repeatMatch

                            # a) First suppl. cluster supporting a repeat bridge -> Initialize list
                            if 'REPEAT' not in self.bridgeClusters:
                                self.bridgeClusters['REPEAT'] = [supplCluster]

                            # b) Add suppl. cluster to pre-existing list 
                            else:
                                self.bridgeClusters['REPEAT'].append(supplCluster)    

        ### Determine bridge type:        
        # a) Partnered transduction (Transduction + repeat) 
        if ('TRANSDUCTION' in self.bridgeClusters) and ('REPEAT' in self.bridgeClusters):
            self.bridgeType = 'partnered'

        # b) Orphan transduction
        elif ('TRANSDUCTION' in self.bridgeClusters):
            self.bridgeType = 'orphan'

        # c) Solo repeat
        elif ('REPEAT' in self.bridgeClusters):
            self.bridgeType = 'repeat'

        # d) No bridge
        else:
            self.bridgeType = None

    def determine_INS_type(self, hits_genome, hits_splicing, hits_viral, repeatsDb, transducedDb, exonsDb):
        '''
        Determine the type of insertion based on the alignments of the inserted sequence on the reference genome

        Input:
            1. hits_genome: PAF object containing inserted sequence alignments on the reference genome
            2. hits_splicing: list of bed entries containing inserted sequence alignments on the reference genome (splice-aware alignment). None if not hit found
            3. hits_viral: PAF object containing inserted sequence alignments on the viral database 
            4. repeatsDb: bin database containing annotated repeats in the reference. None if not available
            5. transducedDb: bin database containing regions transduced by source elements. None if not available
            6. exonsDb: bin database containing annotated exons. None if not available. 

        Output: Add INS type annotation to the attribute SV_features
        ''' 
        ## 0. Abort if consensus event not available 
        if self.consensusEvent is None:
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0
            return 

        ## 1. Assess if input sequence corresponds to a repeat expansion
        is_EXPANSION, self.insertHits = self.is_expansion(hits_genome, repeatsDb)

        # Stop if insertion is a expansion
        if is_EXPANSION:
            return

        ## 2. Assess if input sequence corresponds to duplication 
        is_DUP, self.insertHits = self.is_duplication(hits_genome, 100)

        # Stop if insertion is a duplication
        if is_DUP:
            return

        ## 3. Assess if input sequence corresponds to solo interspersed insertion or transduction
        # Note: return boolean as well specifying if interspersed or not
        is_INTERSPERSED, INS_features, self.insertHits = retrotransposons.is_interspersed_ins(self.consensusEvent.pick_insert(), hits_genome, repeatsDb, transducedDb)

        # Update metacluster with insertion features
        self.SV_features.update(INS_features) 

        # Stop if insertion is a interspersed insertion
        if is_INTERSPERSED:
            return    

        ## 4. Assess if input sequence corresponds to processed pseudogene insertion
        is_PSEUDOGENE, outHits = self.is_processed_pseudogene(hits_splicing, exonsDb)

        # Stop if insertion is a processed pseudogene
        if is_PSEUDOGENE:
            return    

        ## 5. Assess if input sequence corresponds to a viral insertion

    def determine_INS_structure(self, consensusPath, transducedPath, transductionSearch, outDir):
        '''
        Infer inserted sequence structural features

        Input:
            1. consensusPath: path to fasta file containing retrotransposon consensus sequences
            2. transducedPath: path to fasta containing transduced sequences downstream of source elements
            3. transductionSearch: boolean specifying if transduction search is enabled (True) or not (False)
            4. outDir: output directory
    
        Output: 
            structure: dictionary containing insertion structural properties
        '''
        ##  Skip structure inference if consensus event not available
        if self.consensusEvent is None:
            return {}

        ## 1. Read fasta files 
        #  1.1 Consensus sequences
        consensus = formats.FASTA()
        consensus.read(consensusPath)

        #  1.2 Transduced regions
        if transductionSearch: 
            transduced = formats.FASTA()
            transduced.read(transducedPath)

        ## 2. Create fasta object containing database of sequences
        ## The database will contain the following sequences depending on the insertion type:
        ## - Solo      -> consensus sequences for the same family
        ## - Partnered -> consensus sequences for the same family
        #              -> corresponding transduced area
        ## - Orphan    -> corresponding transduced area

        ## Initialize fasta
        fasta = formats.FASTA()

        ## Add to the fasta subfamily consensus sequences for the corresponding family
        if self.SV_features['INS_TYPE'] in ['solo', 'partnered']:
            for seqId, seq in consensus.seqDict.items(): 
                family = seqId.split('|')[1]

                if family in self.SV_features['FAMILY']:
                    fasta.seqDict[seqId] = seq

        ## Add to the fasta transduced region or regions
        if self.SV_features['INS_TYPE'] in ['partnered', 'orphan']:
        
            for seqId, seq in transduced.seqDict.items(): 
                family, srcId = seqId.split('|')[1:3]

                if (family in self.SV_features['FAMILY']) and (srcId in self.SV_features['CYTOBAND']):
                    fasta.seqDict[seqId] = seq

        ##  Skip structure inference if empty database of sequences
        if not fasta.seqDict:
            return {}

        ## 3. Create fasta file
        fastaPath = outDir + '/reference_sequences.fa'
        fasta.write(fastaPath)

        ## 4. Index fasta file
        fileName = 'reference_sequences'  
        indexPath = alignment.index_minimap2(fastaPath, fileName, outDir)

        ## 5. Create fasta file containing consensus inserted sequence
        # Create fasta object
        FASTA = formats.FASTA()
        insert = self.consensusEvent.pick_insert()
        FASTA.seqDict['consensus_insert'] = insert

        # Write fasta
        insertPath = outDir + '/consensus_insert.fa'
        FASTA.write(insertPath)    

        ## 6. Structure inference        
        structure = retrotransposons.retrotransposon_structure(insertPath, indexPath, outDir)

        return structure

    def is_duplication(self, PAF, buffer):
        '''
        Determine if metacluster corresponds to a short tandem duplication
        
        Input:
            1. PAF: PAF object containing consensus inserted sequence alignments on the reference genome
            2. buffer: number of base pairs to extend metacluster interval when assessing overlap

        Output:
            1. DUP: boolean specifying if inserted sequence corresponds to a duplication (True) or not (False)
            2. HITS: PAF object containing inserted sequence alignments supporting a duplication. None if no hit found

        Update SV_features attribute with 'INS_TYPE' and 'PERC_RESOLVED' info

        Note: I need to modify the way PERC_RESOLVED. Now it´s not precise, do it based on the qBeg and qEnd
        alignment coordinates
        '''
        ## 0. Abort if no hit available
        if not PAF.alignments:
            DUP = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

            return DUP, None   

        ## 1. Initialize
        totalPerc = 0
        HITS = formats.PAF()
        
        ## 2. Search hits matching metacluster interval
        # For each hit
        for hit in PAF.alignments: 

            ## Hit in the same ref as the metacluster
            if hit.tName == self.ref: 
                
                overlap, nbBp = gRanges.overlap(self.beg - buffer, self.end + buffer, hit.tBeg, hit.tEnd)
                perc = float(nbBp) / hit.qLen * 100

                ## Hit within metacluster interval
                if overlap:
                    totalPerc += perc
                    HITS.alignments.append(hit)

        ## set upper bound to 100
        if totalPerc > 100:
            totalPerc = 100 

        ## 3. Determine if duplication or not
        # a) Duplication
        if totalPerc >= 40:

            DUP = True
            self.SV_features['INS_TYPE'] = 'duplication'
            self.SV_features['PERC_RESOLVED'] = totalPerc

        # b) Not duplication
        else:

            DUP = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

        return DUP, HITS

    def is_expansion(self, PAF, repeatsDb):
        '''
        Determine if metacluster corresponds to a repeat expansion
        
        Input:
            1. PAF: PAF object containing consensus inserted sequence alignments on the reference genome
            2. repeatsDb: bin database containing annotated repeats in the reference

        Output:
            1. EXPANSION: Boolean specifying if inserted sequence corresponds to an expansion (True) or not (False)
            2. HITS: PAF object containing inserted sequence alignments supporting an expansion. None if no hit found

        Update SV_features attribute with 'INS_TYPE', 'FAMILY', 'SUBFAMILY', 'PERC_RESOLVED' info

        Note: I need to modify the way PERC_RESOLVED. Now it´s not precise, do it based on the qBeg and qEnd
        alignment coordinates        
        '''

        ## 0. Abort if no hit or repeats database not available
        if (not PAF.alignments) or (repeatsDb is None):

            EXPANSION = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

            return EXPANSION, None   
            
        ## 1. Initialize
        totalPerc = 0
        HITS = formats.PAF()

        ## 2. Assess if metacluster located over an annotated repeat sequence 
        # Make list of annotated repeat categories according to repeatmasker
        repeatTypes = ['Low_complexity', 'Simple_repeat', 'Satellite', 'telo', 'centr', 'acro']

        # Selecting those annotated repeats over the target region 
        repeatsFiltered = [ repeat for repeat in self.repeatAnnot if repeat['distance'] == 0 ]

        # Make list of annotated families at the target region
        targetFamilies = [ repeat['family'] for repeat in repeatsFiltered ]
        targetSubfamilies = [ repeat['subfamily'] for repeat in repeatsFiltered ]

        # A) Metacluster over annotated repeat
        if any(family in repeatTypes for family in targetFamilies):

            families = []
            subfamilies = []

            ### Search for inserted sequence hits on repeats of the same family
            # For each hit
            for hit in PAF.alignments: 

                ## Intersect hit coordinates with annotated repeats
                overlaps = annotation.annotate_interval(hit.tName, hit.tBeg, hit.tEnd, repeatsDb)
                
                ## For each repeat overlapping hit coordinates
                for overlap in overlaps:
                    event, nbBp, percBp, coord = overlap
                    
                    # Repeat belonging to the same family as the ones at metacluster interval
                    if event.optional['family'] in targetFamilies:
                        perc = float(nbBp) / hit.qLen * 100
                        totalPerc += perc
                        families.append(event.optional['family'])
                        subfamilies.append(event.optional['subfamily'])

                        if hit not in HITS.alignments:
                            HITS.alignments.append(hit)

            ## set upper bound to 100
            if totalPerc > 100:
                totalPerc = 100 

            ## Determine if expansion or not
            # a) Expansion
            if totalPerc >= 40:

                EXPANSION = True
                self.SV_features['INS_TYPE'] = 'expansion'
                self.SV_features['PERC_RESOLVED'] = totalPerc
                self.SV_features['FAMILY'] = list(set(families))
                self.SV_features['SUBFAMILY'] = list(set(subfamilies))

            # b) Not expansion
            else:

                EXPANSION = False
                self.SV_features['INS_TYPE'] = 'unknown'
                self.SV_features['PERC_RESOLVED'] = 0

        # B) Metacluster outside annotated repeat
        else:
            EXPANSION = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

        return EXPANSION, HITS     

    def is_processed_pseudogene(self, hits, exonsDb):     
        '''
        Determine if metacluster corresponds to a processed pseudogene insertion
        
        Input:
            1. hits: list of bed entries corresponding to inserted sequence hits on the reference 
            2. exonsDb: bin database containing annotated exons in the reference. None if not available

        Output:
            1. PSEUDOGENE: Boolean specifying if inserted sequence corresponds to an expansion (True) or not (False)
            2. outHits: list of bed entries corresponding to hits matching annotated exons
        
        Update SV_features attribute with 'INS_TYPE', 'PERC_RESOLVED', 'NB_EXONS', 'SOURCE_GENE', 'POLYA', 'STRAND'

        Note: I need to modify the way PERC_RESOLVED. Now it´s not precise, do it based on the qBeg and qEnd
        alignment coordinates        
        '''

        ## 0. Abort if no hit available
        if (not hits) or (exonsDb is None):

            PSEUDOGENE = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = 0

            return PSEUDOGENE, None   

        ## 1. Search for polyA/T tails
        insert = self.consensusEvent.pick_insert()
        windowSize = 8
        maxWindowDist = 2
        minMonomerSize = 10
        minPurity = 80 

        # 1.1 Search for poly(A) monomers on the 3' end 
        targetSeq = insert[-50:]
        targetMonomer = 'A'
        monomers3end = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

        ## 1.2 Search for polyT on the 5' end 
        targetSeq = insert[:50]
        targetMonomer = 'T'
        monomers5end = sequences.find_monomers(targetSeq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

        if monomers3end and monomers5end:
            polyA = False
            strand = None

        elif monomers3end:
            polyA = True
            strand = '+'

        elif monomers5end:
            polyA = True
            strand = '-'

        else:
            polyA = False
            strand = None

        ## 2. Intersect each hit with the annotated exons
        outHits = []
        nbResolvedBp = 0
        nbExons = 0

        # For each hit
        for hit in hits: 

            hit.annot = {}

            ## Do intersection
            overlaps = annotation.annotate_interval(hit.ref, hit.beg, hit.end, exonsDb)

            ## Hit intersect with exon
            if overlaps:
                longestOverlap = overlaps[0] # Select exon with longest overlap
                hit.annot['EXON'] = longestOverlap[0] 
                
                outHits.append(hit)
                nbResolvedBp += longestOverlap[1]
                nbExons += 1

        ## 3. Compute percentage of inserted sequence matching on annotated exons
        percResolved = float(nbResolvedBp) / self.consensusEvent.length * 100

        # set upper bound to 100
        if percResolved > 100:
            percResolved = 100 

        ## 4. Determine if pseudogene or not
        # a) Pseudogene
        if percResolved >= 40:

            PSEUDOGENE = True
            self.SV_features['INS_TYPE'] = 'pseudogene'
            self.SV_features['PERC_RESOLVED'] = percResolved
            self.SV_features['NB_EXONS'] = nbExons
            self.SV_features['SOURCE_GENE'] = list(set([hit.annot['EXON'].optional['geneName'] for hit in outHits]))
            self.SV_features['POLYA'] = polyA
            self.SV_features['STRAND'] = strand

        # b) Not pseudogene
        else:
            PSEUDOGENE = False
            self.SV_features['INS_TYPE'] = 'unknown'
            self.SV_features['PERC_RESOLVED'] = percResolved        

        return PSEUDOGENE, outHits
    