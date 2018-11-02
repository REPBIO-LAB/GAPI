'''
Module 'clustering' - Contains functions for clustering sets of objects based on coordinates
'''

## DEPENDENCIES ##
# External
import time

# Internal
import log
import variants
import gRanges

## FUNCTIONS ##
def clusterByDistFast(binHash, maxDist, minClusterSize, clusterType):
    '''
    '''
    clusters = []
    binsInClusters = []

    #print('HOLA: ', binHash)

    # For each bin 
    for binIndex in binHash:

        #print('TIO: ', binIndex)

        # There are elements of the target cluster type in the bin
        if clusterType in binHash[binIndex]:

            #print('MOLON: ', clusterType)

            ### Number of events in the bin
            events = binHash[binIndex][clusterType]
            nbEvents = len(events)

            ### Skip bins fulfilling one of these conditions:
            # a) Bin already incorporated into a cluster 
            if binIndex in binsInClusters:
                continue
    
            # b) Bin without enough number of events to build a root cluster
            elif (nbEvents < minClusterSize):
                continue

            ### Create root cluster containing events on the bin
            cluster = variants.createCluster(events, clusterType)
            clusters.append(cluster)
            binsInClusters.append(binIndex) # now bin incorporated into cluster

            ### Root cluster extension
            ## 2.1 Go backward from current bin (*).  
            #       <---2--- <---1---
            # |---------|--------|----*----|---------
            # Go one bin backward in each iteration. Extend the cluster 
            # if last event in previous bin within 
            # maximum cluster distance or end iteration, otherwise 
            backwardIndex = binIndex - 1

            while True:
        
                # A) There are events in the bin 
                if (backwardIndex in binHash) and (clusterType in binHash[backwardIndex]):
                
                    # Compute the distance between the right most event
                    # in the new bin and the cluster begin 
                    events = binHash[backwardIndex][clusterType]
                    lastEvent = events[-1]
                    dist = cluster.beg - lastEvent.end

                    # a) Last event within maximum distance 
                    if dist <= maxDist:

                        # Add events within bin to the cluster 
                        cluster.add(events, 'left')
                        binsInClusters.append(backwardIndex) # now bin incorporated into cluster
                        backwardIndex -= 1

                    # b) Event outside 
                    else:
                        break

                # B) No events in the bin. Stop iterating
                else:
                    break 

            ## 2.2 Go forward from current bin (*).  
            #                ---1---> ---2--->
            # |---------|----*----|--------|---------
            # Go one bin forward in each iteration. Extend the cluster 
            # if first event in next bin within 
            # maximum cluster distance or end iteration, otherwise 
            forwardIndex = binIndex + 1

            while True:
        
                # A) There are events in the bin 
                if (forwardIndex in binHash) and (clusterType in binHash[forwardIndex]):

                    # Compute the distance between the left most position in the new bin
                    # and the cluster end 
                    events = binHash[forwardIndex][clusterType]
                    firstEvent = events[0]
                    dist = firstEvent.beg - cluster.end 

                    # a) Last event within maximum distance 
                    if dist <= maxDist:

                        # Add events within bin to the cluster 
                        cluster.add(events, 'right')
                        binsInClusters.append(forwardIndex) # now bin incorporated into cluster
                        forwardIndex += 1        

                    # b) Events outside 
                    else:
                        break

                # B) No events in the bin. Stop iterating
                else:
                    break 


def clusterByDist(binDb, maxDist, minClusterSize, clusterType):
    '''
    Cluster events (CLIPPING, INS...) by begin and end coordinates. 

    Two algorithms:
    - maxDist <= bin side -> fast approach
    - maxDist > bin side -> iterative approach

    Input:
        1. binDb: Data structure containing a set of events organized in genomic bins
        2. maxDist: maximum distance between two positions to include them into the same cluster
        3. minClusterSize: minimum number of events required to build a root cluster in a window
        4. clusterType: type of clusters to create (INS: insertion; DEL: deletion; CLIPPING: clipping)

    Output:
    ''' 
    print('INPUT: ', binDb, maxDist, minClusterSize, clusterType, binDb.binSizes)
 
    ## 2. Cluster events ##
    for binSize in sorted(binDb.binSizes, reverse=True):

        binHash = binDb.data[binSize]
        print('TIOO: ', binSize, binHash)


        # a) Fast clustering algorithm
        #if (binSize <= maxDist):
        #    clusterByDistFast(binHash, maxDist, minClusterSize, clusterType)

        # b) Iterative clustering algorithm
        #else:
        #    clusterByDistIterative(binHash, maxDist, minClusterSize, clusterType)
    

    '''
    ## 3. Organize clusters ##
    clustersHash = structures.binsHash(clusters, 'beg', 1000)

    return clustersHash
    '''

def buildMetaClusters(INS_clusters, DEL_clusters, CLIPPING_left_clusters, CLIPPING_right_clusters, targetMetaClusters, confDict):
    '''
    Group clusters into 3 types of metaclusters:

    ** INS: Insertion + clippings **
        beg <-----INS-----> end
                <--CLIP-->
         <--CLIP-->

    ** DEL: Deletion + clippings **
        beg <----DEL----> end
        <--CLIP-->   <--CLIP-->

    ** CLIPPING: Clipping + clipping + ... **
        <--CLIP--> +  <--CLIP--> + ...  
    
    Input:
        1. INS_clusters: insertion clusters organized into a 'binsHash' object. None if not provided
        2. DEL_clusters: deletion '' ''
        3. CLIPPING_left_clusters: left clipping '' ''
        4. CLIPPING_right_clusters: right clipping '' ''
        5. targetMetaClusters: list of types of metaclusters to build (INS: insertion; DEL: deletion; CLIPPING: clipping)
        6. confDict:  
            * maxClusterDist -> maximum cluster distance

    Output:
        1. INS_metaClusters: insertion metaclusters organized into a 'binsHash' object
        2. DEL_metaclusters: deletion metaclusters organized into a '???' object
        3. CLIPPING_metaClusters: clipping metaclusters organized into a 'graph??' object
    ''' 

    if ('INS' in targetMetaClusters):
        build_INS_metaClusters(INS_clusters, CLIPPING_left_clusters, CLIPPING_right_clusters, confDict)


def build_INS_metaClusters(INS_clusters, CLIPPING_left_clusters, CLIPPING_right_clusters, confDict):
    '''
    Group insertion and clipping (left, right) clusters to build insertion metaclusters. 
    
    Note 1: with the current algorithm a single CLIPPING cluster can be assigned to multiple
    INS clusters 
    
    Note 2: with the current algorithm multiple left and right CLIPPING clusters
    can be asigned to the same INS cluster
    '''
    print('INPUT (INS, CLIPPING_LEFT, CLIPPING_RIGHT, maxClusterDist): ', INS_clusters.mapAttribute, INS_clusters.windowSize, INS_clusters.nbEvents(), CLIPPING_left_clusters.nbEvents(), CLIPPING_right_clusters.nbEvents(), confDict['maxClusterDist'])

    ## For each genomic window
    for windowIndex, INS_clusters_window in INS_clusters.data.items():
        
        ## For each insertion cluster in the window
        for INS_cluster in INS_clusters_window:

            ### 1. Search for candidate CLIPPING clusters in nearby windows ###
            ### 1.1 Collect windows within maximum cluster distance
            candidateWindows = [windowIndex]

            ## Go backward from current window (*).  
            #       <---2--- <---1---
            # |---------|--------|----*----|---------
            # Go one window backward in each iteration 
            # while backward window ends within maximum
            # cluster distance or end iteration, otherwise 
            backwardIndex = windowIndex - 1
            
            while True:
            
                windowBeg, windowEnd = INS_clusters.windowBoundaries(backwardIndex)
                dist2beg = INS_cluster.beg - windowEnd
                
                if (dist2beg <= confDict['maxClusterDist']):
                    candidateWindows.append(backwardIndex)
                    backwardIndex -= 1
                else:
                    break

            ## Go forward from current window (*).  
            #                ---1---> ---2--->
            # |---------|----*----|--------|---------
            # Go one window forward in each iteration 
            # while forward window beg within maximum
            # cluster distance or end iteration, otherwise 
            forwardIndex = windowIndex + 1
            
            while True:
            
                windowBeg, windowEnd = INS_clusters.windowBoundaries(forwardIndex)
                dist2end = windowBeg - INS_cluster.end
                
                if (dist2end <= confDict['maxClusterDist']) or (dist2end < 0):
                    candidateWindows.append(forwardIndex)
                    forwardIndex += 1
                else:
                    break   
                
            ### 1.2 Collect candidate CLIPPING clusters from windows
            CLIPPING_left_candidates = []
            CLIPPING_right_candidates = []

            ## For each window
            for windowIndex in candidateWindows:
                
                ## Left clipping clusters available in the window
                if windowIndex in CLIPPING_left_clusters.data:
                    CLIPPING_left_candidates = CLIPPING_left_candidates + CLIPPING_left_clusters.data[windowIndex]

                ## Right clipping clusters available in the window
                if windowIndex in CLIPPING_right_clusters.data:
                    CLIPPING_right_candidates = CLIPPING_right_candidates + CLIPPING_right_clusters.data[windowIndex]
                
            
            ### 2. Create INS metacluster ###
            # Create metacluster if at least one CLIIPING cluster within maximum INS cluster distance

            metaCluster = None

            # Define INS cluster range
            begRange = INS_cluster.beg - confDict['maxClusterDist']
            endRange = INS_cluster.end + confDict['maxClusterDist']

            print('INS-CLUSTER', INS_cluster.beg, INS_cluster.end, begRange, endRange, confDict['maxClusterDist'])

            # For each left CLIPPING cluster assess if within INS cluster range
            for CLIPPING_left_candidate in CLIPPING_left_candidates:

                overlap = gRanges.overlap(CLIPPING_left_candidate.beg, CLIPPING_left_candidate.end, begRange, endRange)
                
                if overlap:
                    print('INTERSECTION-LEFT: ', CLIPPING_left_candidate.beg, CLIPPING_left_candidate.end, gRanges.overlap(CLIPPING_left_candidate.beg, CLIPPING_left_candidate.end, begRange, endRange))
                    
                    ## TILL HERE

                    # a) Create metacluster

                    # b) Add to an existing metacluster 

            # For each right CLIPPING cluster assess if within INS cluster range
            for CLIPPING_right_candidate in CLIPPING_right_candidates:
                overlap = gRanges.overlap(CLIPPING_right_candidate.beg, CLIPPING_right_candidate.end, begRange, endRange)

                if overlap:
                    print('INTERSECTION-RIGHT: ', CLIPPING_right_candidate.beg, CLIPPING_right_candidate.end, gRanges.overlap(CLIPPING_right_candidate.beg, CLIPPING_right_candidate.end, begRange, endRange))

                    ## TILL HERE
                   
                    # a) Create metacluster

                    # b) Add to an existing metacluster 


        print('-------------------')

