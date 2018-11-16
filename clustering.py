'''
Module 'clustering' - Contains functions for clustering sets of objects based on coordinates
'''

## DEPENDENCIES ##
# External
import time
import itertools

# Internal
import log
import variants
import gRanges
import structures

## FUNCTIONS ##
def clusterByDist1D(binDb, maxDist, minClusterSize, eventType):
    '''
    Cluster events (CLIPPING, INS...) by begin positions. 

    Two algorithms:
    - maxDist <= bin side -> fast approach
    - maxDist > bin side -> iterative approach

    Input:
        1. binDb: Data structure containing a set of events organized in genomic bins
        2. maxDist: maximum distance between two positions to include them into the same cluster
        3. minClusterSize: minimum number of events required to build a root cluster in a window
        4. eventType: type of events to cluster (INS: insertion; DEL: deletion; CLIPPING: clipping)

    Output:
    ''' 
    binSize = binDb.binSizes[0]
    
    ## 1. Create clusters ##
    # a) Fast clustering algorithm
    if (binSize <= maxDist):
        clusters = clusterByDist1D_fast(binDb.data[binSize], maxDist, minClusterSize, eventType)

    # b) Iterative clustering algorithm (not implemented yet)
    else:
        print('ITERATIVE ALGORITHM. NOT IMPLEMENTED YET...')
        # clusterByDist1D_iterative(binDb.data[binSize], maxDist, minClusterSize, eventType)
    
    ## 2. Organize clusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    data = [(clusters, eventType + '-CLUSTER')]
    clustersBins = structures.createBinDb(data, binSizes)

    return clustersBins


def clusterByDist1D_fast(binHash, maxDist, minClusterSize, eventType):
    '''
    '''
    clusters = []
    binsInClusters = []

    # For each bin 
    for binIndex in binHash:

        # There are elements of the target cluster type in the bin
        if eventType in binHash[binIndex]:

            ### Number of events in the bin
            binObj = binHash[binIndex][eventType]

            ### Skip bins fulfilling one of these conditions:
            # a) Bin already incorporated into a cluster 
            if binIndex in binsInClusters:
                continue
    
            # b) Bin without enough number of events to build a root cluster
            elif (binObj.nbEvents() < minClusterSize):
                continue

            ### 1. Create root cluster containing events on the bin
            cluster = variants.createCluster(binObj.events, eventType)
            clusters.append(cluster)
            binsInClusters.append(binIndex) # now bin incorporated into cluster

            ### 2. Root cluster extension
            ## 2.1 Go backward from current bin (*).  
            #       <---2--- <---1---
            # |---------|--------|----*----|---------
            # Go one bin backward in each iteration. Extend the cluster 
            # if last event in previous bin within 
            # maximum cluster distance or end iteration, otherwise 
            backwardIndex = binIndex - 1

            while True:
        
                # A) There are events in the bin 
                if (backwardIndex in binHash) and (eventType in binHash[backwardIndex]):
                
                    # Compute the distance between the right most event
                    # in the new bin and the cluster begin 
                    newBinObj = binHash[backwardIndex][eventType]
                    lastEvent = newBinObj.events[-1]
                    dist = cluster.beg - lastEvent.end

                    # a) Last event within maximum distance 
                    if dist <= maxDist:

                        # Add events within bin to the cluster 
                        cluster.add(newBinObj.events, 'left')
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
                if (forwardIndex in binHash) and (eventType in binHash[forwardIndex]):

                    # Compute the distance between the left most position in the new bin
                    # and the cluster end 
                    newBinObj = binHash[forwardIndex][eventType]
                    firstEvent = newBinObj.events[0]
                    dist = firstEvent.beg - cluster.end 

                    # a) Last event within maximum distance 
                    if dist <= maxDist:

                        # Add events within bin to the cluster 
                        cluster.add(newBinObj.events, 'right')
                        binsInClusters.append(forwardIndex) # now bin incorporated into the cluster
                        forwardIndex += 1        

                    # b) Events outside 
                    else:
                        break

                # B) No events in the bin. Stop iterating
                else:
                    break 

    return clusters


def clusterByRcplOverlap(binDb, minPercOverlap, minClusterSize, eventType):
    '''
    '''

    eventsInClusters = []
    clusters = {}

    # For each window size/level
    for windowSize in binDb.binSizes:

        # For each bin in the current window size
        for index in binDb.data[windowSize]:

            ### 1. Collect all the events in the current bin and 
            # in bins located at higher levels of the hierarchy
            events = binDb.traverse(index, windowSize, eventType)

            ### 2. Cluster events based on reciprocal overlap
            ## For each event A
            for idx, eventA in enumerate(events):
                
                ## 2.1. Skip comparisons if A already belongs to a cluster 
                if eventA.id in eventsInClusters:
                    continue

                ## 2.2. Generate 2 lists of clusters and events overlapping A: 
                # - clustersOverlapA: list of clusters overlapping event A
                # - eventsOverlapA: list of events NOT INCLUDED IN A CLUSTER overlapping event A
                clustersOverlapA = [] 
                eventsOverlapA = []

                ## Identify events overlapping A (skip A itself and event pairs already assessed)
                for eventB in events[idx + 1:]:

                    ## Skip comparison if B belongs to a cluster already known to overlap A
                    if (eventB.clusterId in clustersOverlapA):
                        continue
                    
                    overlap = gRanges.rcplOverlap(eventA.beg, eventA.end, eventB.beg, eventB.end, minPercOverlap)

                    # A) Event B overlap A. 
                    if overlap:

                        # a) B already belongs to a cluster. So this cluster overlaps A 
                        if eventB.clusterId != None: 
                            clustersOverlapA.append(eventB.clusterId)

                        # b) B does not belong to any cluster
                        else:
                            eventsOverlapA.append(eventB)

                    # B) Event B NOT overlap A                        
                
                ## 2.3. Finish by adding A and its overlapping events to a cluster or creating a cluster 
                # A) One cluster overlaps A -> Add A and its overlapping events into the cluster
                if len(clustersOverlapA) == 1:

                    # Add events to the list of events already included into clusters
                    events2Cluster = [eventA] + eventsOverlapA 
                    eventsInClusters += [ event.id for event in events2Cluster]

                    # Add events to the cluster
                    clusterId = clustersOverlapA[0]
                    clusters[clusterId].add(events2Cluster, None)

                # B) Multiple clusters overlap A -> Merge clusters and add A and its overlapping events into the merged cluster
                elif len(clustersOverlapA) > 1:

                    ## Add events to the list of events already included into clusters
                    events2Cluster = [eventA] + eventsOverlapA 
                    eventsInClusters += [ event.id for event in events2Cluster]

                    ## Make list of clusters overlapping A
                    clusters2merge = [ clusters[clusterId] for clusterId in clustersOverlapA ]

                    ## Create merged cluster
                    mergedCluster = variants.mergeClusters(clusters2merge, eventType)

                    ## Add events to the merged cluster
                    mergedCluster.add(events2Cluster, None)

                    ## Add merged cluster to the clusters dictionary
                    clusters[mergedCluster.id] = mergedCluster

                    ## Remove events composing the merged cluster from the clusters dictionary
                    for cluster in clusters2merge:
                        clusters.pop(cluster.id, None)
                    
                # C) No cluster overlaps A -> attempt to create a new cluster composed by A and its overlapping events
                else:
                    events2Cluster = [eventA] + eventsOverlapA 
                    clusterSize = len(events2Cluster)

                    # D) Cluster composed by >= X events:
                    if clusterSize >= minClusterSize:
                    
                        # Add events to the list of events already included into clusters
                        eventsInClusters += [ event.id for event in events2Cluster]

                        # Create cluster                        
                        cluster = variants.createCluster(events2Cluster, eventType)
                        clusters[cluster.id] = cluster

                    # Cluster not composed by enough number of events
    
    ## 2. Organize clusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]

    data = [(list(clusters.values()), eventType + '-CLUSTER')]
    clustersBins = structures.createBinDb(data, binSizes)

    return clustersBins

