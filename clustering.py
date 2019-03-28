'''
Module 'clustering' - Contains functions for clustering sets of objects based on coordinates
'''

## DEPENDENCIES ##
# External
import time
import sys

# Internal
import log
import clusters
import gRanges
import structures

## FUNCTIONS ##
def distance_clustering(binDb, binSize, eventTypes, clusterType, maxDist, minClusterSize):
    '''

    '''
    clustersList = []
    binsInClusters = []

    ## Make list with all the available bins
    binIds = list(binDb.data[binSize].keys())
    binIds.sort() # Sort bins in increasing order

    # For each bin 
    for binIndex in binIds:

        ## Skip bin if already incorporated into a cluster 
        if binIndex in binsInClusters:
            continue

        ## 1. Collect all the events of target event types from current bin ##
        events = binDb.collect_bin(binSize, binIndex, eventTypes)

        ## Skip bin if no target event was found
        if not events:
            continue

        ### 2. Initiate root cluster ##
        cluster = clusters.create_cluster(events, clusterType)
        binsInClusters.append(binIndex) # now bin incorporated into cluster

        ### 3. Root cluster extension
        ## Go forward from current bin (*).  
        #                ---1---> ---2--->
        # |---------|----*----|--------|---------
        # Go one bin forward in each iteration. Extend the cluster 
        # if first event in next bin within 
        # maximum cluster distance or end iteration, otherwise 
        forwardIndex = binIndex + 1

        while True:

            ## Collect events in the next bin 
            events = binDb.collect_bin(binSize, forwardIndex, eventTypes)

            # A) There are events in the bin 
            if events:       

                # Compute the distance between the left most event position in the new bin
                # and the cluster end 
                firstEvent = events[0]
                dist = firstEvent.beg - cluster.end 

                # a) First event within maximum distance 
                if dist <= maxDist:

                    # Add events within bin to the cluster 
                    cluster.add(events)
                    binsInClusters.append(forwardIndex) # now bin incorporated into cluster
                    forwardIndex += 1        

                # b) Event outside -> Stop extension
                else:
                    totalNbEvents = cluster.nbEvents()[0]

                    # Filter out cluster if not composed by enough number of events
                    if totalNbEvents >= minClusterSize:
                        clustersList.append(cluster)

                    break

            # B) No events in the bin -> Stop extension
            else:                
                totalNbEvents = cluster.nbEvents()[0]
                    
                # Filter out cluster if not composed by enough number of events
                if totalNbEvents >= minClusterSize:
                    clustersList.append(cluster)

                break 
    
    return clustersList

def reciprocal_clustering(eventsBinDb, minPercOverlap, minClusterSize, eventType, buffer):
    '''
    '''
    eventsInClusters = []
    clustersDict = {}

    # For each window size/level
    for windowSize in eventsBinDb.binSizes:

        # For each bin in the current window size
        for index in eventsBinDb.data[windowSize]:

            ### 1. Collect all the events in the current bin and 
            # in bins located at higher levels of the hierarchy
            events = eventsBinDb.traverse(index, windowSize, eventType)

            ### 2. Cluster events based on reciprocal overlap
            ## For each event A
            for idx, eventA in enumerate(events):
                
                ## 2.1. Skip comparisons if A already belongs to a cluster 
                if eventA.id in eventsInClusters:
                    continue

                ## 2.2. Generate 2 lists containing clusters and events overlapping A: 
                # - clustersOverlapA: list of clusters overlapping event A
                # - eventsOverlapA: list of events NOT INCLUDED IN A CLUSTER overlapping event A
                clustersOverlapA = [] 
                eventsOverlapA = []

                ## Identify events overlapping A (skip A itself and event pairs already assessed)
                for eventB in events[idx + 1:]:

                    ## Skip comparison if B belongs to a cluster already known to overlap A
                    if (eventB.clusterId in clustersOverlapA):
                        continue
                    
                    ## Add buffer to ranges
                    begA = eventA.beg - buffer
                    endA = eventA.end + buffer
                    begB = eventB.beg - buffer
                    endB = eventB.end + buffer
                    
                    overlap = gRanges.rcplOverlap(begA, endA, begB, endB, minPercOverlap)

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
                    clustersDict[clusterId].add(events2Cluster, None)

                # B) Multiple clusters overlap A -> Merge clusters and add A and its overlapping events into the merged cluster
                elif len(clustersOverlapA) > 1:

                    ## Add events to the list of events already included into clusters
                    events2Cluster = [eventA] + eventsOverlapA 
                    eventsInClusters += [ event.id for event in events2Cluster]

                    ## Make list of clusters overlapping A
                    clusters2merge = [ clustersDict[clusterId] for clusterId in clustersOverlapA ]

                    ## Create merged cluster
                    mergedCluster = clusters.mergeClusters(clusters2merge, eventType)

                    ## Add events to the merged cluster
                    mergedCluster.add(events2Cluster, None)

                    ## Add merged cluster to the clusters dictionary
                    clustersDict[mergedCluster.id] = mergedCluster

                    ## Remove events composing the merged cluster from the clusters dictionary
                    for cluster in clusters2merge:
                        clustersDict.pop(cluster.id, None)
                    
                # C) No cluster overlaps A -> attempt to create a new cluster composed by A and its overlapping events
                else:
                    events2Cluster = [eventA] + eventsOverlapA 
                    clusterSize = len(events2Cluster)

                    # D) Cluster composed by >= X events:
                    if clusterSize >= minClusterSize:
                    
                        # Add events to the list of events already included into clusters
                        eventsInClusters += [ event.id for event in events2Cluster]

                        # Create cluster                        
                        cluster = clusters.create_cluster(events2Cluster, eventType)
                        clustersDict[cluster.id] = cluster

                    # Cluster not composed by enough number of events
    
    clustersList = list(clustersDict.values())

    return clustersList
