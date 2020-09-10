'''
Module 'clustering' - Contains functions for clustering sets of objects based on different criteria
'''

## DEPENDENCIES ##
# External
import time
import sys
import operator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np

# Internal
import log
import clusters
import gRanges
import structures


## FUNCTIONS ##

def distance_clustering_targetPos(events, maxDist, pos2cluster):
    '''
    Group events based on their begin or end distances into clusters
    
    Input:
        1. events: list of events to be clustered
        2. maxDist: maximum distance between two events to be clustered together
        3. pos2cluster: 'beg' or 'end' to cluster events based on their begin or end positions, respectively

    Output:
        1. clusterList: nested list containing events clustered together
    '''
    ## 1. Sort events in increasing beg or end positions
    operatorObj = operator.attrgetter(pos2cluster) 
    events.sort(key=operatorObj) 

    ## 2. Make clustering
    clusterList = []

    # For each event
    for event in events: 

        # A) No cluster available -> Initialize first cluster
        if not clusterList:
            currentCluster = clusters.SUPPLEMENTARY_cluster([event])
            currentCluster.bkpSide = pos2cluster
            clusterList.append(currentCluster)

        # B) Clusters available
        else:
    
            # Compare event coordinates with the coordinates from the last event in cluster 
            pos = getattr(event, pos2cluster)
            
            lastEvent = currentCluster.events[-1]
            lastPos = getattr(lastEvent, pos2cluster)

            # a) Add event to cluster as within max distance
            if (pos - lastPos <= maxDist):
                currentCluster.add([event])
        
            # b) Create new cluster
            else:
                currentCluster = clusters.SUPPLEMENTARY_cluster([event])
                currentCluster.bkpSide = pos2cluster                
                clusterList.append(currentCluster)

    return clusterList

def distance_clustering(binDb, binSize, eventTypes, clusterType, maxDist, minClusterSize):
    '''
    Group events located at a given bin size level based on position distance into clusters
    
    Input:
        1. binDb: data structure containing a set of events organized in genomic bins  
        2. binSize: bin size level to do the clustering
        3. eventTypes: list with target event types to be clustered together
        4. clusterType: type of clusters to be created
        5. maxDist: maximum distance between two events to be clustered together
        6. minClusterSize: minimum number of events clustering together for creating a cluster

    Output:
        1. clustersList: list of created clusters
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

def reciprocal_overlap_clustering(binDb, minPercOverlap, minClusterSize, eventTypes, buffer, clusterType):
    '''
    Group events/clusters based on reciprocal overlap into clusters/metaclusters

    Input:
        1. binDb: data structure containing a set of events/clusters organized in genomic bins  
        2. minPercOverlap: minimum percentage of reciprocal overlap to cluster two events together
        3. minClusterSize: minimum number of events clustering together for creating a cluster
        4. eventTypes: list with target event types to be clustered together
        5. buffer: number of nucleotides used to extend cluster begin and end coordinates prior evaluating reciprocal overlap 
        6. clusterType: type of clusters to be created (If "META", metaclustering will be performed)

    Output:
        1. clustersList: list of created clusters/metaclusters
    '''    
    eventsInClusters = []
    clustersDict = {}
    
    # For each window size
    for windowSize in binDb.binSizes:
        
        # For each bin in the current window size
        # It needs to be sorted. Index_bin = beg / windowSize
        for index in sorted(binDb.data[windowSize]):
            
            ### 1. Collect all the events in the current bin and its contiguous
            events = []
            for window in [index, index + 1]:
                events += binDb.traverse(window, windowSize, eventTypes)
            
            events.sort(key=lambda event: event.beg)
            
            ### 2. Cluster events based on reciprocal overlap
            ## For each event A
            for idx, eventA in enumerate(events):
                
                print('idx, eventA')
                print(idx, eventA)
                
                ## 2.1. Skip comparisons if A already belongs to a cluster 
                if eventA.id in eventsInClusters:
                    print('if eventA.id in eventsInClusters')
                    continue

                ## 2.2. Generate 2 lists containing clusters and events overlapping A: 
                # - clustersOverlapA: list of clusters overlapping event A
                # - eventsOverlapA: list of events NOT INCLUDED IN A CLUSTER overlapping event A
                clustersOverlapA = [] 
                eventsOverlapA = []

                ## Identify events overlapping A (skip A itself and event pairs already assessed)
                for eventB in events[idx + 1:]:
                    
                    print('eventB')
                    print(eventB)
                    
                    ## Skip comparison if B belongs to a cluster already known to overlap A
                    if (eventB.clusterId in clustersOverlapA):
                        print('if (eventB.clusterId in clustersOverlapA)')
                        continue
                    
                    ## Add buffer to ranges
                    begA = eventA.beg - buffer
                    endA = eventA.end + buffer
                    begB = eventB.beg - buffer
                    endB = eventB.end + buffer
                    
                    overlap, overlapLen = gRanges.rcplOverlap(begA, endA, begB, endB, minPercOverlap)
                    
                    print('overlap, overlapLen')
                    print(overlap, overlapLen)
                                
                    # A) Event B overlap A. 
                    if overlap:
                        print('if overlap:')
                        
                        # a) B already belongs to a cluster. So this cluster overlaps A 
                        if eventB.clusterId != None:
                            print('if eventB.clusterId != None:')
                            print(eventB)
                            print(eventB.clusterId)
                            
                            clustersOverlapA.append(eventB.clusterId)
                            
                            print(clustersOverlapA)

                        # b) B does not belong to any cluster
                        else:
                            print('eventsOverlapA.append(eventB)')
                            
                            eventsOverlapA.append(eventB)
                            
                            print(eventsOverlapA)

                    # B) Event B NOT overlap A                        
                
                print('# B) Event B NOT overlap A  ')
                print(clustersOverlapA)
                print(eventsOverlapA)
                
                ## 2.3. Finish by adding A and its overlapping events to a cluster or creating a cluster 
                # A) One cluster overlaps A -> Add A and its overlapping events into the cluster
                if len(clustersOverlapA) == 1:

                    # Add events to the list of events already included into clusters
                    events2Cluster = [eventA] + eventsOverlapA 
                    eventsInClusters += [event.id for event in events2Cluster]

                    # Add events to the cluster
                    clusterId = clustersOverlapA[0]
                    clustersDict[clusterId].add(events2Cluster)

                # B) Multiple clusters overlap A -> Merge clusters and add A and its overlapping events into the merged cluster
                elif len(clustersOverlapA) > 1:

                    ## Make list of clusters overlapping A
                    clusters2merge = [ clustersDict[clusterId] for clusterId in clustersOverlapA ]

                    ## Create merged cluster                    
                    mergedCluster = clusters.merge_clusters(clusters2merge, clusterType)

                    ## Add events to the list of events already included into clusters
                    events2Cluster = [eventA] + eventsOverlapA 
                    eventsInClusters += [ event.id for event in events2Cluster]

                    ## Add events to the merged cluster
                    mergedCluster.add(events2Cluster)

                    ## Add merged cluster to the clusters dictionary
                    clustersDict[mergedCluster.id] = mergedCluster

                    ## Remove clusters that were merged from the clusters dictionary 
                    #for cluster in clusters2merge:
                        #clustersDict.pop(cluster.id, None)

                # C) No cluster overlaps A
                else:

                    events2Cluster = [eventA] + eventsOverlapA 

                    ## CAUTION: Can be problematic when using it on clusters!!!
                    # D) A + overlapping events would make a cluster composed by >= minClusterSize:
                    if len(events2Cluster) >= minClusterSize:

                        # Add events to the list of events already included into clusters
                        eventsInClusters += [ event.id for event in events2Cluster]

                        # Create cluster                        
                        cluster = clusters.create_cluster(events2Cluster, clusterType)

                        # Add cluster to the dict
                        clustersDict[cluster.id] = cluster

                    # Cluster not composed by enough number of events
    
    clustersList = list(clustersDict.values())
    
    return clustersList


def KMeans_clustering(events, x, y):
    '''
    Group events into clusters through K-means and one or two attributes 

    Input:
        1. events: List of events to be clustered
        2. x: Attribute used for clustering. If None, the X-axis will not be taken into account for clustering
        3. y: Attribute used for clustering. If None, the Y-axis will not be taken into account for clustering
        4. offset_x: Substract offset to x attribute value (TO DO)
        5. offset_y: Substract offset to y attribute value (TO DO)

    Output:
        1. groups: Dictionary containing grouped events according to K-means clusters (keys)

    NOTE: To transform genomic coordinates into cluster interval coordinates I can use an offset variable for x and y
    KMeans_clustering(events, x, y, x_offset, y_offset) * then I would make the division x-x_offset and/or y-y_offset
    '''
    ## 1. Exit if not enough number of events for clustering
    nbEvents = len(events)

    if nbEvents < 3:
        return {}

    ## 2. Generate nested list with event´s attribute values that will be used for clustering
    data = []

    for event in events:

        ## Define X
        if x is None:
            x_value = 1
        else:
            x_value = getattr(event, x)

        ## Define Y
        if y is None:
            y_value = 1
        else:
            y_value = getattr(event, y)

        ## Add to list
        data.append([x_value, y_value])

    ## 3. Define K values to be tried
    Ks = [k for k in range(2, nbEvents)]

    ## 4. Perform clustering with K-means
    max_coefficient, max_labels = KMeans_multiK(data, Ks)

    if max_coefficient < 0.60:
        return {}

    ## 5. Group events according to K-means clusters
    groups = {}

    for index, label in enumerate(max_labels):

        event = events[index]

        if label not in groups:
            groups[label] = [event]
            
        else:
            groups[label].append(event)

    return groups


def KMeans_multiK(data, Ks):
    '''
    Perform K-means clustering with multiple K-values and return clustering results for the K that maximizes the average silhouette coefficient

    Input:
        1. data: Nested list composed by a single list containing a two elements list per sample with features used for clustering
        2. Ks: List of K values to be used
        
    Output:
        1. max_coefficient: Maximum average Silhouette coefficient obtained with an input K value
        2. max_labels: List containing cluster labels (0, 1, ...) for the samples provided in the input 'data' variable
    '''

    ## Apply K-means clustering for each input K value 
    max_coefficient = 0
    max_labels = []

    for k in Ks:

        ## 1. Do K-means
        kmeans = KMeans(n_clusters=k, random_state=10)
        labels = kmeans.fit_predict(data) 
                    
        ## 2. Compute average silhouette coefficient
        coefficient = silhouette_score(data, labels)

        # a) Coefficient increases with current K value 
        if coefficient > max_coefficient:
            max_coefficient = coefficient
            max_labels = labels

        # b) Stop iterating as silhouette coefficient does not increase
        else:
            break

    return max_coefficient, max_labels

def distance_clustering_SR(events, minPercOverlap, minClusterSize, eventTypes, clusterType, equalOrientBuffer, oppositeOrientBuffer, libraryReadLength):
    '''
    Sequentially group sorted events based on distance clustering into metaclusters

    Input:
        1. events: list of sorted events. They must be sorted in the way they are going to be clustered.
        2. minPercOverlap: minimum percentage of reciprocal overlap to cluster two events together
        3. minClusterSize: minimum number of events clustering together for creating a cluster
        4. eventTypes: list with target event types to be clustered together
        5. clusterType: type of clusters to be created (If "META", metaclustering will be performed)
        6. equalOrientBuffer: Distance between reads that are equally oriented.
        7. oppositeOrientBuffer: Distance between reads that are opposite oriented.
        8. libraryReadLength: Illumina library read length.

    Output:
        1. metaclustersList: list of created clusters/metaclusters
    '''    

    # Inicialize variables and set thresholds
    metaclustersList = []
    equalOrientThreshold = equalOrientBuffer + libraryReadLength
    oppOrientThreshold = oppositeOrientBuffer + libraryReadLength

    # Make numpy array with events beginnings
    eventsStarts=np.array([event.beg for event in events])
    #print ('eventsStarts ' + str(eventsStarts))

    #print ([event.readName for event in events])

    # Calculate differences between event begginings
    differences = np.diff(eventsStarts)
    #print ('differences ' + str(differences))

    # Choose these differences which are greater than clustering threshold and their indexes.
    indices = [(i,differences[i]) for i,v in enumerate(differences >= equalOrientThreshold) if v]
    #print ('indices ' + str(indices))

    # If there are greater differences:
    if indices:
        start = 0
        # Go through differences
        for index, diff in indices:
            #print ('1 ' + str(index) + ' '+ str(diff) + ' '+ str(start))
            # If difference lower than second threhold and events are opposite or is greather than second clustering threshold
            if (diff <= oppOrientThreshold and events[index].orientation == events[index+1].orientation) or (diff > oppOrientThreshold):
                    #print ('2 ' + str(index) + ' '+ str(diff) + ' '+ str(start) +' '+ str(events[index].orientation) +' '+ str(events[index+1].orientation) +' '+ str(events[index].readName) +' '+ str(events[index+1].readName))
                    # And the list has more than one element, cluster events
                    if index-start > 0:
                        #print ('3 ' + str(index) + ' '+ str(diff) + ' '+ str(start) +' '+ str(events[index].orientation) +' '+ str(events[index+1].orientation) +' '+ str(events[index].readName) +' '+ str(events[index+1].readName))
                        events2cluster = events[start:index+1] # mirar bien este rango
                        cluster = clusters.create_cluster(events2cluster, clusterType)
                        metaclustersList.append(clusters.create_cluster([cluster], 'META'))
                    # Re-new start index
                    start = index + 1
                    #print ('4 ' + str(index) + ' '+ str(diff) + ' '+ str(start) +' '+ str(events[index].orientation) +' '+ str(events[index+1].orientation) +' '+ str(events[index].readName) +' '+ str(events[index+1].readName))

        #print ('index ' + str(index))
        #print ('len(events) ' + str(len(eventsStarts)))
        # Cluster events from last difference to last event element:
        finalDiff = len(events) - start
        #print ('finalDiff ' + str(finalDiff))
        if finalDiff > 1: # si la lista es larga
            #make cluster
            #print ('5 ' + str(index+1) + ' '+ str(len(events)))
            events2cluster = events[start:len(events)] # mirar bien este rango
            cluster = clusters.create_cluster(events2cluster, clusterType)
            metaclustersList.append(clusters.create_cluster([cluster], 'META'))
            '''
            for even in cluster.events:
                print ('even.beg ' + str(even.beg) +' '+ str(even.orientation) +' '+ str(even.readName))
            '''
            #print ('2.clusterDiff ' + str(cluster.end - cluster.beg) +' '+ str(cluster.ref) +' '+ str(cluster.beg) +' '+ str(cluster.end))
    
    # If there are no greater differences and the input events list has more elements than 1, cluster them:
    elif len(events) > 1:
        #for even in events:
            #print ('No indice ' + str(even.readName))
        cluster = clusters.create_cluster(events, clusterType)
        '''
        for even in cluster.events:
            print ('even.beg ' + str(even.beg) +' '+ str(even.orientation) +' '+ str(even.readName))
        '''
        #print ('3.clusterDiff ' + str(cluster.end - cluster.beg) +' '+ str(cluster.ref) +' '+ str(cluster.beg) +' '+ str(cluster.end))

        metaclustersList.append(clusters.create_cluster([cluster], 'META'))
        # Add cluster to the dict
        #clustersDict[cluster.id] = cluster    
    
    #clustersList = list(clustersDict.values())

    #for meta in metaclustersList:
        #print ('META events ' + str(meta.ref) + ' ' + str(meta.beg) + ' ' + str(len(meta.events)) +' '+ str(meta.events[0].beg) + ' ' + str(meta.events[-1].beg) +' '+ str(meta.events[-1].beg - meta.events[0].beg))
    return metaclustersList


