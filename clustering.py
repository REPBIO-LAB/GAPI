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

    # For each window size/level
    for windowSize in binDb.binSizes:

        # For each bin in the current window size
        for index in binDb.data[windowSize]:

            ### 1. Collect all the events in the current bin and 
            # in bins located at higher levels of the hierarchy
            events = binDb.traverse(index, windowSize, eventTypes)

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
                    
                    overlap, overlapLen = gRanges.rcplOverlap(begA, endA, begB, endB, minPercOverlap)

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

                    ## Remove clusters that were merged from the clusters dictionary (Investigate later why if enabled the codes break)
                    for cluster in clusters2merge:
                        clustersDict.pop(cluster.id, None)

                    ## Update clustersOverlapA
                    if 'META' in mergedCluster.id:
                        if cluster.id in clustersOverlapA:
                            clustersOverlapA.remove(cluster.id)
                            if mergedCluster.id not in clustersOverlapA:
                                clustersOverlapA.append(mergedCluster.id)

                # C) No cluster overlaps A
                else:

                    events2Cluster = [eventA] + eventsOverlapA 

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



def reciprocal(binDb, minPercOverlap, minClusterSize, buffer):
    '''
    '''
    eventsMinusDict = {} 
    eventsPlusDict = {}
    reciprocalDict = {}

    # For each window size/level
    for windowSize in binDb.binSizes:

        # For each bin in the current window size
        for index in binDb.data[windowSize]:

            # Collect all event types present in the binDb
            commonEventTypes = set([i.split('-', 1)[1] for i in binDb.eventTypes])

            for commonEventType in commonEventTypes:
                # cogo solo las partes comunes de los eventTypes
                #commonEventType =  '-'.join(eventType.split('-')[1:])
                # Hago una lista con los clusters de ese tipo que sean PLUS
                # Hago una lista con los clusters de ese tipo que sean MINUS
                plusEventType =  'PLUS-' + commonEventType
                minusEventType =  'MINUS-' + commonEventType

                #for actualEventType in binDb.data[windowSize][index]:
                # cogo solo las partes comunes de los eventTypes
                #actualCommonEventType =  '-'.join(actualEventType.split('-')[1:])
                #print (actualCommonEventType)
                #plusEvents = []
                #minusEvents = []

                #if actualCommonEventType == commonEventType:
                plusEvents = binDb.traverse(index, windowSize, [plusEventType])
                # Append events from the adjacent left bin
                plusEvents.extend(binDb.collect_bin(windowSize, index-1, plusEventType))
                # Append events from the adjacent right bin
                plusEvents.extend(binDb.collect_bin(windowSize, index+1, plusEventType))

                #if actualCommonEventType == commonEventType:
                minusEvents = binDb.traverse(index, windowSize, [minusEventType])
                # Append events from the adjacent left bin
                minusEvents.extend(binDb.collect_bin(windowSize, index-1, minusEventType))
                # Append events from the adjacent right bin
                minusEvents.extend(binDb.collect_bin(windowSize, index+1, minusEventType))

                # si ninguna de las dos listas esta vacia:
                if len(plusEvents) > 0 and len(minusEvents) > 0:
                    ### 2. Cluster events based on reciprocal overlap
                    ## For each event A
                    for eventPlus in plusEvents:
                        eventsPlusDict[eventPlus] = {}
                        
                        ## [SR CHANGE]: Look at the previous event in order to not split a cluster
                        ## Identify events overlapping A (skip A itself and event pairs already assessed)
                        for eventMinus in minusEvents:

                            overlapDict = {}

                            ## Skip comparison if B belongs to a cluster already known to overlap A
                            #if (eventB.clusterId in clustersOverlapA):
                                #continue
                            
                            ## Add buffer to ranges
                            begA = eventPlus.beg - buffer
                            endA = eventPlus.end + buffer
                            begB = eventMinus.beg - buffer
                            endB = eventMinus.end + buffer
                            
                            overlap, overlapLen = gRanges.rcplOverlap(begA, endA, begB, endB, minPercOverlap)

                            # A) Event B overlap A. 
                            if overlap:
                                overlapDict[eventMinus] = overlapLen
                        # Una vez que tengo mirado, para un eventPlus, los eventMinus que lo solapan (lo he comparado con todos), cojo el que tenga mayor overlapLen (pq significa que estan mas cerca)
                        # Cojo la key del overlap dict cuyo value sea mayor (pq significa que overlapan mas, sera la pareja)
                        eventMinusReciprocal, eventMinusReciprocalOverlapLen = max(overlapDict.items(), key=operator.itemgetter(1))

                        eventsPlusDict[eventPlus][eventMinusReciprocal]=eventMinusReciprocalOverlapLen

                    # VER EXPLANATIOOOOOOOOOOOOON
                    # Una vez que tengo mirados todos los eventPlus, cojo el mejor eventPlus para cada eventMinus.
                    # Para eso, le doy la vuelta al dictionary
                    # TAMBIEN HAY SOLUCION DE UNA SOLA LINEA!!!!! BUSCARLA EN MARCADORES!

                    for eventPlus, nestDict in eventsPlusDict.items():
                        for eventMinus, overlapLen in nestDict.items():
                            eventsMinusDict.setdefault(eventMinus, {}).update({eventPlus:overlapLen})
                    
                    # Y me quedo solo con aquellos eventsPlus que tengan la mayor overlapLen, asi ya tengo los pares
                    for eventMinus,nestDict in eventsMinusDict.items():
                        #reciprocalDict[eventMinus] = {}
                        reciprocalDict[commonEventType] = []
                        #reciprocalDict[commonEventType].append(eventMinus)
                        reciprocalDict[commonEventType].extend(eventMinus.events)
                        if len(nestDict) > 1:
                            eventPlus_to_add, overlapLen_to_add = max(nestDict.items(), key=operator.itemgetter(1))
                            #reciprocalDict[eventMinus][eventPlus_to_add] = overlapLen_to_add
                            #reciprocalDict[commonEventType].append(eventPlus_to_add)
                            reciprocalDict[commonEventType].extend(eventPlus_to_add.events)
                        else:
                            #reciprocalDict[eventMinus] = nestDict
                            #reciprocalDict[commonEventType].append(*nestDict)
                            eventPlus = *nestDict,[0]
                            reciprocalDict[commonEventType].extend(eventPlus[0].events)

                elif len(plusEvents) > 0:

                    for eventPlus in plusEvents:
                        reciprocalDict[commonEventType] = []
                        reciprocalDict[commonEventType].extend(eventPlus.events)

                elif len(minusEvents) > 0:

                    for eventMinus in minusEvents:
                        reciprocalDict[commonEventType] = []
                        reciprocalDict[commonEventType].extend(eventMinus.events)

    # Example of reciprocalDict
    # reciprocalDict = {'DISCORDANT-Hepatitis': [<clusters.DISCORDANT_cluster object at 0x7f32a07230b8>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>], 'DISCORDANT-UNVERIFIED:': [<clusters.DISCORDANT_cluster object at 0x7f32a07230b8>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>], 'DISCORDANT-HBV': [<clusters.DISCORDANT_cluster object at 0x7f32a07230f0>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>]}
    # Aqui podria retornar una lista de clusters, como aparece en el ejemplo de arriba, pero para hacer el metaclustering necesito events, asi que voy a retornar los events directamente
    '''
    for tipo, events in reciprocalDict.items():
        print ('tipo ' + str(tipo))
        for event in events:
            print ('side ' + str(event.side))
            print ('type ' + str(event.type))
            print ('identity ' + str(event.identity))
    '''
    return reciprocalDict

    ## ----------------------------------------- EXPLANATIOOOOOOOOOOOOON -----------------------------------------
    # Hasta aqui en eventsPlusDict se ha hecho lo siguiente: para cada eventPlus se ha cogido el mejor hit de eventMinus, pq para cada eventPlus se han mirado TODOS los eventMinus y nos hemos quedado con el mejor.
    # Dictionario de ejemplo:
    '''
    eventsPlusDict = {'clusterPlus1': {'clusterMinus1': 576}, 'clusterPlus2': {'clusterMinus1': 556}, 'clusterPlus3': {'clusterMinus3': 459}, 'clusterPlus4': {'clusterMinus3': 20}, 'clusterPlus5': {'clusterMinus5': 20}} 
    '''
    # Sin embargo, para cada eventMinus tambien se han mirado todos los eventPlus, pero por separado, de manera que un mismo eventMinus se pudo haber asignado a dos eventPlus, siendo uno mejor que el otro (como en el caso de clusterMinus1 y clusterMinus3).
    # Para chequear si esto ultimo pasa, le damos la vuelta al dictionary:
    '''
    result = {} 
    for k, v in pairsDict.items():
        for a,b in v.items():
            result.setdefault(a, {}).update({k:b})
    print (result)
    # {'clusterMinus1': {'clusterPlus1': 576, 'clusterPlus2': 556}, 'clusterMinus3': {'clusterPlus3': 459, 'clusterPlus4': 20}, 'clusterMinus5': {'clusterPlus5': 20}}
    '''
    # Ahora si que podemos escoger, para cada cluster minus, cual es el mejor clusterPlus, y ya hacemos el dictionary final, que contendra solo los pares reciprocos, junto con la overlapLen
    '''
    final = {}

    for k,v in result.items():
        final[k] = {}
        if len(v) > 1:
            key_to_add, value_to_add = max(v.items(), key=operator.itemgetter(1))
            print (key_to_add)
            print (value_to_add)
            final[k][key_to_add] = value_to_add
        else:
            final[k] = v

    print (final)
    # {'clusterMinus1': {'clusterPlus1': 576}, 'clusterMinus3': {'clusterPlus3': 459}, 'clusterMinus5': {'clusterPlus5': 20}}
    '''
    ## ----------------------------------------- EXPLANATIOOOOOOOOOOOOON -----------------------------------------
