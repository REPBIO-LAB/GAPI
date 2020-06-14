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

import os

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



def reciprocal_DEPRECATED(binDb, minPercOverlap, minClusterSize, buffer):
    '''
    Make reciprocal clusters from PLUS and MINUS discordant clusters. 
    Those clusters that have opposite orientation, are close (according to buffer parameter) and have the same identity are clustered together.

    Input:
        1. binDb: data structure containing a set of events organized in genomic bins 
        2. minPercOverlap: minimum percentage of reciprocal overlap to cluster two events together
        3. minClusterSize: minimum number of events clustering together for creating a cluster
        4. buffer: number of nucleotides used to extend cluster begin and end coordinates prior evaluating reciprocal overlap 
        
    Output:
        1. reciprocalDict: dictionary containing list of clusters grouped according to their type. The type is defined by ORIENTATION-EVENTSTYPE-IDENTITY (i.e. PLUS-DISCORDANT-HBV, RECIPROCAL-DISCORDANT-HBV, MINUS-DISCORDANT-HBV)
    '''

    eventsMinusDict = {} 
    eventsPlusDict = {}
    reciprocalDict = {}

    # Collect all event types present in the binDb
    commonEventTypes = set([i.split('-', 1)[1] for i in binDb.eventTypes])

    # For each window size/level
    for windowSize in binDb.binSizes:

        # For each bin in the current window size
        for index in binDb.data[windowSize]:

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
                ##print (actualCommonEventType)
                #plusEvents = []
                #minusEvents = []

                #if actualCommonEventType == commonEventType:
                plusEvents = binDb.traverse(index, windowSize, [plusEventType])
                # Append events from the adjacent left bin
                plusEvents.extend(binDb.collect_bin(windowSize, index-1, plusEventType))
                # Append events from the adjacent right bin
                plusEvents.extend(binDb.collect_bin(windowSize, index+1, plusEventType))

                #for eventP in plusEvents:
                    #print ('plusEvents ' + str(windowSize)  + ' ' +   str(index-1)  + ' ' +  plusEventType + ' ' + str([event.readName for event in eventP.events])+ '\t' + str(eventP.clusterId) +'\n')

                #if actualCommonEventType == commonEventType:
                minusEvents = binDb.traverse(index, windowSize, [minusEventType])
                # Append events from the adjacent left bin
                minusEvents.extend(binDb.collect_bin(windowSize, index-1, minusEventType))
                # Append events from the adjacent right bin
                minusEvents.extend(binDb.collect_bin(windowSize, index+1, minusEventType))

                #for eventM in minusEvents:
                    #print ('minusEvents ' + str(windowSize)  + ' ' +   str(index-1)  + ' ' +  minusEventType + ' ' + str([event.readName for event in eventM.events])+ '\t' + str(eventM.clusterId) +'\n')

                eventsMinusDict = {} 
                eventsPlusDict = {}

                # si ninguna de las dos listas esta vacia:
                if len(plusEvents) > 0 and len(minusEvents) > 0:
                    ### 2. Cluster events based on reciprocal overlap
                    ## For each event A
                    for eventPlus in plusEvents:
                        overlapDict = {}
                        eventsPlusDict[eventPlus] = {}

                        
                        ## [SR CHANGE]: Look at the previous event in order to not split a cluster
                        ## Identify events overlapping A (skip A itself and event pairs already assessed)
                        for eventMinus in minusEvents:

                            

                            ## Skip comparison if B belongs to a cluster already known to overlap A
                            #if (eventB.clusterId in clustersOverlapA):
                                #continue
                            
                            ## Add buffer to ranges
                            begA = eventPlus.beg - buffer
                            endA = eventPlus.end + buffer
                            begB = eventMinus.beg - buffer
                            endB = eventMinus.end + buffer
                            
                            overlap, overlapLen = gRanges.rcplOverlap(begA, endA, begB, endB, minPercOverlap)
                            #print ('overlap, overlapLen '+ str([event.readName for event in eventMinus.events]) + ' ' +str(overlap) +' '+ str(overlapLen))

                            # A) Event B overlap A. 
                            if overlap:
                                overlapDict[eventMinus] = overlapLen
                        # Una vez que tengo mirado, para un eventPlus, los eventMinus que lo solapan (lo he comparado con todos), cojo el que tenga mayor overlapLen (pq significa que estan mas cerca)
                        # Cojo la key del overlap dict cuyo value sea mayor (pq significa que overlapan mas, sera la pareja)
                        if overlapDict:
                            eventMinusReciprocal, eventMinusReciprocalOverlapLen = max(overlapDict.items(), key=operator.itemgetter(1))
                            #print ('overlap, overlapLen '+ str([event.readName for event in eventMinusReciprocal.events]) +' ' + str([event.readName for event in eventPlus.events]) + ' ' +str(overlap) +' '+ str(overlapLen))

                            eventsPlusDict[eventPlus][eventMinusReciprocal]=eventMinusReciprocalOverlapLen
                            #print ('eventsPlusDict')
                            #for k,v in eventsPlusDict.items():
                                #print (str([event.readName for event in k.events]))
                                #for l,m in v.items():
                                    #print (str([event.readName for event in l.events]))
                                    #print (m)

                    # VER EXPLANATIOOOOOOOOOOOOON
                    # Una vez que tengo mirados todos los eventPlus, cojo el mejor eventPlus para cada eventMinus.
                    # Para eso, le doy la vuelta al dictionary
                    # TAMBIEN HAY SOLUCION DE UNA SOLA LINEA!!!!! BUSCARLA EN MARCADORES!

                    for eventPlus, nestDict in eventsPlusDict.items():
                        for eventMinus, overlapLen in nestDict.items():
                            eventsMinusDict.setdefault(eventMinus, {}).update({eventPlus:overlapLen})
                    
                    # Y me quedo solo con aquellos eventsPlus que tengan la mayor overlapLen, asi ya tengo los pares
                    for eventMinus,nestDict in eventsMinusDict.items():
                        #print ('eventMinus.events')
                        #print ([event.readName for event in eventMinus.events])
                        #print ('nestDict')
                        #for key in nestDict.keys():
                            #print ([event.readName for event in key.events])
                        #reciprocalDict[eventMinus] = {}
                        if 'RECIPROCAL-' + commonEventType in reciprocalDict.keys():
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventMinus)
                        else:
                            reciprocalDict['RECIPROCAL-' + commonEventType] = []
                            #reciprocalDict[commonEventType].append(eventMinus)
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventMinus)
                        if len(nestDict) > 1:
                            eventPlus_to_add, overlapLen_to_add = max(nestDict.items(), key=operator.itemgetter(1))
                            #reciprocalDict[eventMinus][eventPlus_to_add] = overlapLen_to_add
                            #reciprocalDict[commonEventType].append(eventPlus_to_add)
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventPlus_to_add)
                        else:
                            #reciprocalDict[eventMinus] = nestDict
                            #reciprocalDict[commonEventType].append(*nestDict)
                            eventPlus = *nestDict,[0]
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventPlus[0])

                elif len(plusEvents) > 0:
                    if 'PLUS-' + commonEventType in reciprocalDict.keys():
                        reciprocalDict['PLUS-' + commonEventType].extend(plusEvents)
                    else:
                        reciprocalDict['PLUS-' + commonEventType] = []
                        reciprocalDict['PLUS-' + commonEventType].extend(plusEvents)

                elif len(minusEvents) > 0:
                    if 'MINUS-' + commonEventType in reciprocalDict.keys():
                        reciprocalDict['MINUS-' + commonEventType].extend(minusEvents)
                    else:
                        reciprocalDict['MINUS-' + commonEventType] = []
                        reciprocalDict['MINUS-' + commonEventType].extend(minusEvents)

    # Example of reciprocalDict
    # reciprocalDict = {'DISCORDANT-Hepatitis': [<clusters.DISCORDANT_cluster object at 0x7f32a07230b8>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>], 'DISCORDANT-UNVERIFIED:': [<clusters.DISCORDANT_cluster object at 0x7f32a07230b8>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>], 'DISCORDANT-HBV': [<clusters.DISCORDANT_cluster object at 0x7f32a07230f0>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>]}
    # Aqui podria retornar una lista de clusters, como aparece en el ejemplo de arriba, pero para hacer el metaclustering necesito events, asi que voy a retornar los events directamente
    '''
    for tipo, events in reciprocalDict.items():
        #print ('tipo ' + str(tipo))
        for event in events:
            #print ('side ' + str(event.side))
            #print ('type ' + str(event.type))
            #print ('identity ' + str(event.identity))
    '''
    #print (reciprocalDict)
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
    #print (result)
    # {'clusterMinus1': {'clusterPlus1': 576, 'clusterPlus2': 556}, 'clusterMinus3': {'clusterPlus3': 459, 'clusterPlus4': 20}, 'clusterMinus5': {'clusterPlus5': 20}}
    '''
    # Ahora si que podemos escoger, para cada cluster minus, cual es el mejor clusterPlus, y ya hacemos el dictionary final, que contendra solo los pares reciprocos, junto con la overlapLen
    '''
    final = {}

    for k,v in result.items():
        final[k] = {}
        if len(v) > 1:
            key_to_add, value_to_add = max(v.items(), key=operator.itemgetter(1))
            #print (key_to_add)
            #print (value_to_add)
            final[k][key_to_add] = value_to_add
        else:
            final[k] = v

    #print (final)
    # {'clusterMinus1': {'clusterPlus1': 576}, 'clusterMinus3': {'clusterPlus3': 459}, 'clusterMinus5': {'clusterPlus5': 20}}
    '''
    ## ----------------------------------------- EXPLANATIOOOOOOOOOOOOON -----------------------------------------

def reciprocal_DEPRECATED2(binDb, minPercOverlap, minClusterSize, buffer):
    '''
    Make reciprocal clusters from PLUS and MINUS discordant clusters. 
    Those clusters that have opposite orientation, are close (according to buffer parameter) and have the same identity are clustered together.

    Input:
        1. binDb: data structure containing a set of events organized in genomic bins 
        2. minPercOverlap: minimum percentage of reciprocal overlap to cluster two events together
        3. minClusterSize: minimum number of events clustering together for creating a cluster
        4. buffer: number of nucleotides used to extend cluster begin and end coordinates prior evaluating reciprocal overlap 
        
    Output:
        1. reciprocalDict: dictionary containing list of clusters grouped according to their type. The type is defined by ORIENTATION-EVENTSTYPE-IDENTITY (i.e. PLUS-DISCORDANT-HBV, RECIPROCAL-DISCORDANT-HBV, MINUS-DISCORDANT-HBV)
    '''

    eventsMinusDict = {} 
    eventsPlusDict = {}
    reciprocalDict = {}

    # Collect all event types present in the binDb
    commonEventTypes = set([i.split('-', 1)[1] for i in binDb.eventTypes])

    # For each window size/level
    for windowSize in binDb.binSizes:

        # For each bin in the current window size
        for index in binDb.data[windowSize]:

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
                ##print (actualCommonEventType)
                #plusEvents = []
                #minusEvents = []

                #if actualCommonEventType == commonEventType:
                plusEvents = binDb.traverse(index, windowSize, [plusEventType])
                # Append events from the adjacent left bin
                plusEvents.extend(binDb.collect_bin(windowSize, index-1, plusEventType))
                # Append events from the adjacent right bin
                plusEvents.extend(binDb.collect_bin(windowSize, index+1, plusEventType))

                #for eventP in plusEvents:
                    #print ('plusEvents ' + str(windowSize)  + ' ' +   str(index-1)  + ' ' +  plusEventType + ' ' + str([event.readName for event in eventP.events])+ '\t' + str(eventP.clusterId) + '\t' + str(eventP)+ '\t' + str(eventP.id) +'\n')

                #if actualCommonEventType == commonEventType:
                minusEvents = binDb.traverse(index, windowSize, [minusEventType])
                # Append events from the adjacent left bin
                minusEvents.extend(binDb.collect_bin(windowSize, index-1, minusEventType))
                # Append events from the adjacent right bin
                minusEvents.extend(binDb.collect_bin(windowSize, index+1, minusEventType))

                #for eventM in minusEvents:
                    #print ('minusEvents ' + str(windowSize)  + ' ' +   str(index-1)  + ' ' +  minusEventType + ' ' + str([event.readName for event in eventM.events])+ '\t' + str(eventM.clusterId) + '\t' + str(eventM)+ '\t' + str(eventM.id) +'\n')

                eventsMinusDict = {} 
                eventsPlusDict = {}

                # si ninguna de las dos listas esta vacia:
                if len(plusEvents) > 0 and len(minusEvents) > 0:
                    ### 2. Cluster events based on reciprocal overlap
                    ## For each event A
                    for eventPlus in plusEvents:
                        highestOverlapLen = 0
                        eventsPlusDict[eventPlus] = {}

                        ## [SR CHANGE]: Look at the previous event in order to not split a cluster
                        ## Identify events overlapping A (skip A itself and event pairs already assessed)
                        for eventMinus in minusEvents:

                            ## Skip comparison if B belongs to a cluster already known to overlap A
                            #if (eventB.clusterId in clustersOverlapA):
                                #continue
                            
                            ## Add buffer to ranges
                            begA = eventPlus.beg - buffer
                            endA = eventPlus.end + buffer
                            begB = eventMinus.beg - buffer
                            endB = eventMinus.end + buffer
                            
                            overlap, overlapLen = gRanges.rcplOverlap(begA, endA, begB, endB, minPercOverlap)
                            #print ('overlap, overlapLen '+ str([event.readName for event in eventMinus.events]) + ' ' +str(overlap) +' '+ str(overlapLen))

                            # A) Event B overlap A. 
                            if overlap:
                                if highestOverlapLen > 0:
                                    if overlapLen > highestOverlapLen:
                                        highestOverlapLen = overlapLen
                                        closestEventMinus = eventMinus
                                else:
                                    highestOverlapLen = overlapLen
                                    closestEventMinus = eventMinus

                        if highestOverlapLen > 0:
                            eventsPlusDict[eventPlus][closestEventMinus] = highestOverlapLen


                    # VER EXPLANATIOOOOOOOOOOOOON
                    # Una vez que tengo mirados todos los eventPlus, cojo el mejor eventPlus para cada eventMinus.
                    # Para eso, le doy la vuelta al dictionary
                    # TAMBIEN HAY SOLUCION DE UNA SOLA LINEA!!!!! BUSCARLA EN MARCADORES!

                    for eventPlus, nestDict in eventsPlusDict.items():
                        for eventMinus, overlapLen in nestDict.items():
                            eventsMinusDict.setdefault(eventMinus, {}).update({eventPlus:overlapLen})
                    
                    # Y me quedo solo con aquellos eventsPlus que tengan la mayor overlapLen, asi ya tengo los pares
                    for eventMinus,nestDict in eventsMinusDict.items():
                        #print ('eventMinus.events')
                        #print ([event.readName for event in eventMinus.events])
                        #print ('nestDict')
                        #for key in nestDict.keys():
                            #print ([event.readName for event in key.events])
                        #reciprocalDict[eventMinus] = {}
                        if 'RECIPROCAL-' + commonEventType in reciprocalDict.keys():
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventMinus)
                        else:
                            reciprocalDict['RECIPROCAL-' + commonEventType] = []
                            #reciprocalDict[commonEventType].append(eventMinus)
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventMinus)
                        if len(nestDict) > 1:
                            eventPlus_to_add, overlapLen_to_add = max(nestDict.items(), key=operator.itemgetter(1))
                            #reciprocalDict[eventMinus][eventPlus_to_add] = overlapLen_to_add
                            #reciprocalDict[commonEventType].append(eventPlus_to_add)
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventPlus_to_add)
                        else:
                            #reciprocalDict[eventMinus] = nestDict
                            #reciprocalDict[commonEventType].append(*nestDict)
                            eventPlus = *nestDict,[0]
                            reciprocalDict['RECIPROCAL-' + commonEventType].append(eventPlus[0])

                elif len(plusEvents) > 0:
                    if 'PLUS-' + commonEventType in reciprocalDict.keys():
                        reciprocalDict['PLUS-' + commonEventType].extend(plusEvents)
                    else:
                        reciprocalDict['PLUS-' + commonEventType] = []
                        reciprocalDict['PLUS-' + commonEventType].extend(plusEvents)

                elif len(minusEvents) > 0:
                    if 'MINUS-' + commonEventType in reciprocalDict.keys():
                        reciprocalDict['MINUS-' + commonEventType].extend(minusEvents)
                    else:
                        reciprocalDict['MINUS-' + commonEventType] = []
                        reciprocalDict['MINUS-' + commonEventType].extend(minusEvents)

    # Example of reciprocalDict
    # reciprocalDict = {'DISCORDANT-Hepatitis': [<clusters.DISCORDANT_cluster object at 0x7f32a07230b8>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>], 'DISCORDANT-UNVERIFIED:': [<clusters.DISCORDANT_cluster object at 0x7f32a07230b8>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>], 'DISCORDANT-HBV': [<clusters.DISCORDANT_cluster object at 0x7f32a07230f0>, <clusters.DISCORDANT_cluster object at 0x7f32a079ad68>]}
    # Aqui podria retornar una lista de clusters, como aparece en el ejemplo de arriba, pero para hacer el metaclustering necesito events, asi que voy a retornar los events directamente
    '''
    for tipo, events in reciprocalDict.items():
        #print ('tipo ' + str(tipo))
        for event in events:
            #print ('side ' + str(event.side))
            #print ('type ' + str(event.type))
            #print ('identity ' + str(event.identity))

    '''
    for key, value in reciprocalDict.items():
        reciprocalDict[key] = set(value)
    #print ('reciprocalDict ' + str(reciprocalDict))
    return reciprocalDict


def reciprocal(binDb, minPercOverlap, minClusterSize, buffer):
    '''
    Make reciprocal clusters from PLUS and MINUS discordant clusters. 
    Those clusters that have opposite orientation, are close (according to buffer parameter) and have the same identity are clustered together.

    Input:
        1. binDb: data structure containing a set of events organized in genomic bins 
        2. minPercOverlap: minimum percentage of reciprocal overlap to cluster two events together
        3. minClusterSize: minimum number of events clustering together for creating a cluster
        4. buffer: number of nucleotides used to extend cluster begin and end coordinates prior evaluating reciprocal overlap 
        
    Output:
        1. finalDict: dictionary containing list of clusters grouped according to their type. The type is defined by ORIENTATION-EVENTSTYPE-IDENTITY (i.e. PLUS-DISCORDANT-HBV, RECIPROCAL-DISCORDANT-HBV, MINUS-DISCORDANT-HBV)
    '''

    eventsPlusDict = {}
    reciprocalDict = {}
    plusDict = {}
    minusDict = {}

    # Collect all event types present in the binDb
    commonEventTypes = set([i.split('-', 1)[1] for i in binDb.eventTypes])

    # For each window size/level
    for windowSize in binDb.binSizes:

        # For each bin in the current window size
        for index in binDb.data[windowSize]:

            for commonEventType in commonEventTypes:
                # Set PLUS event type
                plusEventType =  'PLUS-' + commonEventType
                # Set MINUS event type
                minusEventType =  'MINUS-' + commonEventType

                # Collect PLUS events of current event type
                plusEvents = binDb.traverse(index, windowSize, [plusEventType])
                # Append events from the adjacent left bin
                plusEvents.extend(binDb.collect_bin(windowSize, index-1, plusEventType))
                # Append events from the adjacent right bin
                plusEvents.extend(binDb.collect_bin(windowSize, index+1, plusEventType))

                # Collect MINUS events of current event type
                minusEvents = binDb.traverse(index, windowSize, [minusEventType])
                # Append events from the adjacent left bin
                minusEvents.extend(binDb.collect_bin(windowSize, index-1, minusEventType))
                # Append events from the adjacent right bin
                minusEvents.extend(binDb.collect_bin(windowSize, index+1, minusEventType))

                # If both lists are filled:
                if len(plusEvents) > 0 and len(minusEvents) > 0:
                    ### 2. Cluster events based on reciprocal overlap
                    ## For each event A
                    for eventPlus in plusEvents:

                        ## Identify events overlapping A (skip A itself and event pairs already assessed)
                        for eventMinus in minusEvents:
                            
                            ## Add buffer to ranges
                            begA = eventPlus.beg - buffer
                            endA = eventPlus.end + buffer
                            begB = eventMinus.beg - buffer
                            endB = eventMinus.end + buffer
                            
                            overlap, overlapLen = gRanges.rcplOverlap(begA, endA, begB, endB, minPercOverlap)

                            # A) Event B overlap A.
                            # Dictionary containing eventsPlusDict[eventPlus][eventMinus] = overlapLen. Keeping only the pair with lowest overlapLen
                            if overlap:
                                try:
                                    for value in eventsPlusDict[eventPlus].values():    
                                        actualOverlapLen = int(value)
                                    if overlapLen > actualOverlapLen:
                                        eventsPlusDict[eventPlus] = {}
                                        eventsPlusDict[eventPlus][eventMinus] = overlapLen
                                except KeyError:
                                    eventsPlusDict[eventPlus] = {}
                                    eventsPlusDict[eventPlus][eventMinus] = overlapLen


    # If eventMinus is repeated in previous dictionary (paired with different events plus), pick the pair with the highest overlapLen:
    for eventPlus, nestDict in eventsPlusDict.items():
        maximumEventMinus = max(nestDict, key=nestDict.get)
        # Save the best pair in reciprocalDict:
        if 'RECIPROCAL-' + str(eventPlus.clusterType) + '-' + str(eventPlus.identity) in reciprocalDict.keys():
            reciprocalDict['RECIPROCAL-' + str(eventPlus.clusterType) + '-' + str(eventPlus.identity)].append(maximumEventMinus)
            reciprocalDict['RECIPROCAL-' + str(eventPlus.clusterType) + '-' + str(eventPlus.identity)].append(eventPlus)
        else:
            reciprocalDict['RECIPROCAL-' + str(eventPlus.clusterType) + '-' + str(eventPlus.identity)] = []
            reciprocalDict['RECIPROCAL-' + str(eventPlus.clusterType) + '-' + str(eventPlus.identity)].append(maximumEventMinus)
            reciprocalDict['RECIPROCAL-' + str(eventPlus.clusterType) + '-' + str(eventPlus.identity)].append(eventPlus)

    # Collect those plus and minus events that have no reciprocal cluster:
    plusEventTypes = []
    minusEventTypes = []
    for commonEventType in commonEventTypes:
        plusEventTypes.append('PLUS-' + commonEventType)
        minusEventTypes.append('MINUS-' + commonEventType)

    for plusEvent in binDb.collect(plusEventTypes):
        if not any(plusEvent in value1 for value1 in reciprocalDict.values()) and not any(plusEvent in value2 for value2 in plusDict.values()):
            if 'PLUS-' + str(plusEvent.clusterType) + '-' + str(plusEvent.identity) in plusDict.keys():
                plusDict['PLUS-' + str(plusEvent.clusterType) + '-' + str(plusEvent.identity)].append(plusEvent)
            else:
                plusDict['PLUS-' + str(plusEvent.clusterType) + '-' + str(plusEvent.identity)] = []
                plusDict['PLUS-' + str(plusEvent.clusterType) + '-' + str(plusEvent.identity)].append(plusEvent)

    for minusEvent in binDb.collect(minusEventTypes):
        if not any(minusEvent in value3 for value3 in reciprocalDict.values()) and not any(minusEvent in value4 for value4 in minusDict.values()):
            if 'MINUS-' + str(minusEvent.clusterType) + '-' + str(minusEvent.identity) in minusDict.keys():
                minusDict['MINUS-' + str(minusEvent.clusterType) + '-' + str(minusEvent.identity)].append(minusEvent)
            else:
                minusDict['MINUS-' + str(minusEvent.clusterType) + '-' + str(minusEvent.identity)] = []
                minusDict['MINUS-' + str(minusEvent.clusterType) + '-' + str(minusEvent.identity)].append(minusEvent)
                        

    # Merge the three dictionaries:
    plusDict.update(minusDict)
    plusDict.update(reciprocalDict)

    finalDict = plusDict

    return finalDict

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

                    ## TODO SR: Remove this KeyError message at reciprocal_overlap_clustering when debugging is done!!!
                    try:
                        clustersDict[clusterId].add(events2Cluster)
                    except KeyError:
                        time.sleep(60)
                        sys.exit('ERROR 13. clusterId = clustersOverlapA[0] ' +str(clusterId) +' ' + str(eventsInClusters) +' '+ str(idx) +' '+ str(eventA) +' '+ str(events) +' '+ str(index) +' '+ str(windowSize) +' '+ str(eventsInClusters) +' '+ str(clustersDict) +' '+ str(os.getpid()))


                # B) Multiple clusters overlap A -> Merge clusters and add A and its overlapping events into the merged cluster
                elif len(clustersOverlapA) > 1:

                    ## Make list of clusters overlapping A
                    ## TODO SR: Remove this KeyError message at reciprocal_overlap_clustering when debugging is done!!!
                    try:
                        clusters2merge = [ clustersDict[clusterId] for clusterId in clustersOverlapA ]
                    except KeyError:
                        time.sleep(60)
                        sys.exit('ERROR 15. clusters2merge ' + str(clusters2merge) +' ' + str(eventsInClusters) +' '+ str(idx) +' '+ str(eventA) +' '+ str(events) +' '+ str(index) +' '+ str(windowSize) +' '+ str(eventsInClusters) +' '+ str(clustersDict) +' '+ str(os.getpid()))

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
                    for cluster in clusters2merge:                        
                        # TODO SR: Remove this piece of commented code
                        '''
                        FUNCIONO
                        if 'META' in cluster.id:
                            cluster.id = mergedCluster.id
                        
                        if 'CLUSTER_' in cluster.id:
                            #print ('AAAAAA ' + cluster.id)
                            cluster.clusterId = mergedCluster.id
                        '''
                        

                        ## Sanity check. If cluster is not in the dict, do not attempt to remove it:
                        
                        if cluster.id in clustersDict:
                            clustersDict.pop(cluster.id, None)
                        
                        #if cluster.id in clustersOverlapA:
                            #clustersOverlapA.remove(cluster.id)
                        
                        # NOTE SR: NEEDED TO NOT CRASH!!!!!
                        if 'META' in mergedCluster.id:
                            # TODO SR: Remove this piece of commented code, as it was copied to clusters.py
                                #clusterNew.clusterId = mergedCluster.id
                            if cluster.id in clustersOverlapA:
                                clustersOverlapA.remove(cluster.id)
                                if mergedCluster.id not in clustersOverlapA:
                                    clustersOverlapA.append(mergedCluster.id)
                            # TODO SR: Remove this commented code
                            #del cluster
	
                    # TODO SR: Remove this piece of commented code
                    '''
                    if 'META' in mergedCluster.id:
                        for clusterNew in mergedCluster.subclusters.values():
                            clusterNew.clusterId = mergedCluster.id
                        for evento in mergedCluster.events:
                            evento.clusterId = mergedCluster.id
                    elif 'CLUSTER_' in mergedCluster.id:
                        for evento in mergedCluster.events:
                            evento.clusterId = mergedCluster.id
                    '''
                    #clustersOverlapA.append(mergedCluster.id)


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
