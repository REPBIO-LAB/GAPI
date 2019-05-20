'''
Module 'clustering' - Contains functions for clustering sets of objects based on coordinates
'''

## DEPENDENCIES ##
# External
import time
import sys
import operator

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

## SR CHANGE
def reciprocal_clustering(eventsBinDb, minPercOverlap, minClusterSize, eventType, buffer, clusterType):
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
            events = eventsBinDb.traverse(index, windowSize, [eventType])

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

                ## [SR CHANGE]: Look at the previous event in order to not split a cluster
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

                        # [SR CHANGE]
                        # Create cluster                        
                        cluster = clusters.create_cluster(events2Cluster, clusterType)
                        clustersDict[cluster.id] = cluster

                    # Cluster not composed by enough number of events
    
    clustersList = list(clustersDict.values())

    return clustersList


def reciprocal(eventsBinDb, minPercOverlap, minClusterSize, buffer):
    '''
    '''
    eventsMinusDict = {} 
    eventsPlusDict = {}
    reciprocalDict = {}

    # For each window size/level
    for windowSize in eventsBinDb.binSizes:

        # For each bin in the current window size
        for index in eventsBinDb.data[windowSize]:

            # Collect all event types present in the binDb
            commonEventTypes = set([i.split('-', 1)[1] for i in eventsBinDb.eventTypes])

            for commonEventType in commonEventTypes:
                # cogo solo las partes comunes de los eventTypes
                #commonEventType =  '-'.join(eventType.split('-')[1:])
                # Hago una lista con los clusters de ese tipo que sean PLUS
                # Hago una lista con los clusters de ese tipo que sean MINUS
                plusEventType =  'PLUS-' + commonEventType
                minusEventType =  'MINUS-' + commonEventType

                #for actualEventType in eventsBinDb.data[windowSize][index]:
                # cogo solo las partes comunes de los eventTypes
                #actualCommonEventType =  '-'.join(actualEventType.split('-')[1:])
                #print (actualCommonEventType)
                #plusEvents = []
                #minusEvents = []

                #if actualCommonEventType == commonEventType:
                plusEvents = eventsBinDb.traverse(index, windowSize, [plusEventType])
                # Append events from the adjacent left bin
                plusEvents.extend(eventsBinDb.collect_bin(windowSize, index-1, plusEventType))
                # Append events from the adjacent right bin
                plusEvents.extend(eventsBinDb.collect_bin(windowSize, index+1, plusEventType))

                #if actualCommonEventType == commonEventType:
                minusEvents = eventsBinDb.traverse(index, windowSize, [minusEventType])
                # Append events from the adjacent left bin
                minusEvents.extend(eventsBinDb.collect_bin(windowSize, index-1, minusEventType))
                # Append events from the adjacent right bin
                minusEvents.extend(eventsBinDb.collect_bin(windowSize, index+1, minusEventType))

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
