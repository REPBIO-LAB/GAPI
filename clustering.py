'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import time

# Internal
import log
import structures
import variants

## FUNCTIONS ##
def clusterByPos(events, maxDist, minClusterSize):
    '''
    Cluster events (CLIPPING, INS...) by breakpoint position coordinates. Organize data into a dictionary with genomic windows for clustering

    Input:
        1. events: list of objects. Every object must have a 'pos' argument (position)
        2. maxDist: maximum distance between two positions to include them into the same cluster
        3. minClusterSize: minimum number of events required to build a root cluster in a window

    Output:
        1. clusterDict: dictionary containing the clusters organized in genomic windows
    ''' 

    ## 1. Organize events by their breakpoint position in genomic windows ##
    eventsDict = structures.mapEvents2windows(events, 'pos', maxDist)

    ## 2. Cluster events ##
    # Initialize list with windows already incorporated into clusters
    windowsInClusters = []
    clusterList = []

    # For each genomic window
    for windowIndex, events in eventsDict.items():
    
        ### Number of events in the window
        nbEvents = len(events)

        ### Skip windows fulfilling one of these conditions:
        # a) Window already incorporated into a cluster 
        if windowIndex in windowsInClusters:
            continue
    
        # b) Window without enough number of events to build a root cluster
        elif (nbEvents < minClusterSize):
            continue

        ### Create root cluster containing events on the window
        cluster = variants.cluster(events)
        clusterList.append(cluster)
        windowsInClusters.append(windowIndex) # now window incorporated into cluster

        ### Root cluster extension
        ## 2.1 Go backward from current window (*).  
        #       <---2--- <---1---
        # |---------|--------|----*----|---------
        # Go one window backward in each iteration. Extend the cluster 
        # if last breakpoint position in previous window within 
        # maximum cluster distance or end iteration, otherwise 
        backwardIndex = windowIndex - 1

        while True:
        
            # A) There are breakpoints in the window 
            if backwardIndex in eventsDict:
                
                # Compute the distance between the right most breakpoint position
                # in the new window and the cluster begin 
                events = eventsDict[backwardIndex]
                lastEvent = events[-1]
                posDist = cluster.beg - lastEvent.pos

                # a) Last event within maximum distance 
                if posDist <= maxDist:

                    # Add events within window to the cluster 
                    cluster.add(events, 'left')
                    windowsInClusters.append(backwardIndex) # now window incorporated into cluster
         
                # b) Event outside 
                else:
                    break

            # B) No events in the window. Stop iterating
            else:
                break 
            
            backwardIndex -= 1
        
        ## 2.2 Go forward from current window (*).  
        #                ---1---> ---2--->
        # |---------|----*----|--------|---------
        # Go one window forward in each iteration. Extend the cluster 
        # if first breakpoint position in next window within 
        # maximum cluster distance or end iteration, otherwise 
        forwardIndex = windowIndex + 1

        while True:
        
            # A) There are breakpoints in the window 
            if forwardIndex in eventsDict:
                
                # Compute the distance between the left most position in the new window
                # and the cluster end 
                events = eventsDict[forwardIndex]
                firstEvent = events[0]
                posDist = firstEvent.pos - cluster.end 

                # a) Last event within maximum distance 
                if posDist <= maxDist:

                    # Add events within window to the cluster 
                    cluster.add(events, 'right')
                    windowsInClusters.append(forwardIndex) # now window incorporated into cluster

                # b) Events outside 
                else:
                    break

            # B) No events in the window. Stop iterating
            else:
                break 
            
            forwardIndex += 1        

    ## 3. Organize clusters into a dictionary ##
    clustersDict = structures.mapEvents2windows(clusterList, 'beg', 1000)
    nbClusters = len(clusterList)

    return clustersDict, nbClusters

def makeMetaClusters(events, maxDist, minClusterSize):
    '''
    Group clusters into metaclusters. D
    Cluster events (CLIPPING, INS, DEL...) by breakpoint position coordinates. Organize data into a dictionary with genomic windows for clustering

    Input:
        1. events: list of objects. Every object must have a 'pos' argument (position)
        2. maxDist: maximum distance between two positions to include them into the same cluster
        3. minClusterSize: minimum number of events required to build a root cluster in a window

    Output:
        1. clusterDict: dictionary containing the clusters organized in genomic windows
    ''' 

## Pending, create these functions:
#clusterDEL

## CLASSES ##
