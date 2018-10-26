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
def clusterCLIPPING(CLIPPING_list, maxBkpDist, minClusterSize):
    '''
    Cluster CLIPPING events based on breakpoint position coordinates

    Input:
        1. CLIPPING_list: list of CLIPPING objects
        2. maxBkpDist: maximum distance between two breakpoints to include them into the same cluster
        3. minClusterSize: minimum number of clippings required to build a root cluster in a window

    Output:
        1. clusterDict: positions dictionary containing the clusters
    ''' 
    step = 'CLUSTER-CLIPPING'
    msg = 'Input (CLIPPING_list, maxBkpDist, minClusterSize): ' +  "\t".join([str(len(CLIPPING_list)), str(maxBkpDist), str(minClusterSize)])    
    log.step(step, msg)

    ## 1. Organize CLIPPING events by their breakpoint position ##
    posDict = structures.buildPosDict(CLIPPING_list, maxBkpDist)

    ## 2. Cluster CLIPPING events ##
    # Initialize list with windows already incorporated into clusters
    windowsInClusters = []
    clusterList = []

    # For each genomic window
    for windowIndex, CLIPPINGS in posDict.items():
    
        ### Number of CLIPPING events in the window
        nbCLIPPINGS = len(CLIPPINGS)

        ### Skip window fulfilling one of these conditions:
        # a) Window already incorporated into a cluster 
        if windowIndex in windowsInClusters:
            continue
    
        # b) Window without enough number of CLIPPING events to build a root cluster
        elif (nbCLIPPINGS < minClusterSize):
            continue

        ### Create root cluster containing CLIPPINGS on the window
        cluster = variants.CLIPPING_cluster(CLIPPINGS)
        clusterList.append(cluster)
        windowsInClusters.append(windowIndex) # now window incorporated into cluster

        ### Root cluster extension
        ## 2.1 Go backward from current window (*).  
        #       <---2--- <---1---
        # |---------|--------|----*----|---------
        # Go one window backward in each iteration. Extend the cluster 
        # if last clipping breakpoint in previous window within 
        # maximum cluster distance or end iteration, otherwise 
        backwardIndex = windowIndex - 1

        while True:
        
            # A) There are CLIPPINGS in the window 
            if backwardIndex in posDict:
                
                # Compute the distance between the right most CLIPPING in the new window
                # and the cluster begin 
                CLIPPINGS = posDict[backwardIndex]
                lastCLIPPING = CLIPPINGS[-1]
                bkpDist = cluster.pos - lastCLIPPING.pos

                # a) Last CLIPPING within maximum distance 
                if bkpDist <= maxBkpDist:

                    # Add CLIPPINGS within window to the cluster 
                    cluster.add(CLIPPINGS, 'left')
                    windowsInClusters.append(backwardIndex) # now window incorporated into cluster
         
                # b) CLIPPING outside 
                else:
                    break

            # B) No CLIPPINGS in the window. Stop iterating
            else:
                break 
            
            backwardIndex -= 1
        
        ## 2.2 Go forward from current window (*).  
        #                ---1---> ---2--->
        # |---------|----*----|--------|---------
        # Go one window forward in each iteration. Extend the cluster 
        # if first clipping breakpoint in next window within 
        # maximum cluster distance or end iteration, otherwise 
        forwardIndex = windowIndex + 1

        while True:
        
            # A) There are CLIPPINGS in the window 
            if forwardIndex in posDict:
                
                # Compute the distance between the left most CLIPPING in the new window
                # and the cluster end 
                CLIPPINGS = posDict[forwardIndex]
                firstCLIPPING = CLIPPINGS[0]
                bkpDist = firstCLIPPING.pos - cluster.end 

                # a) Last CLIPPING within maximum distance 
                if bkpDist <= maxBkpDist:

                    # Add CLIPPINGS within window to the cluster 
                    cluster.add(CLIPPINGS, 'right')
                    windowsInClusters.append(forwardIndex) # now window incorporated into cluster

                # b) CLIPPING outside 
                else:
                    break

            # B) No CLIPPINGS in the window. Stop iterating
            else:
                break 
            
            forwardIndex += 1        

    ## 3. Organize CLIPPING clusters into a dictionary ##
    clustersDict = structures.buildPosDict(clusterList, 1000)
    nbClusters = len(clusterList)

    return clustersDict, nbClusters

## Pending, create these functions:
#clusterINS
#clusterDEL

## CLASSES ##
