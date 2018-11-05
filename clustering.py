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
import structures

## FUNCTIONS ##
def clusterByDist1D(binDb, maxDist, minClusterSize, clusterType):
    '''
    Cluster events (CLIPPING, INS...) by begin positions. 

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
    binSize = binDb.binSizes[0]
    
    ## 1. Create clusters ##
    # a) Fast clustering algorithm
    if (binSize <= maxDist):
        clusters = clusterByDist1D_fast(binDb.data[binSize], maxDist, minClusterSize, clusterType)

    # b) Iterative clustering algorithm (not implemented yet)
    else:
        print('ITERATIVE ALGORITHM. NOT IMPLEMENTED YET...')
        # clusterByDist1D_iterative(binDb.data[binSize], maxDist, minClusterSize, clusterType)
    

    ## 2. Organize clusters into bins ##    
    binSizes = [100, 1000, 10000, 100000, 1000000]
    data = [(clusters, clusterType + '-CLUSTER')]
    clustersBins = structures.createBinDb(data, binSizes)

    return clustersBins


def clusterByDist1D_fast(binHash, maxDist, minClusterSize, clusterType):
    '''
    '''
    clusters = []
    binsInClusters = []

    # For each bin 
    for binIndex in binHash:

        # There are elements of the target cluster type in the bin
        if clusterType in binHash[binIndex]:

            ### Number of events in the bin
            binObj = binHash[binIndex][clusterType]

            ### Skip bins fulfilling one of these conditions:
            # a) Bin already incorporated into a cluster 
            if binIndex in binsInClusters:
                continue
    
            # b) Bin without enough number of events to build a root cluster
            elif (binObj.nbEvents() < minClusterSize):
                continue

            ### Create root cluster containing events on the bin
            cluster = variants.createCluster(binObj.events, clusterType)
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
                    newBinObj = binHash[backwardIndex][clusterType]
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
                if (forwardIndex in binHash) and (clusterType in binHash[forwardIndex]):

                    # Compute the distance between the left most position in the new bin
                    # and the cluster end 
                    newBinObj = binHash[forwardIndex][clusterType]
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


