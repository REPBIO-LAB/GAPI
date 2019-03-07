'''
Module 'filters' - Contains functions for filtering clusters
'''

###############
## FUNCTIONS ##
###############

def filterClusters(clusters, clusterType, confDict):
    '''
    Function to apply filters to each cluster according to its type. It does not produce any output just modify cluster's attribute filters.

    Input:
        1. clusters: bin database containing a set of cluster objects
        2. clusterType: type of cluster (INS-CLUSTER: insertion; DEL-CLUSTER: deletion; LEFT-CLIPPING-CLUSTER: left clipping; RIGHT-CLIPPING-CLUSTER: right clipping)
        3. confDict
    '''

    ## Get a list containing the filters to apply
    filters2Apply = confDict['clusterFilters'].split(',')

    ## For each cluster
    for cluster in clusters.collect(clusterType):
        
        ## a) Filter INS cluster
        if (clusterType == 'INS-CLUSTER'):
            cluster.filters = filterINS(cluster, filters2Apply, confDict)

        ## b) Filter DEL cluster
        elif (clusterType == 'DEL-CLUSTER'):
            cluster.filters = filterDEL(cluster, filters2Apply, confDict)

        ## c) Filter CLIPPING cluster
        elif (clusterType == 'LEFT-CLIPPING-CLUSTER') or (clusterType == 'RIGHT-CLIPPING-CLUSTER'):
            cluster.filters = filterCLIPPING(cluster, filters2Apply, confDict)

        ## d) Unexpected cluster type
        else:
            log.info('Error at \'filterClusters\'. Unexpected cluster type')
            sys.exit(1)

def filterINS(cluster, filters2Apply, confDict):
    '''
    Apply appropriate filters to each INS-CLUSTER.

    Input:
        1. cluster: cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict
    Output:
        1. filterInsResults -> keys: name of filters; values: True if the cluster pass the filter, False if it doesn't pass.
    '''
    filterInsResults = {}

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: # check if the filter is selected
        filterInsResults['MIN-NBREADS'] = minNbEventsFilter(cluster,  confDict['minClusterSize'])

    ## 2. FILTER 2: Maximum number of reads per cluster
    if 'MAX-NBREADS' in filters2Apply: # check if the filter is selected
        filterInsResults['MAX-NBREADS'] = maxNbEventsFilter(cluster,  confDict['maxClusterSize'])

    ## 3. FILTER 3: Maximum Coefficient of Variance per cluster
    if 'CV' in filters2Apply: # check if the filter is selected
        filterInsResults['CV'] = maxCvFilter(cluster, confDict['maxClusterCV'])

    ## 4. FILTER 4: Maximum percentage of removed outliers per cluster
    if 'OUTLIERS' in filters2Apply: # check if the filter is selected
        filterInsResults['OUTLIERS'] = maxPercOutliers(cluster, confDict['maxOutliers'])

    return filterInsResults

def filterDEL(cluster, filters2Apply, confDict):
    '''
    Apply appropriate filters to each DEL-CLUSTER.

    Input:
        1. cluster: cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict
    Output:
        1. filterDelResults -> keys: name of filters; values: True if the cluster pass the filter, False if it doesn't pass.
    '''

    filterDelResults = {}

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: # check if the filter is selected
        filterDelResults['MIN-NBREADS'] = minNbEventsFilter(cluster, confDict['minClusterSize'])

    ## 2. FILTER 2: Maximum number of reads per cluster
    if 'MAX-NBREADS' in filters2Apply: # check if the filter is selected
        filterDelResults['MAX-NBREADS'] = maxNbEventsFilter(cluster, confDict['maxClusterSize'])

    ## 3. FILTER 3: Maximum Coefficient of Variance per cluster
    if 'CV' in filters2Apply: # check if the filter is selected
        filterDelResults['CV'] = maxCvFilter(cluster, confDict['maxClusterCV'])

    ## 4. FILTER 4: Maximum number of removed outliers per cluster
    if 'OUTLIERS' in filters2Apply: # check if the filter is selected
        filterDelResults['OUTLIERS'] = maxPercOutliers(cluster, confDict['maxOutliers'])

    return filterDelResults

def filterCLIPPING(cluster, filters2Apply, confDict):
    '''
    Apply appropriate filters to each CLIPPING-CLUSTER.

    Input:
        1. cluster: cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict
    Output:
        1. filterClippingResults -> keys: name of filters; values: True if the cluster pass the filter, False if it doesn't pass.
    '''

    filterClippingResults = {}

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: # check if the filter is selected
        filterClippingResults['MIN-NBREADS'] = minNbEventsFilter(cluster, confDict['minClusterSize'])

    ## 2. FILTER 2: Maximum number of reads per cluster
    if 'MAX-NBREADS' in filters2Apply: # check if the filter is selected
        filterClippingResults['MAX-NBREADS'] = maxNbEventsFilter(cluster, confDict['maxClusterSize'])

    return filterClippingResults

def minNbEventsFilter(cluster, minNbEvents):
    '''
    Filter cluster by comparing the number of cluster supporting events with a minimum treshold

    Input:
        1. cluster: cluster object
        2. minNbEvents: min number of events threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute number of events supporting the cluster 
    nbEvents = cluster.nbEvents()[0]

    ## 2. Compare the number of events supporting the cluster against the minimum required
    if nbEvents >= minNbEvents:
        PASS = True
    else:
        PASS = False
    
    return PASS

def maxNbEventsFilter(cluster, maxNbEvents):
    '''
    Filter cluster by comparing the number of cluster supporting events with a maximum treshold

    Input:
        1. cluster: cluster object
        2. maxNbEvents: maximum number of events threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute number of events supporting the cluster 
    nbEvents = cluster.nbEvents()[0]

    ## 2. Compare the number of events supporting the cluster against the maximum required
    if nbEvents <= maxNbEvents:
        PASS = True
    else:
        PASS = False
    
    return PASS


def maxCvFilter(cluster, maxClusterCV):
    '''
    Filter cluster by comparing its Coefficient of Variation with a maximum threshold.

    Input:
        1. cluster: cluster object
        2. maxClusterCV: maximum Coefficient of Variation threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute CV of the cluster 
    mean, std, cv = cluster.meanLen()

    ## 2. Compare the cluster CV against the maximum required
    if cv <= maxClusterCV:
        PASS = True
    else:
        PASS = False

    return PASS

def maxPercOutliers(cluster, maxOutliers):
    '''
    Filter cluster by comparing its percentage of outliers with a maximum threshold.

    Input:
        1. cluster: cluster object
        2. maxOutliers: maximum percentage of outlier events threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute the percentage of outliers: nbOutliers / nb events before polishing the cluster
    nbEventsUnpolished = cluster.nbOutliers + cluster.nbEvents()[0]
    percOutliers = cluster.nbOutliers / nbEventsUnpolished

    ## 2. Compare the percentage of outliers against the maximum required
    if percOutliers <= maxOutliers:
        PASS = True
        
    else:
        PASS = False

    return PASS