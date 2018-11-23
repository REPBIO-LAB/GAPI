'''
Module 'filters' - Contains functions for filtering clusters
'''

###############
## FUNCTIONS ##
###############

def filterClusters(clusters, clusterType, confDict):
    '''
    Function to apply filters to each cluster. It does not produce any output just modify cluster objects attribute filters.

    Input:
        1. clusters: bin database containing a set of cluster objects
        2. clusterType: type of cluster (INS-CLUSTER: insertion; DEL-CLUSTER: deletion; LEFT-CLIPPING-CLUSTER: left clipping; RIGHT-CLIPPING-CLUSTER: right clipping)
    '''

    ## For each cluster
    for cluster in clusters.collect(clusterType):
        
        ## Filter
        cluster.filters = filter(cluster,confDict)


def filter(cluster,confDict):
    '''
    Apply filters different filters to each cluster.

    Input:
        1. cluster: cluster object
        2. confDict
    Output:
        1. filtersResults -> keys: name of filters; values: True if the cluster pass the filter, False if it doesn't pass. 
    '''

    filtersResults = {}

    ## 0. Get a list containing the filters to apply
    filters2Apply = confDict['clusterFilters'].split(',')

    ## 1. FILTER 1: Minimum number of reads per cluster
    ## 1.0 Check if the filter is selected
    if "NBREADS" in filters2Apply:

        ## 1.1 Compute number of events supporting the cluster 
        nbEvents = cluster.nbEvents()[0]

        ## 1.2 Compare the number of events supporting the cluster against the minimum required
        minNbEvents = confDict['minClusterSize']

        if nbEvents >= minNbEvents:
            nbEventsFilter = True
        else:
            nbEventsFilter = False
        
        filtersResults["NBREADS"] = nbEventsFilter

    ## 2. FILTER 2: Maximum Coefficient of Variance per cluster
    ## 2.0 Check if the filter is selected and if it makes sense to apply the filter according to cluster type
    if "CV" in filters2Apply:

        # Check if length attribute is available for each event
        lengthsBool = [hasattr(event, 'length') for event in cluster.events]

        # Attemp polishing if length attribute available for all the events 
        if False not in lengthsBool:

            ## 2.1 Compute CV of the cluster 
            mean, std, cv = cluster.meanLen()
            ## 2.2 Compare the cluster CV against the maximum required
            maxClusterCV = confDict['maxClusterCV']

            if cv != "NA" and cv <= maxClusterCV:
                cvFilter = True
            else:
                cvFilter = False

            filtersResults["CV"] = cvFilter


    ## 3. FILTER 3: Maximum number of removed outliers per cluster
    ## 3.0 Check if the filter is selected and if it makes sense to apply the filter according to cluster type
    if "OUTLIERS" in filters2Apply:

        # Check if length attribute is available for each event
        lengthsBool = [hasattr(event, 'length') for event in cluster.events]

        # Attemp polishing if length attribute available for all the events 
        if False not in lengthsBool:

            ## 3.1 Compute the percentage of outliers: nbOutliers / nb events before polishing the cluster
            nbEventsUnpolished = cluster.nbOutliers + cluster.nbEvents()[0]
            percentageOutliers = cluster.nbOutliers / nbEventsUnpolished

            ## 3.1 Compare the percentage of outliers against the maximum required
            maxOutliers = confDict['maxOutliers']

            if percentageOutliers <= maxOutliers:
                outliersFilter = True
            else:
                outliersFilter = False

            filtersResults["OUTLIERS"] = outliersFilter

    
    return filtersResults