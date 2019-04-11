'''
Module 'filters' - Contains functions for filtering clusters
'''

###############
## FUNCTIONS ##
###############

## [SR CHANGE]: bam files added as argument
## TODO: ADD NORMAL BAM!!!
def filterClusters(clusters, clusterType, confDict, tumourBam):
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

        ## [SR CHANGE]
        ## TODO: CHANGE THIS EVENTYPE!!!!!!
        ## c) Filter DISCORDANT cluster
        elif 'DISCORDANT' in clusterType:
            cluster.filters = filterDISCORDANT(cluster, filters2Apply, confDict, tumourBam)

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

## [SR CHANGE]
def filterDISCORDANT(cluster, filters2Apply, confDict, bam):
    '''
    Apply appropriate filters to each DISCORDANT-CLUSTER.

    Input:
        1. cluster: cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict
    Output:
        1. filterDiscordantResults -> keys: name of filters; values: True if the cluster pass the filter, False if it doesn't pass.
    '''

    filterDiscordantResults = {}

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: # check if the filter is selected
        filterDiscordantResults['MIN-NBREADS'] = minNbEventsFilter(cluster, confDict['minClusterSize'])

    ## 2. FILTER 2: Maximum number of reads per cluster
    if 'MAX-NBREADS' in filters2Apply: # check if the filter is selected
        filterDiscordantResults['MAX-NBREADS'] = maxNbEventsFilter(cluster, confDict['maxClusterSize'])

    ## 3. FILTER 3: Area mapping quality
    if "AREAMAPQ" in filters2Apply:
        filterDiscordantResults["AREAMAPQ"] = area(cluster,confDict,bam)[0]

    ## 4. FILTER 4: Area clipping SMS
    if "AREASMS" in filters2Apply:
        filterDiscordantResults["AREASMS"] = area(cluster,confDict,bam)[1]

    return filterDiscordantResults

# [SR CHANGE]
def area(cluster,confDict,bam):
    '''
    Apply filters to cluster (SR bam) based on the characteristics of its region.

    Input:
        1. cluster: cluster object
        2. confDict
        3. bam
    Output:
        1. percMAPQFilter -> boolean: True if the cluster pass the filter, False if it doesn't
        2. percSMSReadsFilter -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## Set region coordinates
    if cluster.beg == cluster.end:
        if cluster.beg > 100:
            binBeg = cluster.beg - 100
        else:
            binBeg = cluster.beg
        binEnd = cluster.beg + 100
    else:
        if cluster.beg > 50:
            binBeg = cluster.beg - 50
        else:
            binBeg = cluster.beg
        binEnd = cluster.end + 50

    ref = cluster.ref

    ## Extract filter parameters from config dict
    minReadsRegionMQ = confDict['minReadsRegionMQ']
    maxRegionlowMQ = confDict['maxRegionlowMQ']
    maxRegionSMS = confDict['maxRegionSMS']

    # Set counts to 0
    lowMAPQ = 0
    SMSReads = 0
    nbReads = 0

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)

    for alignmentObj in iterator:
        
        # Check if aligment pass minimum mapq for reads within the cluster region
        passMAPQ = areaMAPQ(alignmentObj, minReadsRegionMQ)

        # If it doesnt pass, add 1 to the counts of low mapping quality reads within the cluster region
        if passMAPQ == False:
            lowMAPQ += 1

        # Check if aligment is mapped this way: Soft Match Soft (SMS)
        SMSRead = areaSMS(alignmentObj)
        
        # If it is mapped SMS, add 1 to the counts of SMS reads within the cluster region
        if SMSRead == True:
            SMSReads += 1
        # Count total number of reads in the region
        nbReads += 1
    
    ## Calculate percentages
    percMAPQ = percentage(lowMAPQ, nbReads)
    percSMSReads = percentage(SMSReads, nbReads)

    ## If the percentage of low MQ reads is lower than the threshold pass the filter.
    if percMAPQ < float(maxRegionlowMQ):
        percMAPQFilter = True
    else:
        percMAPQFilter = False

    ## If the percentage of SMS reads is lower than the threshold pass the filter.
    if percSMSReads < float(maxRegionSMS):
        percSMSReadsFilter = True
    else:
        percSMSReadsFilter = False
    
    #print (str(cluster.ref) + ' ' + str(cluster.beg) + ' ' + str (cluster.end) + ' ' + str(percMAPQFilter) + ' ' + str(percSMSReadsFilter))

    return percMAPQFilter, percSMSReadsFilter

# [SR CHANGE]
def areaMAPQ(alignmentObj, minReadsRegionMQ):
    '''
    Check if the MAPQ of a read pass the minReadsRegionMQ threshold

    Input:
        1. alignmentObj
        2. minReadsRegionMQ
    Output:
        1. percMAPQFilter -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    MAPQ = int(alignmentObj.mapping_quality)

    if MAPQ > int(minReadsRegionMQ):
        passMAPQ = True
    else:
        passMAPQ = False

    return passMAPQ

# [SR CHANGE]
def areaSMS(alignmentObj):
    '''
    Check if aligment is mapped this way: Soft Match Soft (SMS)

    Input:
        1. alignmentObj
        2. maxRegionSMS
    Output:
        1. percSMSReadsFilter -> boolean: True if the cluster pass the filter, False if it doesnt
    '''

    SMSRead = False

    # Select first and last operation from cigar to search for clipping
    firstOperation, firstOperationLen = alignmentObj.cigartuples[0]
    lastOperation, lastOperationLen = alignmentObj.cigartuples[-1]
    ## Clipping >= X bp at the left
    #  Note: soft (Operation=4) or hard clipped (Operation=5)     
    if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation == 4) or (lastOperation == 5)):
        SMSRead = True
    if ((lastOperation == 4) or (lastOperation == 5)) and ((firstOperation != 4) and (firstOperation != 5)):
        SMSRead = True

    return SMSRead

## TODO: WHERE CAN WE PUT THIS FUNCTION??
def percentage(counts, total):
    
    ## Percentage of SMS reads
    if counts > 0:
        perc = counts/total
    else:
        perc = 0
    
    return perc

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