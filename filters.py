'''
Module 'filters' - Contains functions for filtering clusters
'''

## External
import pysam

## Internal
import gRanges

###############
## FUNCTIONS ##
###############

def filter_metaclusters(metaclustersDict, filters2Apply, confDict):
    '''
    Function to apply filters all metaclusters. 

    Input:
        1. metaclustersDict: dictionary with the following structure: keys -> SV_type, value -> list of metaclusters corresponding to this SV_type.
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict

    Output:
        1. metaclustersPassDict: Dictionary with same structure as the input one, containing those metaclusters that passed all the filters.
        2. metaclustersFailDict: Dictionary with same structure as the input one, containig those metaclusters that failed one or more filters.
    '''

    metaclustersPassDict = {}
    metaclustersFailDict = {}

    ## For each type of SV
    for SV_type, metaclusters in metaclustersDict.items():

        ## 1. Make list with the indexes of the metaclusters do not passing some filter
        filteredIndexes = []

        ## For each metacluster
        for index, metacluster in enumerate(metaclusters):

            ## Apply filters
            metacluster.failedFilters = filter_metacluster(metacluster, filters2Apply, confDict)

            # Metacluster fails some filter
            if metacluster.failedFilters:
                filteredIndexes.append(index)

        ## 2. Divide metaclusters in those passing and failing filtering
        for index, metacluster in enumerate(metaclusters):
            
            ## a) Failing some filter
            if index in filteredIndexes:

                ## Initialize list
                if SV_type not in metaclustersFailDict:
                    metaclustersFailDict[SV_type] = []
                
                metaclustersFailDict[SV_type].append(metacluster)

            ## b) Passing all the filters
            else:

                ## Initialize list
                if SV_type not in metaclustersPassDict:
                    metaclustersPassDict[SV_type] = []

                metaclustersPassDict[SV_type].append(metacluster)

    return metaclustersPassDict, metaclustersFailDict

def filter_metacluster(metacluster, filters2Apply, confDict):
    '''
    Apply selected filters to one metacluster.

    Input:
        1. metacluster: metacluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict

    Output:
        1. failedFilters -> list containing those filters that the metacluster doesn't pass.
    '''        
    failedFilters = []

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: 
        if not minNbEventsFilter(metacluster, confDict['minSupportingReads'], confDict['minNormalSupportingReads']):
            failedFilters.append('MIN-NBREADS')

    ## 2. FILTER 2: Maximum number of reads per cluster
    if 'MAX-NBREADS' in filters2Apply: 
        if not maxNbEventsFilter(metacluster, confDict['maxClusterSize']):
            failedFilters.append('MAX-NBREADS')

    ## 3. FILTER 3: Maximum Coefficient of Variance per cluster
    if ('CV' in filters2Apply) and ('INS' in metacluster.subclusters): 
        if not maxCvFilter(metacluster, confDict['maxClusterCV']):
            failedFilters.append('CV')

    ## 4. FILTER 4: Whether a metacluster has a SV_type assigned or not
    if 'SV-TYPE' in filters2Apply: 
        if not SVTypeFilter(metacluster):
            failedFilters.append('SV-TYPE')

    ## 5. FILTER 5: Minimum percentage of inserted sequence resolved
    if ('PERC-RESOLVED' in filters2Apply) and (metacluster.SV_type == 'INS') and ('PERC_RESOLVED' in metacluster.SV_features): 

        if not percResolvedFilter(metacluster, confDict['minPercResolved']):
            failedFilters.append('PERC-RESOLVED')

    return failedFilters

def minNbEventsFilter(metacluster, minSupportingReads, minNormalSupportingReads):
    '''
    Filter metacluster by comparing the number of supporting events with a minimum treshold

    Input:
        1. metacluster: metacluster object
        2. minSupportingReads: min number of events threshold
        3. minSupportingReads: min number of events threshold for normal sample

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute number of events supporting the cluster 
    nbTotal, nbTumour, nbNormal, nbINS, nbDEL, nbCLIPPING = metacluster.nbEvents()

    ## 2. Compare the number of events supporting the cluster against the minimum required
    # 2.1 Paired mode:
    if nbTumour != None:

        if nbTumour >= minSupportingReads and nbNormal >= minNormalSupportingReads:
            metacluster.mutOrigin = 'germline'
            PASS = True

        elif nbTumour >= minSupportingReads and not nbNormal >= minNormalSupportingReads:
            metacluster.mutOrigin = 'somatic-tumour'
            PASS = True

        elif not nbTumour >= minSupportingReads and nbNormal >= minNormalSupportingReads:
            metacluster.mutOrigin = 'somatic-normal'
            PASS = True

        else:
            PASS = False

    # 2.1 Single mode:
    else:
        # If running in single mode (do no set mutation origin because it sdoesn't make sense)
        if nbTotal >= minSupportingReads:
            PASS = True
        else:
            PASS = False
    
    return PASS

def maxNbEventsFilter(metacluster, maxNbEvents):
    '''
    Filter metacluster by comparing the number of cluster supporting events with a maximum treshold

    Input:
        1. metacluster: metacluster object
        2. maxNbEvents: maximum number of events threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute number of events supporting the cluster 
    nbTotal = metacluster.nbEvents()[0]

    ## 2. Compare the number of events supporting the cluster against the maximum required
    if nbTotal <= maxNbEvents:
        PASS = True
    else:
        PASS = False
    
    return PASS

def maxCvFilter(metacluster, maxClusterCV):
    '''
    Filter metacluster by comparing its Coefficient of Variation with a maximum threshold.

    Input:
        1. metacluster: metacluster object
        2. maxClusterCV: maximum Coefficient of Variation threshold
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute CV of the cluster 
    cv = metacluster.subclusters['INS'].cv_len()[1]

    ## 2. Compare the cluster CV against the maximum required
    if cv <= maxClusterCV:
        PASS = True
    else:
        PASS = False

    return PASS

def percResolvedFilter(metacluster, minPercResolved):
    '''
    Filter metacluster by comparing the % of inserted sequence resolved with a minimum threshold

    Input:
        1. metacluster: metacluster object
        2. minPercResolved: minimum % resolved

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    if (metacluster.SV_features['PERC_RESOLVED'] >= minPercResolved):
        PASS = True

    else:
        PASS = False

    return PASS

def SVTypeFilter(metacluster):
    '''
    Filter metacluster by checking if it has a SV type assigned.

    Input:
        1. metacluster: metacluster object
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 2. Compare the percentage of outliers against the maximum required
    if metacluster.SV_type != None:
        PASS = True 
    else:
        PASS = False

    return PASS


def filter_discordant_mate_position(discordants, ref, beg, end):
    '''
    Filter out discordant cluster if mates align within provided region as input

    Input:
        1. discordants: list of discordant clusters (clustering of discordant done by proximity and then by mate position)
        2. ref: reference for excluded region
        3. beg: begin for excluded region 
        4. end: end for excluded region

    Output:
        1. filteredDiscordant: list of filtered discordant clusters
    '''
    filteredDiscordant = []

    for cluster in discordants:

        ## Retrieve mates interval of mate positions
        matesBeg, matesEnd = cluster.mates_start_interval()

        ## Check if mates interval overlaps with excluded region
        overlap, overlapLen = gRanges.overlap(matesBeg, matesEnd, beg, end)

        if not overlap:
            filteredDiscordant.append(cluster)
    
    return filteredDiscordant

# --------------- SHORT READS -----------------------
# HACER OTRA PARECIDA A LA QUE ESTABA PARA SHORT READS

## [SR CHANGE]
def applyFilters(clusters):
    '''
    Remove those clusters that fail in one or more filters
    '''
    newClusterDict = {}
    clusterTypes = clusters.eventTypes

    ## TODO: aqui seria mejor qitarlo de la lista en vez de hacer una nueva
    ## TODO: aqui mirar de coger el clusterType de otra manera sin tener que hacer dos loops:
    for clusterType in clusterTypes:
        for cluster in clusters.collect([clusterType]):
            newClusterDict[clusterType]=[]
            if False not in [value for value in cluster.filters.values()]:
                newClusterDict[clusterType].append(cluster)

    return newClusterDict


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
        
        if alignmentObj.cigartuples != None:
        
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
    percMAPQ = fraction(lowMAPQ, nbReads)
    percSMSReads = fraction(SMSReads, nbReads)

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
def fraction(counts, total):
    
    ## Percentage of SMS reads
    if counts > 0:
        perc = counts/total
    else:
        perc = 0
    
    return perc