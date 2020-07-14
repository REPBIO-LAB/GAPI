'''
Module 'filters' - Contains functions for filtering clusters
'''
## External
import pysam
import statistics

## Internal
import gRanges
import bamtools
import stats


###############
## FUNCTIONS ##
###############

def filter_clippings(clippings, filters2Apply, confDict):
    '''
    Function to apply filters all metaclusters. 

    Input:
        1. clippings: list of clipping clusters
        2. filters2Apply: list containing the filters to apply       
        3. confDict

    Output:
        1. clippingsPass: List of clippings clusters passing all the filters    
    '''
    clippingsPass = []

    # For discordant cluster
    for clipping in clippings:
        
        ## Apply filters
        failedFilters = filter_clipping(clipping, filters2Apply, confDict)

        # Metacluster pass all the filters
        if not failedFilters: 
            clippingsPass.append(clipping)
        
        else:
            print("Discarded clipping cluster")
            print(clipping.ref, clipping.beg, clipping.end)
            print(failedFilters)
            print([event.readName for event in clipping.events])

    return clippingsPass


def filter_clipping(clipping, filters2Apply, confDict):
    '''
    Apply selected filters to a clipping cluster provided as input

    Input:
        1. clipping: clipping cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        5. confDict:

    Output:
        1. failedFilters -> list containing those filters that the clipping cluster didn´t pass.
    '''        
    failedFilters = []

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: 
        if not filter_min_nb_reads(clipping, confDict['minNbCLIPPING'], confDict['minNormalReads']):
            failedFilters.append('MIN-NBREADS')
        
    ## 2. FILTER 2: Filter out those clusters with suppl outside target reference ##
    if 'SUPPL-REF' in filters2Apply: 

        if not filter_clipping_suppl_ref(clipping, confDict['targetRefs']):
            failedFilters.append('SUPPL-REF')

    ## 3. FILTER 3: Filter out those clusters whose supplementary alignments map over any source element downstream region 
    if 'SUPPL-SRC' in filters2Apply:

        if not filter_clipping_suppl_position(clipping, confDict['rangesDict'], 10000):
            failedFilters.append('SUPPL-SRC')
    
    ## 4. FILTER 4: Average mapping quality of supplementary alignments 
    if 'SUPPL-MAPQ' in filters2Apply:
    
        if not filter_suppl_MAPQ(clipping, 20):
            failedFilters.append('SUPPL-MAPQ')

    ## 5. FILTER 5: filter out clusters formed by tumour and normal reads. Discard germline variation
    if 'GERMLINE' in filters2Apply:

        if not filter_germline(clipping, confDict['minNormalReads']):
            failedFilters.append('GERMLINE')

    ## 6. FILTER 6: Filter out clusters based on duplicate percentage (Ex: 40%) 
    if 'READ-DUP' in filters2Apply:

        if not filter_highDup_clusters(clipping, 75):
            failedFilters.append('READ-DUP')
    
    ## 7. FILTER 7: Filter out clusters based on cluster coordinates ##
    if 'CLUSTER-RANGE' in filters2Apply:
        
        if not filter_clusterRange_clipping(clipping):
            failedFilters.append('CLUSTER-RANGE')

    return failedFilters


def filter_discordants(discordants, filters2Apply, bam, normalBam, confDict):
    '''
    Function to apply filters all discordant clusters. 

    Input:
        1. discordants: list of discordant clusters
        2. filters2Apply: list containing the filters to apply 
        3. bam: path to bam file. None if not needed
        4. normalBam: path to matched normal bam file. None if not needed        
        5. confDict

    Output:
        1. discordantsPass: List of discordant clusters passing all the filters
    '''
    discordantsPass = []
    
    # For discordant cluster
    for discordant in discordants:

        ## Apply filters
        failedFilters = filter_discordant(discordant, filters2Apply, bam, normalBam, confDict)

        # Cluster pass all the filters
        if not failedFilters: 
            discordantsPass.append(discordant)
            
        else:
            print("Discarded discordant cluster")
            print(discordant.ref, discordant.beg, discordant.end)
            print(failedFilters)
            print([event.readName for event in discordant.events])
    
    return discordantsPass

def filter_discordant(discordant, filters2Apply, bam, normalBam, confDict):
    '''
    Apply selected filters to a discordant cluster provided as input

    Input:
        1. discordant: discordant read pair cluster object
        2. filters2Apply: list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. bam: path to bam file. None if not needed
        4. normalBam: path to matched normal bam file. None if not needed
        5. confDict

    Output:
        1. failedFilters -> list containing those filters that the discordant cluster doesn't pass.
    '''        
    failedFilters = []

    ## 1. FILTER 1: Minimum number of reads per cluster
    if 'MIN-NBREADS' in filters2Apply: 
        if not filter_min_nb_reads(discordant, confDict['minNbDISCORDANT'], confDict['minNormalReads']):
            failedFilters.append('MIN-NBREADS')

    ## 2. FILTER 2: Filter out those clusters with mates not over NOT target reference ##
    if 'MATE-REF' in filters2Apply: 

        if not filter_discordant_mate_ref(discordant, confDict['targetRefs']):
            failedFilters.append('MATE-REF')

    ## 3. FILTER 3: Filter out those clusters whose mates aligns over any source element downstream region ##
    if 'MATE-SRC' in filters2Apply:

        if not filter_discordant_mate_position(discordant, confDict['rangesDict'], 10000):
            failedFilters.append('MATE-SRC')
        
    ## 4. FILTER 4: Filter out clusters based on average MAPQ for mate alignments ##
    if 'MATE-MAPQ' in filters2Apply:
    
        if not filter_discordant_mate_MAPQ(discordant, 20, bam, normalBam):
            failedFilters.append('MATE-MAPQ')
        
    ## 5. FILTER 5: filter out clusters formed by tumour and normal reads. Discard germline variation (TILL HERE) 
    if 'GERMLINE' in filters2Apply:

        if not filter_germline(discordant, confDict['minNormalReads']):
            failedFilters.append('GERMLINE')
            
    ## 6. FILTER 6: Filter out clus in unspecific regions ##
    if 'UNESPECIFIC' in filters2Apply:

        if not filter_discordant_mate_unespecific(discordant, 0.2, bam):
            failedFilters.append('UNESPECIFIC')

    ## 7. FILTER 7: Filter out clusters based on duplicate percentage (Ex: 50%) ##
    if 'READ-DUP' in filters2Apply:

        if not filter_highDup_clusters(discordant, 75):
            failedFilters.append('READ-DUP')
            
    ## 8. FILTER 8: Filter out clusters based on mates beg coordinates ##
    if 'CLUSTER-RANGE' in filters2Apply:
        
        if not filter_clusterRange_discordant(discordant):
            failedFilters.append('CLUSTER-RANGE')
    
    return failedFilters

def filter_metaclusters_SR(metaclusters, filters2Apply, confDict, bam):
    '''
    Function to apply filters all metaclusters. 

    Input:
        1. metaclustersDict: dictionary with the following structure: keys -> SV_type, value -> list of metaclusters corresponding to this SV_type.
        2. filters2Apply: dictionary containing lists as values list containing the filters to apply (only those filters that make sense with the cluster type will be applied)
        3. confDict
        4. bam

    Output:
        1. metaclustersPassDict: Dictionary with same structure as the input one, containing those metaclusters that passed all the filters.
        2. metaclustersFailDict: Dictionary with same structure as the input one, containig those metaclusters that failed one or more filters.
    '''

    metaclustersPass = []
    metaclustersFail = []

    ## For each type of SV
    #for metacluster in metaclusters:
    ## 1. Make list with the indexes of the metaclusters do not passing some filter
    filteredIndexes = []

    ## For each metacluster
    for index, metacluster in enumerate(metaclusters):

        # Set meacluster element:
        element = metacluster.setElement() if metacluster.setElement() else 'GENERIC'
        ## Apply filters
        metacluster.failedFilters = filter_metacluster(metacluster, filters2Apply[element], confDict, bam)

        # Metacluster fails some filter
        if metacluster.failedFilters:
            filteredIndexes.append(index)

    ## 2. Divide metaclusters in those passing and failing filtering
    for index, metacluster in enumerate(metaclusters):
        
        ## a) Failing some filter
        if index in filteredIndexes:

            metaclustersFail.append(metacluster)

        ## b) Passing all the filters
        else:

            metaclustersPass.append(metacluster)

    return metaclustersPass, metaclustersFail


def filter_metaclusters(metaclustersDict, filters2Apply, confDict):
    '''
    Function to apply filters to a set of metaclusters organized in a dictionary

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
            metacluster.failedFilters = filter_metacluster(metacluster, filters2Apply, confDict, None)

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

def filter_metacluster(metacluster, filters2Apply, confDict, bam):
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
        if not filter_min_nb_reads(metacluster, confDict['minReads'], confDict['minNormalReads']):
            failedFilters.append('MIN-NBREADS')

    ## 2. FILTER 2: Maximum number of reads per cluster
    if 'MAX-NBREADS' in filters2Apply: 
        if not filter_max_nb_reads(metacluster, confDict['maxClusterSize']):
            failedFilters.append('MAX-NBREADS')

    ## 3. FILTER 3: Maximum Coefficient of Variance per cluster
    if ('CV' in filters2Apply) and ('INS' in metacluster.subclusters): 
        if not filter_max_cv(metacluster, confDict['maxClusterCV']):
            failedFilters.append('CV')

    ## 4. FILTER 4: Whether a metacluster has a SV_type assigned or not
    if 'SV-TYPE' in filters2Apply: 
        if not filter_SV_type(metacluster, confDict['targetSV']):
            failedFilters.append('SV-TYPE')

    ## 5. FILTER 5: Minimum percentage of inserted sequence resolved
    if ('PERC-RESOLVED' in filters2Apply) and (metacluster.SV_type == 'INS') and ('PERC_RESOLVED' in metacluster.SV_features): 

        if not filter_perc_resolved(metacluster, confDict['minPercResolved']):
            failedFilters.append('PERC-RESOLVED')

    # NOTE 2020: NEw 2020
    # TODO 2020: Put in another way!!
    '''
    ## 5. FILTER 5: Minimum percentage of inserted sequence resolved
    if ('PERC-RESOLVED' in filters2Apply) and ('PERC_RESOLVED' in metacluster.SV_features): 

        if not filter_perc_resolved(metacluster, confDict['minPercResolved']):
            failedFilters.append('PERC-RESOLVED')
    '''

    ## 3. FILTER 3: Area mapping quality
    if "AREAMAPQ" in filters2Apply:
        if not area(metacluster,confDict,bam)[0]:
            failedFilters.append('AREAMAPQ')

    ## 4. FILTER 4: Area clipping SMS
    if "AREASMS" in filters2Apply:
        if not area(metacluster,confDict,bam)[1]:
            failedFilters.append('AREASMS')

    ## 5. FILTER 5: Whether a metacluster has a SV_type assigned or not
    if 'IDENTITY' in filters2Apply: 
        if not identityFilter(metacluster):
            failedFilters.append('IDENTITY')

    return failedFilters

def filter_min_nb_reads(cluster, minReads, minNormalReads):
    '''
    Filter cluster by comparing the number of supporting events with a minimum treshold

    Input:
        1. cluster: cluster object
        2. minReads: min number of events threshold
        3. minReads: min number of events threshold for normal sample

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 1. Compute number of events supporting the cluster 
    nbTotal, nbTumour, nbNormal = cluster.supportingReads()[0:3]

    ## 2. Compare the number of events supporting the cluster against the minimum required
    # 2.1 Paired mode:
    if nbTumour != None:

        if nbTumour >= minReads and nbNormal >= minNormalReads:
            cluster.mutOrigin = 'germline'
            PASS = True

        elif nbTumour >= minReads and not nbNormal >= minNormalReads:
            cluster.mutOrigin = 'somatic-tumour'
            PASS = True

        elif nbNormal >= minReads:
            cluster.mutOrigin = 'somatic-normal'
            PASS = True

        else:
            PASS = False

    # 2.1 Single mode:
    else:
        # If running in single mode (do no set mutation origin because it sdoesn't make sense)
        if nbTotal >= minReads:
            PASS = True
        else:
            PASS = False
    
    return PASS

def filter_max_nb_reads(metacluster, maxNbEvents):
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

def filter_max_cv(metacluster, maxClusterCV):
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

def filter_perc_resolved(metacluster, minPercResolved):
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

def filter_SV_type(metacluster, targetSV):
    '''
    Filter metacluster by checking its SV type

    Input:
        1. metacluster: metacluster object
        2. targetSV: list containing list of target SV types
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    if metacluster.SV_type in targetSV:
        PASS = True 
        
    else:
        PASS = False

    return PASS

def identityFilter(cluster):
    '''
    Filter metacluster by checking its SV type

    Input:
        1. metacluster: metacluster object
        2. targetSV: list containing list of target SV types
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## 2. Compare the percentage of outliers against the maximum required
    if cluster.identity != None:
        PASS = True 
        
    else:
        PASS = False

    return PASS

def filter_suppl_MAPQ(clipping, minMAPQ):
    '''
    Filter out clipping cluster based on the average mapping quality of its supplementary alignments

    Input:
        1. clipping: clipping cluster object
        2. minMAPQ: minimum mapping quality

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't    
    '''

    ## Compute average mapping quality for supplementary alignments
    MAPQs = [event.mapQ for event in clipping.supplCluster.events]
    meanMAPQ = statistics.mean(MAPQs)

    ## Apply filter
    if meanMAPQ >= minMAPQ:
        PASS = True

    else:
        PASS = False
    
    return PASS

def filter_clipping_suppl_ref(clipping, targetRefs):
    '''
    Filter out clipping cluster if its supplementary cluster not located over a target reference

    Input:
        1. clipping: clipping cluster
        2. targetRefs: List of target references

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Retrieve suppl. alignment reference 
    if clipping.supplCluster.ref in targetRefs:
        PASS = True 

    else:
        PASS = False

    return PASS

def filter_discordant_mate_ref(discordant, targetRefs):
    '''
    Filter out discordant cluster if its mate not located over a target reference

    Input:
        1. discordant: discordant cluster object
        2. targetRefs: List of target references

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Retrieve mates referene 
    matesRef = discordant.events[0].mateRef

    if matesRef in targetRefs:
        PASS = True 
        
    else:
        PASS = False

    return PASS


def filter_clipping_suppl_position(clipping, ranges, buffer):
    '''
    Filter out clipping cluster if suppl. alignment within one of the provided regions
    
    Input:
        1. clipping: clipping cluster object
        2. ranges: Dictionary with reference ids as keys and the list of ranges on each reference as values
        3. buffer: Extend each range at their begin and end coordinate by a number of nucleotides == buffer length

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Do not filter out clipping cluster if no input range on that particular reference 
    if clipping.supplCluster.ref not in ranges:
        PASS = True
        return PASS
        
    ## Assess overlap between mates interval and provided regions. 
    PASS = True

    for interval in ranges[clipping.supplCluster.ref]:
        rangeBeg, rangeEnd = interval

        # Add buffer
        rangeBeg = rangeBeg - buffer
        rangeEnd = rangeEnd + buffer

        # Assess overlap
        overlap, overlapLen = gRanges.overlap(clipping.supplCluster.beg, clipping.supplCluster.end, rangeBeg, rangeEnd)

        if overlap:
            PASS = False

    return PASS

def filter_discordant_mate_position(discordant, ranges, buffer):
    '''
    Filter out discordant cluster if mates align within one of the provided regions

    Input:
        1. discordant: discordant cluster object
        2. ranges: Dictionary with reference ids as keys and the list of ranges on each reference as values
        3. buffer: Extend each range at their begin and end coordinate by a number of nucleotides == buffer length

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''

    ## Retrieve mates position and interval
    matesRef = discordant.events[0].mateRef
    matesBeg, matesEnd = discordant.mates_start_interval()

    ## Do not filter out discordant cluster if no input range on that particular reference 
    if matesRef not in ranges:
        PASS = True
        return PASS
        
    ## Assess overlap between mates interval and provided regions. 
    PASS = True

    for interval in ranges[matesRef]:
        rangeBeg, rangeEnd = interval

        # Add buffer
        rangeBeg = rangeBeg - buffer
        rangeEnd = rangeEnd + buffer

        # Assess overlap
        overlap, overlapLen = gRanges.overlap(matesBeg, matesEnd, rangeBeg, rangeEnd)

        if overlap:
            PASS = False

    return PASS

def filter_discordant_mate_MAPQ(discordant, minMAPQ, bam, normalBam):
    '''
    Filter out discordant clusters based on average MAPQ for mate alignments

    Input:
        1. discordant: discordant cluster object
        2. minMAPQ: minimum average of mapping quality for mate alignments
        3. bam: path to bam file containing alignments for discordant cluster supporting reads and their mate
        4. normalBam: path to the matched normal bam file. If running in single mode, set to 'None' 

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    ## Open BAM files for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    if normalBam != None:
        bamFile_normal = pysam.AlignmentFile(normalBam, "rb")
        
    ## Define interval to search for mate alignment objects
    matesBeg, matesEnd = discordant.mates_start_interval()

    intervalRef = discordant.events[0].mateRef 
    intervalBeg = matesBeg - 500
    intervalEnd = matesEnd + 500

    ## Collect cluster supporting reads
    nbTotal, nbTumour, nbNormal, readIds, readIdsTumour, readIdsNormal = discordant.supportingReads()

    ## Compute average mapping quality for mates of cluster supporting reads
    # if running in single mode or running in paired but no reads from the normal are found in this interval
    if (normalBam == None or (normalBam != None and readIdsNormal == [])):
            
        avMAPQ = bamtools.average_MAPQ_reads_interval(intervalRef, intervalBeg, intervalEnd, readIds, bamFile)
            
        if avMAPQ >= minMAPQ:
            PASS = True

        else:
            PASS = False
                
    # if running in paired mode and there is reads in the interval belonging to the normal bam
    else:
            
        avMAPQ = bamtools.average_MAPQ_reads_interval(intervalRef, intervalBeg, intervalEnd, readIdsTumour, bamFile)
        avMAPQ_normal = bamtools.average_MAPQ_reads_interval(intervalRef, intervalBeg, intervalEnd, readIdsNormal, bamFile_normal)
            
        # tumor and normal MAPQ average
        avMAPQ_pair = (avMAPQ * len(readIdsTumour) + avMAPQ_normal * len(readIdsNormal))/len(readIdsTumour + readIdsNormal)          
             
        if avMAPQ_pair >= minMAPQ:
            PASS = True

        else:
            PASS = False

    ## Close 
    bamFile.close()

    if normalBam != None:
        bamFile_normal.close()
        
    return PASS


def filter_germline(cluster, minNormal):
    '''
    Filter out those clusters formed by tumour and normal reads
    
    Input:
        1. cluster: cluster object
        2. minNormal: minimum number of reads supporting a SV in normal sample
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    count = 0
        
    for event in cluster.events:
            
        if event.sample == 'NORMAL':
            count += 1
                
    if count < minNormal:
        PASS = True
        
    else:
        PASS = False

    return PASS


def filter_discordant_mate_unespecific(discordant, threshold, bam):
    '''
    Filter out discordant whose mates are located in regions captured unspecifically
    Insertion points where there is more than discordant reads are likely to be false positives
    Example of recurrent false positive filtered: chr6:29765954 (hg19)
    
    Input:
        1. discordant: discordant cluster instance
        2. threshold: ratio of nbDiscordant reads between nbTotal reads
        3. bam: path to bam file

    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    # Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")
    
    nbProperPair = 0
    nbDiscordants = len(discordant.events)
        
    buffer = 200
    ref = discordant.events[0].mateRef
    beg, end = discordant.mates_start_interval()
        
    # Extract alignments
    iterator = bamFile.fetch(ref, beg - buffer, end + buffer)
        
    # Count properly paired reads in region
    for alignmentObj in iterator:
            
        if alignmentObj.is_proper_pair:
                
            nbProperPair += 1
        
    # if there are properly paired reads around insertion
    if nbProperPair > 0:
            
        # if the ratio discordants/properly paired reads is greater than threshold
        if nbDiscordants/nbProperPair > threshold:        
            PASS = True

        else:  
            PASS = False

    else:
            
        PASS = True
            
    return PASS

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
        filterDiscordantResults['MAX-NBREADS'] = filter_max_nb_reads(cluster, confDict['maxClusterSize'])

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
    percMAPQ = stats.fraction(lowMAPQ, nbReads)
    percSMSReads = stats.fraction(SMSReads, nbReads)

    ## If the percentage of low MQ reads is lower than the threshold pass the filter.
    if percMAPQ != None:
        if percMAPQ < float(maxRegionlowMQ):
            percMAPQFilter = True
        else:
            percMAPQFilter = False
    else:
        percMAPQFilter = False

    ## If the percentage of SMS reads is lower than the threshold pass the filter.
    if percSMSReads != None:
        if percSMSReads < float(maxRegionSMS):
            percSMSReadsFilter = True
        else:
            percSMSReadsFilter = False
    else:
        percSMSReadsFilter = False
    
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
        
    #if ((lastOperation == 4) or (lastOperation == 5)) and ((firstOperation != 4) and (firstOperation != 5)):
        #SMSRead = True

    return SMSRead


def filter_highDup_clusters(cluster, maxDupPerc):
    '''
    Filter out those clusters formed by more than a percentage of duplicates
    
    Input:
        1. cluster: list of discordant clusters formed by DISCORDANT events
        2. maxDupPerc: Maximun % of duplicates allowed (not included)
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    dupPerc = cluster.dupPercentage()
        
    if dupPerc < maxDupPerc:
            
        PASS = True

    else:
        PASS = False
                
    return PASS



def filter_clusterRange_discordant(cluster):
    '''
    Filter out those discordant clusters in which all reads are piled up.
    This filter is only applied to clusters formed by more than a discordant alignment
    
    ----------****>                       ---------******>
       -------*******>                    ---------******>
     ---------*****>                      ---------******>
         -----*********>                  ---------******>
    |         |                          |         |   
    beg      end                         beg      end
    ----------       clusterRange         ----------
         -----      min(alignRanges)      ----------
       True             PASS                 False      
    
    Input:
        1. cluster: cluster formed by DISCORDANT events
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    PASS = True
        
    # if there is more than a discordant alignment
    if len(cluster.events) > 1:
               
        # define cluster range
        clusterRange = cluster.end - cluster.beg
        
        # define minimum alignment range of all reads supporting the cluster
        readRanges = []
        
        for event in cluster.events:
            beg, end = event.readCoordinates()
            readRange = abs(abs(end) - abs(beg))
            readRanges.append(readRange)
        
        alignRange = min(readRanges)
        
        # if the cluster range is smaller or equal to the minimum alignment range
        if (clusterRange <= alignRange):
            
            # discard the cluster
            PASS = False
    
    return PASS



def filter_clusterRange_clipping(cluster):
    '''
    Filter out those clipping clusters in which all clippings have the same 
    coordinates relative to the read. This filter is only applied to clusters 
    formed by more than a clipping alignment
    
    ----------****>                       ---------******>
             55   75                              40    75
       -------*******>                    ---------******>
              38    75                            40    75
     ---------*****>                      ---------******>
              61   75                             40    75
                          clusterCoord 
    [55, 75, 38, 75, 61, 75]        [40, 75, 40, 75, 40, 75]
          True              PASS                False      
    
    Input:
        1. cluster: cluster formed by CLIPPING events
    
    Output:
        1. PASS -> boolean: True if the cluster pass the filter, False if it doesn't
    '''
    
    PASS = True
       
    # if there is more than a clipping
    if len(cluster.events) > 1:
               
        # with clipping events beg and end coordinates are the same. However, readCoordinates()
        # returns the aligment coordinates relative to the read. If more than 2 coordinates, 
        # it is not a piled up cluster of clippings        
        clusterCoord = [event.readCoordinates() for event in cluster.events]
        clusterCoord_flatList = [item for sublist in clusterCoord for item in sublist]
        n_clusterCoord = len(set(clusterCoord_flatList))
              
        if(n_clusterCoord <= 2):
            PASS = False
        
    return PASS
