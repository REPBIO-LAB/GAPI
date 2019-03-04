'''
Module 'gRanges' - Contains functions to do operations with genomic ranger or coordinates (i.e. overlaps...)
'''

## DEPENDENCIES ##

## FUNCTIONS ##
def overlap(begA, endA, begB, endB):
    '''
    Check if two ranges overlap and return the number of overlapping bases. 
    Coordinate ranges are expected to be 1-based
    '''    
    maxBeg = max([begA, begB])
    minEnd = min([endA, endB])

    # a) Overlap
    if maxBeg <= minEnd:
        overlapLen = minEnd - maxBeg + 1

    # d) No overlap
    else:
        overlapLen = 0

    return overlapLen
        
def rcplOverlap(begA, endA, begB, endB, minPercOverlap):
    '''
    Check if two ranges overlap with a minimum percentage of reciprocal overlap. 
    '''    
    overlapLen = overlap(begA, endA, begB, endB)
    lenA = (endA - begA)
    lenB = (endB - begB)
    percA = (overlapLen / lenA) * 100
    percB = (overlapLen / lenB) * 100

    # a) Both ranges overlap with a minimum X percentage of reciprocal overlap
    if (percA >= minPercOverlap) and (percB >= minPercOverlap):
        boolean = True

    # b) No overlap
    else:
        boolean = False

    return boolean


def complementary(begA, endA, begB, endB, maxDist, maxPercOverlap):
    '''
    Check if two intervals are complementary. Two ranges are considered complementary when:
    1) Ranges located less than X bp of distance AND
    2) Overlapping region does not span X% or more of any of the intervals

    Input:
        1. begA: begin position of range A
        2. endA: end position of range A
        3. begB: begin position of range B
        4. endB: end position of range B
        5. maxDist: maximum distance between both ranges 
        6. maxPercOverlap: maximum percentage of overlap between ranges

    Output:
        1. boolean: True (intervals are complementary) or False (intervals not complementary)
        2. orientation: 'LEFT' (B on the left of A), 'RIGHT' (B on the right of A) or None (not complementary)
    '''  
    ## 1. Redefine intervals by adding the maximum distance to their begin and end coordinates
    newBegA = begA - maxDist
    newEndA = endA + maxDist
    newBegB = begB - maxDist
    newEndB = endB + maxDist

    ## 2. Assess if redefined intervals do overlap
    nbBases = overlap(newBegA, newEndA, newBegB, newEndB)
    
    # A) No overlap
    if nbBases == 0:
        boolean = False
        orientation = None
    
    # B) Overlap 
    else:

        # 3. Discard those cases with an overlapping region spanning X% or more of at least one of the intervals
        # --------A---------
        #              <---> Overlapping region
        #              -------B-------
        nbBases = overlap(begA, endA, begB, endB)   # Compute real degree of overlap

        lenA = (endA - begA)
        lenB = (endB - begB)
        percA = (nbBases / lenA) * 100
        percB = (nbBases / lenB) * 100

        # a) One of the intervals overlaps the other by more than the maximum allowed % cutoff 
        if (percA > maxPercOverlap) or (percB > maxPercOverlap):
            boolean = False
            orientation = None

        # b) Overlapping region does not span X% or more of any of the intervals
        else:
            boolean = True

            ## Determine complementariety orientation (use A as reference)
            # a) B located on the left of A
            #                   begA <------A------> 
            #  begB <------B------>
            if begA > begB:
                orientation = 'LEFT'

            # b) B located on the right of A
            #  begA <------A------>
            #                  begB <------B------> 
            else:
                orientation = 'RIGHT'

    return boolean, orientation