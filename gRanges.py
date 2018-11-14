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






    