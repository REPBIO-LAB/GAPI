'''
Module 'structures' - Contains functions organize data into more complex data structures
'''

## DEPENDENCIES ##
# External
# Internal
import log


## FUNCTIONS ##
def buildPosDict(event_list, windowSize):
    '''
    '''

    ## 1. Organize the objects into the dictionary 
    posDict = {}

    # For each event
    for event in event_list:
        
        # Determine to what window index the event belongs
        windowIndex = int(event.pos / windowSize)

        # a) First event into that window -> Initialize window
        if windowIndex not in posDict:
            posDict[windowIndex] = [ event ]
            
        # b) Add event into the existing window 
        else:
            posDict[windowIndex].append(event)

    ## 2. For each window sort the events in increasing coordinates order
    for windowIndex, events in posDict.items():

        events.sort(key=lambda event: event.pos)
        
    return posDict


