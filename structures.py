'''
Module 'structures' - Contains functions organize data into more complex data structures
'''

## DEPENDENCIES ##
# External
# Internal
import log


## FUNCTIONS ##
def mapEvents2windows(events, mapAttribute, windowSize):
    '''
    Organize events into windows of coordinates. 

    Input:
        1. events: list of objects. Every object must have a "mapAttribute" argument 
        2. mapAttribute: Object attribute used for mapping the objects to the corresponding window
        3. windowSize: Size of each window

    Output:
        1. windowDict: dictionary containing the clusters organized in genomic windows:
            windowIndex -> Objects list belonging to that window sorted by "mapAttribute"
    '''

    ## 1. Organize the objects into the dictionary 
    windowDict = {}

    # For each event
    for event in events:
        
        # Determine to what window index the event belongs
        windowIndex = int(getattr(event, mapAttribute)/ windowSize)

        # a) First event into that window -> Initialize window
        if windowIndex not in windowDict:
            windowDict[windowIndex] = [ event ]
            
        # b) Add event into the existing window 
        else:
            windowDict[windowIndex].append(event)

    ## 2. For each window sort the events in increasing coordinates order
    for windowIndex, events in windowDict.items():

        events.sort(key=lambda event: getattr(event, mapAttribute))
        
    return windowDict


