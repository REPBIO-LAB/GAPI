'''
Module 'structures' - Contains functions and classes to organize data into more complex data structures
'''

## DEPENDENCIES ##
# External
# Internal
import log


## FUNCTIONS ##

## CLASSES ##

class windowsHash():
    '''
    Hash structure to organize a set of events into genomic windows/intervals
    '''
    def __init__(self, events, mapAttribute, windowSize):
        '''
        Organize events into windows of coordinates. 

        Input:
            1. events: list of objects. Every object must have the attribute defined in "mapAttribute" variable
            2. mapAttribute: Object attribute used for mapping the objects to the corresponding window
            3. windowSize: Size of each window
        '''
        self.mapAttribute = mapAttribute
        self.windowSize = windowSize
        self.data = {}
        self.mapEvents2windows(events)
    
    def mapEvents2windows(self, events):
        '''
        Organize events into windows of coordinates. 

        Input:
            1. events: list of objects. Every object must have the attribute defined in "mapAttribute" variable

        Output:
            1. Add events to 'data' attribute. 'data' is a dictionary containing events organized in genomic windows:
            windowIndex -> List of events belonging to that window sorted by "mapAttribute"
        '''

        ## 1. Organize the objects into the dictionary 
        # For each event
        for event in events:
        
            # Determine to what window index the event belongs
            windowIndex = int(getattr(event, self.mapAttribute)/ self.windowSize)

            # a) First event into that window -> Initialize window
            if windowIndex not in self.data:
                self.data[windowIndex] = [ event ]
            
            # b) Add event into the existing window 
            else:
                self.data[windowIndex].append(event)

        ## 2. For each window sort the events in increasing coordinates order
        for windowIndex, events in self.data.items():

            events.sort(key=lambda event: getattr(event, self.mapAttribute))
        
    def nbEvents(self):
        '''
        Compute the number of events composing the hash structure
        '''
        nbEvents = 0

        for events in self.data.values():
            nbEvents += len(events)     

        return nbEvents



