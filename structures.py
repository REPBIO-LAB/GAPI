'''
Module 'structures' - Contains functions and classes to organize data into more complex data structures
'''

## DEPENDENCIES ##
# External
# Internal
import log

## FUNCTIONS ##
def create_bin_database(ref, beg, end, eventsDict, binSizes):
    '''

    Input:
        1. ref: reference/chromosome
        2. beg: bin begin coordinate
        3. end: bin end coordinate
        4. eventsDict: dictionary containing:
            * KEY_1 -> list of objects
            * KEY_2 -> list of objects
    
        5. binSizes: list of bin sizes that will be used to create a bin database 

    Output:
        1. binDb: 'bin_database' instance containing all the input events organized in genomic bins
    '''          

    # Initiate bin database
    binDb = bin_database(ref, beg, end, binSizes)

    # For each type of input event
    for eventType, events in eventsDict.items():

        # Add all the events from the given event type to the bin database
        binDb.addEvents(events, eventType)

    return binDb

## CLASSES ##
class bin_database():
    '''
    Database to organize a set of events into a hierarchy of genomic bins
    '''
    def __init__(self, ref, beg, end, binSizes):
        '''
        Organize events into a hierarchy of genomic bins

        Input:
            1. binSizes: list of bin sizes 
        '''
        self.ref = ref
        self.beg = int(beg)  
        self.end = int(end)
        self.binSizes = sorted(binSizes)
        self.data = {}
        
        ## Initialize one bin dictionary per size
        for binSize in self.binSizes:
            self.data[binSize] = {}

    def addEvents(self, events, eventType):
        '''
        Add events into a hierarchy of genomic bins

        Input:
            1. events: list of objects. Every object must have 'beg' and 'end' attributes
            2. eventType: type of events (DEL, INS, CLIPPING, ...)
        '''        
        ## 1. Allocate each event into a genomic bin 
        # For each event
        for event in events:
        
            # For each bin size (from smaller to bigger bin sizes)
            for binSize in self.binSizes:

                # Determine to what bin index the event belongs
                binIndexBeg = int(event.beg / binSize)
                binIndexEnd = int(event.end / binSize)

                # A) Event fits in one bin
                if (binIndexBeg == binIndexEnd):

                    # a) First event into that bin -> Initialize dict and bin  
                    if binIndexBeg not in self.data[binSize]:
                        self.data[binSize][binIndexBeg] = {}
                        self.data[binSize][binIndexBeg][eventType] = events_bin([event])                

                    # b) First event of that type in this bin -> Initialize bin
                    elif eventType not in self.data[binSize][binIndexBeg]:
                        self.data[binSize][binIndexBeg][eventType] = events_bin([event]) 

                    # c) There are already events of this type in this bin -> Add event to the bin
                    else:
                        self.data[binSize][binIndexBeg][eventType].add([event])
       
                    # Do not check other bin sizes once event allocated in a bin
                    break

                ## B) Event spans several bins. Try with the next bin size 

        ## 2. For each bin sort the events in increasing coordinates order
        for binSize in self.data.keys():
            for binIndex in self.data[binSize].keys():
                if eventType in self.data[binSize][binIndex]:
                    self.data[binSize][binIndex][eventType].sort()
        
    def collect(self, eventTypes):
        '''
        Collect all the events of target event types that are stored 
        in the bin database structure
        
         Input:
            1. eventTypes: list containing target event types

         Output:
            2. events. List of events
        '''  
        events = []

        # For each bin size
        for binSize in self.binSizes:

            # For each bin
            for binIndex in self.data[binSize].keys():
                
                # For each target event type
                for eventType in eventTypes:

                    # There are events of the target event type in the bin 
                    if eventType in self.data[binSize][binIndex]:

                        # Add events to the list
                        events = events + self.data[binSize][binIndex][eventType].events

        ## Sort events by begin coordinate
        events.sort(key=lambda event: event.beg)

        return events

    def collectEventTypes(self):

        eventTypes = []

        for binSize in self.binSizes:
            for dicti in self.data[binSize].values():
                for eventType in dicti.keys():
                    eventTypes.append(eventType)

        return eventTypes

    def collect_bin(self, binSize, binIndex, eventTypes):
        '''
        Collect all the events of target event types that are stored 
        in a particular bin 

         Input:
            1. binSize: bin size corresponding to the target bin index
            2. binIndex: target bin index
            3. eventTypes: list containing target event types

         Output:
            2. events. List of events
        '''  
        events = []

        ## Check if bin database contains target bin
        if (binSize in self.data) and (binIndex in self.data[binSize]):
    
            # For each target event type
            for eventType in eventTypes:

                # There are events of the target event type in the bin 
                if eventType in self.data[binSize][binIndex]:

                    # Add events to the list
                    events = events + self.data[binSize][binIndex][eventType].events

        ## Sort events by begin coordinate
        events.sort(key=lambda event: event.beg)

        return events


    def traverse(self, rootIndex, rootSize, eventType):
        '''
        Traverse bin structure starting in a root bin and going through all the bins located at upper levels
        in the hierarchy. Collect events from all the visited bins. E.g:

        <-----------------1---------------->
        <-------2--------><-------3-------->
        <---4---><---5---><---6---><---7---> 
        # BinIndex=4; Output: events in bins (4, 2, 1)
        # BinIndex=7; Output: events in bins (7, 3, 1)
        # BinIndex=2; Output: events in bins (2, 1)

        Input:
            1. rootIndex: root bin index
            2. rootSize: window size/level where the root index is located 
            3. eventType: type of events (DEL, INS, CLIPPING, ...)

        Output:
            1. events: list of events 
        '''      
        ### Initialize events list adding the events from the root bin
        rootBin = self.data[rootSize][rootIndex][eventType]
        events = rootBin.events

        ### Select upper windows sizes/levels
        upperSizes = [ binSize for binSize in self.binSizes if binSize > rootSize]

        # For each upper window size
        for upperSize in upperSizes:
                
            # Compute corresponding upper bin index
            upperIndex = int(rootIndex * rootSize / upperSize)

            # The upper bin contains events:
            if (upperIndex in self.data[upperSize]):

                # There are events of the type of interest on the bin 
                if (eventType in self.data[upperSize][upperIndex]):

                    ## Add upper bin´s events to the list
                    upperBin = self.data[upperSize][upperIndex][eventType]
                    events = events + upperBin.events

        return events

    def nbEvents(self):
        '''
        Compute the number of events composing the bin dictionary structure
        '''
        totalNbEvents = 0
        nbEventsBinSizes = {}

        for binSize in self.binSizes:
            nbEventsBinSizes[binSize] = 0

            for binIndex in self.data[binSize].keys():
                for eventType in self.data[binSize][binIndex].keys():
                    totalNbEvents += self.data[binSize][binIndex][eventType].nbEvents()
                    nbEventsBinSizes[binSize] += self.data[binSize][binIndex][eventType].nbEvents()

        return totalNbEvents, nbEventsBinSizes


class events_bin():
    '''
    Contain a set of events of the same type in a genomic bin
    '''
    def __init__(self, events):
        self.events = events

    def add(self, events):
        '''
        Add input events to the bin´s list of events

        Input:
            1. events: List of events to add to the bin 
        '''
        self.events = self.events + events
        
    def nbEvents(self):
        '''
        Compute the number of events composing the bin
        '''
        return len(self.events)

    def sort(self):
        '''
        Sort events in increasing coordinates order
        '''
        self.events.sort(key=lambda event: event.beg)