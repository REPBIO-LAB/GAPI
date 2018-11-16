'''
Module 'structures' - Contains functions and classes to organize data into more complex data structures
'''

## DEPENDENCIES ##
# External
# Internal
import log

## FUNCTIONS ##
def createBinDb(data, binSizes):
    '''
    Input:
        1. data: tuple list containing the list of events corresponding to each provided event type. E.g: 
                 [ tuple1(event_list, event_type), ..., tupleN(event_list, event_type)]
        2. binSizes: list of bin sizes 

    Output:
        1. binDbObj: 'binDb' instance containing all the input events organized in genomic bins
    '''            
    binDbObj = binDb(binSizes)
    
    # For each type of event add the corresponding 
    # event instances to the data structure 
    for events, eventType in data:
        binDbObj.addEvents(events, eventType)
    
    return binDbObj


## CLASSES ##
class binDb():
    '''
    Database to organize a set of events into a hierarchy of genomic bins
    '''
    def __init__(self, binSizes):
        '''
        Organize events into a hierarchy of genomic bins

        Input:
            1. binSizes: list of bin sizes 
        '''
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
                        self.data[binSize][binIndexBeg][eventType] = eventsBin([event])                

                    # b) First event of that type in this bin -> Initialize bin
                    elif eventType not in self.data[binSize][binIndexBeg]:
                        self.data[binSize][binIndexBeg][eventType] = eventsBin([event]) 

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
        
    def collect(self, eventType):
        '''
        Collect all the events from the provided event type that are stored in the bin structure
        
         Input:
            1. eventType. Target event type

         Output:
            2. events. List of events
        '''  
        events = []

        for binSize in self.binSizes:

            for binIndex in self.data[binSize].keys():
                
                events = events + self.data[binSize][binIndex][eventType].events

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
        Compute the number of events composing the hash structure
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


class eventsBin():
    '''
    Contain a set of events of the same type in a genomic bin
    '''
    def __init__(self, events):
        self.events = events

    def add(self, events):
        '''
        Contain a set of events of the same type in a genomic bin
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