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
        binDbObj.addEvents2bins(events, eventType)
    
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

    
    def addEvents2bins(self, events, eventType):
        '''
        Organize events into a hierarchy of genomic bins

        Input:
            1. events: list of objects. Every object must have 'beg' and 'end' attributes
            2. eventType: type of events (DEL, INS, CLIPPING, ...)

        Output:
            1. Add events to 'data' attribute. 'data' is a dictionary containing events organized in genomic bins:
            binIndex -> List of events belonging to that bin sorted by "mapAttribute"
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

                    # a) First event into that bin -> Initialize bin dict 
                    if binIndexBeg not in self.data[binSize]:
                        self.data[binSize][binIndexBeg] = {}
                        self.data[binSize][binIndexBeg][eventType] = [event]
                
                    # b) First event of that type in this bin -> Initialize events list
                    elif eventType not in self.data[binSize][binIndexBeg]:
                        self.data[binSize][binIndexBeg][eventType] = [event]
            
                    # c) There are already events of this type in this bin -> Add event to the list
                    else:
                        self.data[binSize][binIndexBeg][eventType].append(event)
        
                    # Do not check other bin sizes once event allocated in a bin
                    break

                ## B) Event spans several bins. Try with the next bin size 

        ## 2. For each bin sort the events in increasing coordinates order
        for binSize in self.data.keys():
            for binIndex in self.data[binSize].keys():
                if eventType in self.data[binSize][binIndex]:
                    self.data[binSize][binIndex][eventType].sort(key=lambda event: event.beg)

        
    def nbEvents(self):
        '''
        Compute the number of events composing the hash structure
        '''
        nbEvents = 0

        for binSize in self.data.keys():
            for binIndex in self.data[binSize].keys():
                for eventType in self.data[binSize][binIndex].keys():
                    events = self.data[binSize][binIndex][eventType]
                    nbEvents += len(events)     

        return nbEvents

    '''
    def binBoundaries(self, binIndex):
        Compute the begin and end boundaries for a given bin index
        beg = binIndex * self.binSize
        end = ((binIndex + 1) * self.binSize) - 1

    #    return beg, end
    '''


