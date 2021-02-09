'''
Module 'singleCell' - Contains classes for with single cell data (SVs, for now...)
'''

# External
import multiprocessing as mp

# Internal 
import depth 

class population():
    '''
    Class representing a population of single cells
    '''
    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.cells = {}

    def create_cells(self, cellIds, refLen):
        '''
        Create set of cells 

        Input:
            1. cellIds: list of cell identifiers
            2. refLen: dictionary with a key per reference/chromosome and its length as value

        Output:
            Update 'cells' argument with a dictionary containing cellIds as keys and cell instances as values
        '''
        ## For each cell
        for cellId in cellIds:
            
            ## Initialize cell
            cellInstance = cell(cellId, refLen)

            ## Add to the dict
            self.cells[cellId] = cellInstance

    def nbCells(self):
        '''
        Return the number of cells in the population
        '''
        return len(self.cells)

    def cells2list(self):
        '''
        Return list with all the cells in the population
        '''
        return list(self.cells.values())

    def bam2cells(self, bamDict):
        '''
        Add bam file to each corresponding cell

        Input:
            1. bamDict: dictionary with cellIds as keys and bam paths as values

        Output:
            Update 'bam' attribute for each cell in the population with the path to their corresponding bam files
        '''
        # Iterave over the bam files
        for cellId, bamPath in bamDict.items():

            # Add bam path
            self.cells[cellId].bam = bamPath
    
    def read_depth(self, binSize, minMAPQ, filterDup, targetRefs, processes):
        '''
        Compute read depth genome wide for all cells in the population

        Input:
            1. binSize: segments size 
            2. minMAPQ: minimum read mapping quality
            3. filterDup: filter read duplicates (True) or not (False)
            4. targetRefs: list with target references
            5. processes: number of processes for parallelization
        '''
        ## 1. Make cell list
        cells = self.cells2list()
        cells = [[cell, binSize, minMAPQ, filterDup, targetRefs, 1] for cell in cells]

        ## 2. Compute read depth per cell
        # Cells will be distributed into X processes
        pool = mp.Pool(processes=processes)
        pool.starmap(self.read_depth_sc, cells)
        pool.close()
        pool.join()


    def read_depth_sc(self, cell, binSize, minMAPQ, filterDup, targetRefs, processes):
        '''
        Compute read depth genome wide a single cell

        Input:
            1. cell: cell instance
            2. binSize: segments size 
            3. minMAPQ: minimum read mapping quality
            4. filterDup: filter read duplicates (True) or not (False)
            5. targetRefs: list with target references
            6. processes: number of processes for parallelization
        '''
        print('read_depth_sc: ', cell.id, cell, binSize, minMAPQ, filterDup, targetRefs, processes)

        ## Create configuration dictionary
        confDict = {}
        confDict['binSize'] = binSize
        confDict['minMAPQ'] = minMAPQ        
        confDict['filterDup'] = filterDup
        confDict['targetRefs'] = targetRefs
        confDict['processes'] = processes

        ## Do depth calling
        caller = depth.read_depth_caller(cell.bam, cell.id, confDict)
        segments = caller.read_depth_wg()

        print('segments: ', len(segments), segments)
        
        ## Organize segments by chromosomes
        for segment in segments:
            cell.chromosomes[segment.ref].segments.append(segment)

class cell():
    '''
    Class representing a single cells
    '''
    def __init__(self, cellId, chrLen):
        '''
        Initialize class instance
        '''
        self.id = cellId

        ## Create chromosomes
        self.chromosomes = self.create_chromosomes(chrLen)

        ## 
        self.bam = None

    def chrom_lengths(self):
        '''
        Return dictionary with chromosome ids as keys and lengths as values
        '''        
        chromLen = dict([(chrom.id, chrom.len) for chrom in self.chromosomes.values()])
        return chromLen
            
    def create_chromosomes(self, chrLen):
        '''
        Create chromosome instances for the cell

        Input:
            1. chrLen: dictionary with a key per chromosome and its length as value
        
        Output:
            1. chromosomes: Dictionary containing chromosome ids as keys and the corresponding chromosome instance as valuew
        '''
        chromosomes = {}

        ## For each chromosome
        for chrId, length in chrLen.items():

            # Create chromosome
            chrom = chromosome(chrId, length)
            chromosomes[chrom.id] = chrom  

            # Add cell id to the chromosome
            chromosomes[chrom.id].cellId = self.id
        
        return chromosomes

class chromosome():
    '''
    Class representing a chromosome
    '''
    def __init__(self, chrId, length):
        '''
        Initialize class instance
        '''
        self.id = chrId
        self.beg = 0
        self.end = length
        self.len = length
        self.segments = []
        
    