'''
Module 'singleCell' - Contains classes for with single cell data (SVs, for now...)
'''

# External

# Internal 


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
        
    