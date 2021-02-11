'''
Module 'singleCell' - Contains classes for with single cell data (SVs, for now...)
'''

# External
import multiprocessing as mp
import sys
import math
import pandas as pd

# Internal 
from GAPI import depth 
from GAPI import log

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

    def add_cells(self, cells):
        '''
        Add cells to the population. 

        Input:
            1. cells: list of cell instances
            
        Output: Updated cells attribute containing the incorporated cells
        '''
        for cell in cells:
            self.cells[cell.id] = cell

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
    
    def depth_ratio(self, controlPop, binSize, minMAPQ, filterDup, targetRefs, processes, outDir):
        '''
        Compute read depth ratios between each cell and an average from a control cell population 

        Input:
            1. controlPop: control cell population for ratio computation
            2. binSize: segments size 
            3. minMAPQ: minimum read mapping quality
            4. filterDup: filter read duplicates (True) or not (False)
            5. targetRefs: list with target references
            6. processes: number of processes for parallelization
            7. outDir: output directory
        
        Output: 
            1. Update TDP_fc, WDP_fc, CDP_fc for each segment in every single cell. 
            2. Write tsv with genomic segments as rows, cells as columns and log2(depth_ratio) as values
        '''
        ## 1. Compute average read depths for the control population
        msg = '1. Compute average read depths for the control population'
        log.header(msg)
        c_segments = controlPop.average_read_depth(binSize, minMAPQ, filterDup, targetRefs, processes)

        ## 2. Compute read depth per cell in the population
        msg = '2. Compute read depth per cell in the population'
        log.header(msg)        
        self.read_depth(binSize, minMAPQ, filterDup, targetRefs, processes)

        ## 3. Compute depth ratios 
        msg = '3. Compute depth ratios'
        log.header(msg)            
        for cell in self.cells2list():
            cell.depth_ratio(c_segments)

        ## 4. Write read depth ratios into an output file
        msg = '4. Write read depth ratios into an output file'
        log.header(msg)           
        seriesList = []

        for cell in self.cells2list():
            
            ## Collect segment coordinates and total depth log2ratio
            tuples = [(segment.coordId(), segment.TDP_fc) for segment in cell.collect_segments()]
            indexes, values = zip(*tuples)

            ## Convert data into series
            series = pd.Series(values, index=indexes, name=cell.id)
            seriesList.append(series)

        ## Organize data into a dataframe
        df = pd.concat(seriesList, axis=1)

        ## Write to output file
        outFile = outDir + '/log2_depth_ratios.tsv'
        df.to_csv(outFile, sep='\t')


    def average_read_depth(self, binSize, minMAPQ, filterDup, targetRefs, processes):
        '''
        Compute the average read depth for each genomic bin taking into account all cells in the population

        Input:
            1. binSize: segments size 
            2. minMAPQ: minimum read mapping quality
            3. filterDup: filter read duplicates (True) or not (False)
            4. targetRefs: list with target references
            5. processes: number of processes for parallelization

        Output:
            1. consensusList: list of consensus genomic segments containing average depth across all cells
        '''
        ## 1. Compute read depth genome wide per each cell
        self.read_depth(binSize, minMAPQ, filterDup, targetRefs, processes)

        ## 2. Collect all segments across all cells
        segments = self.collect_segments()

        ## 3. Reorganize segments by grouping together segments spanning the same genomic bin across all cells 
        segment_groups = [list(i) for i in zip(*list(segments.values()))] 

        ## 4. Compute average read depth per genomic interval
        # NOTE: this may need to be parallelized
        consensusList = []

        for group in segment_groups:

            ## Create consensus segment
            consensus = depth.consensus_segment(group[0].ref, group[0].beg, group[0].end, group)
            
            ## Compute average read depths and std
            consensus.consensus_read_depth()

            ## Add consensus to the list
            consensusList.append(consensus)

        return consensusList

    def read_depth(self, binSize, minMAPQ, filterDup, targetRefs, processes):
        '''
        Compute read depth for genomic bins genome wide for each cell in the population

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
        cells = pool.starmap(self.read_depth_sc, cells)
        pool.close()
        pool.join()

        ## 3. Replace cells by new instances containing read depth counts
        self.add_cells(cells)


    def read_depth_sc(self, cell, binSize, minMAPQ, filterDup, targetRefs, processes):
        '''
        Compute read depth for genomic bins genome wide in a single cell

        Input:
            1. cell: cell instance
            2. binSize: segments size 
            3. minMAPQ: minimum read mapping quality
            4. filterDup: filter read duplicates (True) or not (False)
            5. targetRefs: list with target references
            6. processes: number of processes for parallelization
        '''
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
        
        ## Organize segments by chromosomes
        for segment in segments:
            cell.chromosomes[segment.ref].segments.append(segment)
        
        return cell

    def collect_segments(self):
        '''
        Collect genomic segments across all cells

        NOTE: this may need to be parallelized

        Output: 
            1. segments: dictionary containing cell ids as keys and the list of segments per each cell as value
        '''
        segments = {}

        ## For each cell
        for cell in self.cells.values():
            segments[cell.id] = cell.collect_segments()
        
        return segments

        
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

    def collect_segments(self):
        '''
        Collect genomic segments across all chromosomes
        '''
        segments = []

        for chrom in self.chromosomes.values():
            segments = segments + chrom.segments 

        return segments

    def depth_ratio(self, c_segments):
        '''
        Make read depth ratio between each genomic segment in the cell and the corresponding control segment
        
        Input:
            1. c_segment: list of control segments.
        
        Output: Update TDP_fc, WDP_fc and CDP_fc attributes for each segment in the cell

        ''' 
        ## 1. Collect segments 
        segments = self.collect_segments()

        ## 2. Compute read depth ratios for each segment
        for index, segment in enumerate(segments):

            ## Abort if segments don´t match
            if segment.coordId() != c_segments[index].coordId():
                print('[ERROR] Non-matching segments at depth_ratio: ', segment.coordId(), '!=', c_segments[index].coordId())
                sys.exit(1)

            ## Go ahead
            pseudocount = 0.0001
            segment.TDP_fc = math.log((segment.TDP + pseudocount) / (c_segments[index].TDP + pseudocount), 2)
            segment.WDP_fc = math.log((segment.WDP + pseudocount) / (c_segments[index].WDP + pseudocount), 2)
            segment.CDP_fc = math.log((segment.CDP + pseudocount) / (c_segments[index].CDP + pseudocount), 2)


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
        
    