'''
Module 'formats' - Contains classes for dealing with file formats such as fasta, bed, vcf, etc...
'''

## DEPENDENCIES ##
# External
import itertools
import sys
import formats

# Internal
import log
import gRanges
import structures

## FUNCTIONS ##
def chrom_lengths_index(index):
    '''
    Read FASTA index and build a dictionary containing chromosome lengths

    Input:
        1. index: FASTA file index generated with samtools faidx
    Output:
        1. chromLen: Dictionary containing the length for each chromosome 
    '''

    chromLengths = {}

    with open(index) as indexFile:
        for line in indexFile:
            line = line.rstrip().split("\t")
            ref = line[0]
            length = line[1]
            chromLengths[ref] = int(length)

    return chromLengths

def merge_FASTA(FASTA_list):
    '''
    Merge a list of FASTA objects into a single one

    Input:
        1. FASTA_list: List of FASTA objects to be merged
    Output:
        1. FASTA_merged: FASTA object resulting from merging
    '''
    ## Initialize output FASTA
    FASTA_merged = FASTA()

    ## For each FASTA object
    for FASTA_obj in FASTA_list:

        # For each sequence in the FASTA
        for seqId, seq in FASTA_obj.seqDict.items():
                    
            # Add sequence to the merged FASTA
            FASTA_merged.seqDict[seqId] = seq 

    return FASTA_merged

def bed2binDb(bedPath, refLengths, threads):
    '''
    Organize features in a bed file into a whole genome bin database. 
    
    Bed file must contain at least 4 fields: ref, beg, end and name (feature name). Extra fields
    will not be considered

    Input:
        1. bedPath: path to bed file
        2. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        3. threads: number of threads used to parallelize the bin database creation

    Output:
        1. wgBinDb: dictionary containing references as keys and the corresponding 'bin_database' as value
    '''
    ## Read bed
    bed = formats.BED()
    targetRefs = list(refLengths.keys())
    bed.read(bedPath, 'nestedDict', targetRefs)

    ## Create bin database
    wgBinDb = structures.create_bin_database(refLengths, bed.lines, threads)

    return wgBinDb

## CLASSES ##
class FASTA():
    '''
    Class for dealing with files in FASTA format
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.seqDict = {}

    def read(self, filePath):
        '''
        FASTA file reader. Read and add data to a dictionary with the following format:
            - Key. Sequence identifier
            - Value. Sequence 
        '''
        fastaFile = open(filePath)

        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in itertools.groupby(fastaFile, lambda line: line[0] == '>'))

        for header in faiter:
            # drop the >
            header = header.__next__()[1:].strip()
            
            # drop the info
            header = header.split(' ')[0]

            # join all sequence lines to one.
            seq = ''.join(s.strip() for s in next(faiter))
            self.seqDict[header] = seq

    def write(self, filePath):
        '''
        FASTA file writer. Write data stored in the dictionary into a FASTA file
        '''
        fastaFile = open(filePath, 'w')

        for header, seq in self.seqDict.items():
            header = '>' + header

            fastaFile.write("%s\n" % header)
            fastaFile.write("%s\n" % seq)

        # Close output fasta file
        fastaFile.close()
        

class FASTQ():
    '''
    Class for dealing with files in FASTQ format
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.seqDict = {}

    def write(self, filePath):
        '''
        FASTQ file writer. Write data stored in the dictionary into a FASTA file
        '''
        # Open fastq file
        fastqFile = open(filePath, 'w')

        # Write FASTQ lines
        for FASTQ_entry in self.seqDict.values(): 

            seqId = '@' + FASTQ_entry.seqId
            description = '+' + FASTQ_entry.description

            fastqFile.write("%s\n" % seqId)
            fastqFile.write("%s\n" % FASTQ_entry.seq)
            fastqFile.write("%s\n" % description)
            fastqFile.write("%s\n" % FASTQ_entry.qual)

        # Close output fasta file
        fastqFile.close()

    def add(self, fastqLine):
        '''
        Add sequence to the fastq object
        '''
        self.seqDict[fastqLine.seqId] = fastqLine
        

class FASTQ_entry():
    '''
    FASTQ entry class 
    '''

    def __init__(self, seqId, seq, description, qual):
        '''
        Initialize fastq line
        '''
        self.seqId = seqId
        self.seq = seq
        self.description = description
        self.qual = qual

class BED():
    '''
    Class for dealing with files in BED format. 
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.lines = None
        self.structure =  None

    def read(self, filePath, structure, targetRefs):
        '''
        BED file reader. Read and store bed lines into a data structure

        Input:
            1) filePath: path to bed file
            2) structure: data structure where BED lines will be stored. 3 Possible structures:
                - 'List': lines saved in a list
                - 'Dict': lines saved in a dictionary where each key will correspond to a reference and the corresponding value will be the list of lines in that reference
                - 'nestedDict': lines saved in a nested dictionary where first level keys correspond to the references and second level keys to bed entry names

            3) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded
        Initialize lines attribute as output
        '''
        self.structure = structure

        # a) Organize BED entries into a list
        if (self.structure == 'List'):
            self.lines = self.organize_list(filePath, targetRefs)

        # b) Organize BED entries into a dict
        elif (self.structure == 'Dict'):
            self.lines = self.organize_dict(filePath, targetRefs)

        # c) Organize BED entries into a nested dict
        elif (self.structure == 'nestedDict'):
            self.lines = self.organize_nestedDict(filePath, targetRefs)

        # d) Unkown data type structure provided
        else:
            print('[ERROR] Bed file reader. Unknown structure provided: ', structure)
            sys.exit(1)

    def write(self, outPath):
        '''
        BED file writer. Write bed entries into bed file

        Input:
            1) outPath: path to output bed file

        NOTE: Update to include optional fields
        '''
        ## Collect all the entries into a list
        # a) Entries organized into a list
        if (self.structure == 'List'):
            outEntries = self.lines

        # b) Entries organized into a dict (TO TEST LATER)
        elif (self.structure == 'Dict'):
            outEntries = [line for line in self.lines.values()]

        # c) Entries organized into a nested dict (TO IMPLEMENT LATER)
        #elif (self.structure == 'nestedDict'):

        ## Write entries into output bed file
        with open(outPath, 'w') as outFile:
            for entry in outEntries:

                fields = [entry.ref, str(entry.beg), str(entry.end)]

                # Add name if available
                if hasattr(entry, 'name'):
                    fields.append(entry.name)

                # Create row
                row = "\t".join(fields)
                outFile.write(row + '\n')

    def write_annovar(self, outPath):
        '''
        BED file writer. Write bed entries into a file format that can be used as input for annovar

        Input:
            1) outPath: path to output bed file
        '''
        ## Collect all the entries into a list
        # a) Entries organized into a list
        if (self.structure == 'List'):
            outEntries = self.lines

        # b) Entries organized into a dict (TO TEST LATER)
        elif (self.structure == 'Dict'):
            outEntries = [line for line in self.lines.values()]

        # c) Entries organized into a nested dict (TO IMPLEMENT LATER)
        #elif (self.structure == 'nestedDict'):

        ## Write entries into output bed file
        with open(outPath, 'w') as outFile:
            for entry in outEntries:

                fields = [entry.ref, str(entry.beg), str(entry.end), '0', '0']

                # Add name if available
                if 'name' in entry.optional:
                    fields.append('comments: ' + entry.optional['name'])

                # Create row
                row = "\t".join(fields)
                outFile.write(row + '\n')

    def organize_list(self, filePath, targetRefs):
        '''
        Organize bed file lines in a list

        Input:
            1) filePath: path to bed file
            2) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded
        
        Output:
            1) lines: list containing bed entries
        '''
        bedFile = open(filePath)
        lines = []
        header = []
        
        # For line in the file
        for line in bedFile:
            
            # Skip blank lines
            if not line:
                continue
            
            # Split data line into fields
            fields = line.split()                 

            # A) Header
            if line.startswith('#'):
                header = fields

            # B) Data line
            else:
                line = BED_line(fields, header)

                if (targetRefs is None) or (line.ref in targetRefs):
                    lines.append(line)
        
        return lines

    def organize_dict(self, filePath, targetRefs):
        '''
        Organize bed file lines in a dictionary where each key will correspond to a reference and the corresponding value will be the list of lines in that reference
        
        Input:
            1) filePath: path to bed file
            2) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded

        Output:
            1) lines: dictionary containing bed entries
        '''

        bedFile = open(filePath)
        lines = {}
        header = []

        # For line in the file
        for line in bedFile:
            
            # Skip blank lines
            if not line:
                continue
            
            # Split data line into fields
            fields = line.split()       

            # A) Header
            if line.startswith('#'):
                header = fields

            # B) Data line
            else:
                line = BED_line(fields, header)

                if (targetRefs is None) or (line.ref in targetRefs):

                    # a) Initialize reference and add line
                    if line.ref not in lines:
                        lines[line.ref] = [line]

                    # b) Add to preexisting reference
                    else:
                        lines[line.ref].append(line)

        return lines

    def organize_nestedDict(self, filePath, targetRefs):
        '''
        Organize bed file lines in a nested dictionary where first level keys correspond to the references and second level keys to bed entry names  

        Input:
            1) filePath: path to bed file
            2) targetRefs: list containing target references. Entries with refs not included in the list will not be loaded. 
                            If 'None' all the entries will be loaded

        Output:
            1) lines: nested dictionary containing bed entries

        NOTE: I need to update the code to take into account the bed header and optional fields
        '''

        bedFile = open(filePath)
        lines = {}
        header = []

        # For line in the file
        for line in bedFile:
            
            # Skip blank lines
            if not line:
                continue
            
            # Split data line into fields
            fields = line.split()       

            # A) Header
            if line.startswith('#'):
                header = fields

            # B) Data line
            else:
                line = BED_line(fields, header)
                
                if (targetRefs is None) or (line.ref in targetRefs):

                    # A) Initialize reference and add line
                    if line.ref not in lines:
                        lines[line.ref] = {}
                        lines[line.ref][line.optional['name']] = [line]

                    # B) Add to preexisting reference
                    else:

                        # a) Initialize entry name
                        if line.optional['name'] not in lines[line.ref]:
                            lines[line.ref][line.optional['name']] = [line]

                        # b) Add to preexisting entry name
                        else:
                            lines[line.ref][line.optional['name']].append(line)

        return lines

class BED_line():
    '''
    BED line class 
    '''
    number = 0 # Number of instances

    def __init__(self, fields, header):
        '''
        Initialize bed line

        Input:
            1. fields: list containing a bed feature (== data line)
            2. header: list containing bed header (required for parsing optional fields)
        '''
        BED_line.number += 1 # Update instances counter
        self.id = 'BED_line_' + str(BED_line.number)

        ## Mandatory fields
        self.ref = str(fields[0])
        self.beg = int(fields[1])
        self.end = int(fields[2])
        self.clusterId = None

        ## Optional fields dictionary (Optional fields only considered if header provided)
        self.optional = {}

        for i in range(3, len(header), 1):
            self.optional[header[i]] = fields[i]
        

class PAF():
    '''
    Class for dealing with files in PAF format. 
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.alignments = []

    def read(self, filePath):
        '''
        PAF file reader. Read and store data line objects into a list:
        '''
        pafFile = open(filePath)

        # For line in the file
        for line in pafFile:
            
            # Skip comments and blank lines
            if line.startswith('#') or not line:
                continue

            fields = line.split() 
            line = PAF_line(fields)
            self.alignments.append(line)

    def sortByLen(self):
        '''
        Sort alignments by query alignment length in descending order
        '''
        sortedAlignments = sorted(self.alignments, key=lambda x: x.alignmentLen(), reverse=True)
        return sortedAlignments

    ## [SR CHANGES]
    def sortNbMatches(self):
        '''
        Sort alignments by query alignment length
        '''
        sortedAlignments = sorted(self.alignments, key=lambda x: x.nbMatches, reverse=True)
        return sortedAlignments

    def chain(self, maxDist, maxPercOverlap):
        '''
        Chain PAF alignments based on alignment complementariety

        Input:
            1. maxDist: maximum distance between both ranges 
            2. maxPercOverlap: maximum percentage of overlap between ranges

        Output:
            1. chain: PAF_chain object instance
        '''
        ## 1. Sort alignments by decreasing query alignment length 
        sortedAlignments = self.sortByLen()

        ## 2. Pick longest alignment and initiate chain
        longest = sortedAlignments[0]
        chain = PAF_chain([longest])
        
        # remove alignment 
        del sortedAlignments[0]

        roundCounter = 1

        ## 3. Attemp to extend the chain with complementary alignments
        while True:

            # START ALIGNMENT CHAIN EXTENSION ROUND
            # Initialize boolean as not complementary alignment found
            complBool = False

            # Go through all the available alignments 
            for index, alignment in enumerate(sortedAlignments):

                ## Assess if alignment complementary to the chain
                chainBeg, chainEnd = chain.interval()
                complBool, orientation = gRanges.complementary(chainBeg, chainEnd, alignment.qBeg, alignment.qEnd, maxDist, maxPercOverlap)

                ## Complementary alignment found 
                if complBool:

                    ## Add alignment to the chain 
                    # a) Add to the chain begin
                    if orientation == "LEFT":
                        chain.alignments.insert(0,alignment)

                    # b) Add to the chain end
                    else:
                        chain.alignments.append(alignment)

                    ## Remove from list
                    del sortedAlignments[index]

                    ## Stop once complementary found
                    break

            roundCounter += 1

            # STOP CHAIN EXTENSION IF: 
            # a) No complementary alignment found in the last round OR
            # b) Mo alignments left
            if complBool == False or not sortedAlignments:
                break

        return chain

class PAF_line():
    '''
    PAF line class 
    '''
    number = 0 # Number of instances

    def __init__(self, fields):
        '''
        Initialize paf line
        '''
        PAF_line.number += 1 # Update instances counter
        self.id = 'PAF_line_' + str(PAF_line.number)
        self.qName = str(fields[0])
        self.qLen = int(fields[1])
        self.qBeg = int(fields[2])
        self.qEnd = int(fields[3])
        self.strand = str(fields[4])
        self.tName = str(fields[5])
        self.tLen = int(fields[6])
        self.tBeg = int(fields[7])
        self.tEnd = int(fields[8])
        self.nbMatches = int(fields[9])
        self.blockLen = int(fields[10])
        self.MAPQ = int(fields[11])   

    def alignmentLen(self):
        '''
        Compute the query alignment length
        '''

        return self.qEnd - self.qBeg

class PAF_chain():
    '''    
    Chain of complementary PAF alignments  
    '''

    def __init__(self, alignments):
        '''
        Initialize chain instance. 
        
        Input:
            1. alignments. List of PAF_line instances
        '''
        self.alignments = alignments

    def interval(self):
        '''
        Return query interval covered by the chain
        '''
        firstAlignment = self.alignments[0]
        lastAlignment = self.alignments[-1]
        
        return firstAlignment.qBeg, lastAlignment.qEnd
    
    def perc_query_covered(self):
        '''
        Compute the percentage of the query sequence covered by the chain of alignments
        '''
        # a) No alignments available
        if len(self.alignments) == 0:
            percCovered = 0

        # b) Alignments available
        else:

            ## Compute the number of bases covered
            beg, end = self.interval()
            alignmentLen = end - beg

            ## Compute the percentage of bases covered
            percCovered = float(alignmentLen)/self.alignments[0].qLen*100

        return percCovered

