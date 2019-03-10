'''
Module 'formats' - Contains classes for dealing with file formats such as fasta, bed, vcf, etc...
'''

## DEPENDENCIES ##
# External
import itertools

# Internal
import log
import gRanges


## FUNCTIONS ##
def chromLengths_FASTA(index):
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
        self.lines = []

    def read(self, filePath):
        '''
        BED file reader. Read and store data line objects into a list:
        '''
        bedFile = open(filePath)

        # For line in the file
        for line in bedFile:
            
            # Skip comments and blank lines
            if line.startswith('#') or not line:
                continue

            fields = line.split() 
            line = BED_line(fields)
            self.lines.append(line)

class BED_line():
    '''
    BED line class 
    '''
    def __init__(self, fields):
        '''
        Initialize bed line
        '''
        self.ref = str(fields[0])
        self.beg = int(fields[1])
        self.end = int(fields[2])
        

class PAF():
    '''
    Class for dealing with files in PAF format. 
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.lines = []

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
            self.lines.append(line)

    def sortByLen(self):
        '''
        Sort alignments by query alignment length
        '''
        sortedAlignments = sorted(self.lines, key=lambda x: x.alignmentLen(), reverse=True)
        return sortedAlignments

    def chain(self):
        '''
        Chain PAF alignments based on alignment complementariety

        Output:
            1. chain
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
                maxDist = 50
                maxPercOverlap = 20
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
    def __init__(self, fields):
        '''
        Initialize paf line
        '''
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

