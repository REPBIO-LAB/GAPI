'''
Module 'formats' - Contains classes for dealing with file formats such as fasta, bed, vcf, etc...
'''

## DEPENDENCIES ##
import log


## CLASSES ##
class fasta():
    '''
    Class for dealing with files in FASTA format
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.fastaDict = {}

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
            header = header.next()[1:].strip()
            
            # drop the info
            header = header.split(' ')[0]

            # join all sequence lines to one.
            seq = ''.join(s.strip() for s in next(faiter))
            self.fastaDict[header] = seq

    def write(self, filePath):
        '''
        FASTA file writer. Write data stored in the dictionary into a FASTA file
        '''
        fastaFile = open(filePath, 'w')

        for header, seq in self.fastaDict.items():
            header = '>' + header

            fastaFile.write("%s\n" % header)
            fastaFile.write("%s\n" % seq)

        # Close output fasta file
        fastaFile.close()


class bed():
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

            ref, beg, end = line.split()[:3]
            line = bedLine(ref, beg, end)
            self.lines.append(line)


class bedLine():
    '''
    BED line class 
    '''

    def __init__(self, ref, beg, end):
        '''
        Initialize bed line
        '''
        self.ref = str(ref)
        self.beg = int(beg)
        self.end = int(end)
