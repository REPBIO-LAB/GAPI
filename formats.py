'''
Module 'formats' - Contains classes for dealing with file formats such as fasta, bed, vcf, etc...
'''

## DEPENDENCIES ##
import log


## CLASSES ##
class fasta():
    '''
    Class for dealing with files in FASTA format. Read and writes fastas
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.fastaDict = {}

    def read(self, fastaFile):
        '''
        FASTA file reader. Read fasta file and store data into a dictionary with the following format:
            - Key. Sequence identifier
            - Value. Sequence 
        '''
        fastaDict = {}

        log.subHeader('FASTA reader')

        fh = open(fastaFile)
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == '>'))
        for header in faiter:
            # drop the >
            header = header.next()[1:].strip()

            # drop the info
            header = header.split(' ')[0]

            log.info('Reading ' + header + '...')

            # join all sequence lines to one.
            seq = ''.join(s.strip() for s in next(faiter))
            fastaDict[header] = seq

        self.fastaDict = fastaDict

    def write(self, outFilePath):
        '''
        FASTA file writer. Write data stored in the dictionary into a FASTA file
        '''
        outFile = open(outFilePath, 'w')

        for header, seq in self.fastaDict.items():
            header = '>' + header

            outFile.write("%s\n" % header)
            outFile.write("%s\n" % seq)

        # Close output fasta file
        outFile.close()

