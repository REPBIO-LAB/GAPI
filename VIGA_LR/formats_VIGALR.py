'''
Module 'formats_VIGALR' - Contains classes for dealing with file formats such as fasta, bed, vcf, etc...
'''

## DEPENDENCIES ##
# External
import itertools
import sys
import time
import re

# Internal
import log
import gRanges
import structures
from VIGA_LR import call_VIGALR

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

    def write(self, filePath, mode = 'write', safetyLock = False):
        '''
        FASTA file writer. Write data stored in the dictionary into a FASTA file
        Mode: write -> write new file. append -> append to existing file or create if tit doesnt exist.
        '''
        openMode = 'a' if mode == 'append' else 'w'
        
        if safetyLock:
            call_VIGALR.lock.acquire()

        fastaFile = open(filePath, openMode)

        for header, seq in self.seqDict.items():
            header = '>' + header

            fastaFile.write("%s\n" % header)
            fastaFile.write("%s\n" % seq)

        # Close output fasta file
        fastaFile.close()
        
        if safetyLock:
            call_VIGALR.lock.release()

    def retrieve_seqs(self, targetNames):
        '''
        Retrieve set of sequences from fasta file

        Input:
            1. targetNames: list of read ids to be retrieved
        
        Output:
            2. outDict: dictionary containing sequences
        '''
        outDict = {readName: self.seqDict[readName] for readName in targetNames if readName in self.seqDict}
        
        if safetyLock:
            callers.lock.release()

    def retrieve_seqs(self, targetNames):
        '''
        Retrieve set of sequences from fasta file

        Input:
            1. targetNames: list of read ids to be retrieved
        
        Output:
            2. outDict: dictionary containing sequences
        '''
        outDict = {readName: self.seqDict[readName] for readName in targetNames if readName in self.seqDict}
        
        return outDict