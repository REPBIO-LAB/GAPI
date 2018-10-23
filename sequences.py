'''
Module 'sequences' - Contains classes and functions for the manipulation of sequences
'''

## DEPENDENCIES ##
import pysam


## FUNCTIONS ##
def rev_complement(seq):
    '''
        Make the reverse complementary of a dna sequence
        Input:
            - seq: DNA sequence
        Output:
            - revComplementSeq: Reverse complementary of input DNA sequence
    '''
    baseComplementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    seq = seq.upper()
    revSeq = seq[::-1] # Make reverse sequence
    letters = list(revSeq)
    letters = [baseComplementDict[base] for base in letters]
    revComplementSeq = ''.join(letters) # Make complement of reverse sequence

    return revComplementSeq

