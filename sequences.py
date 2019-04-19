'''
Module 'sequences' - Contains functions for the manipulation and extracting information from nucleotidic sequences
'''

## DEPENDENCIES ##
# External
import os
import subprocess

# Internal
import formats


## FUNCTIONS ##
def rev_complement(seq):
    '''
    Make the reverse complementary of a dna sequence

    Input:
        1. seq: DNA sequence
    Output:
        1. revComplementSeq: Reverse complementary of input DNA sequence
    '''
    baseComplementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    seq = seq.upper()
    revSeq = seq[::-1] # Make reverse sequence
    letters = list(revSeq)
    letters = [baseComplementDict[base] for base in letters]
    revComplementSeq = ''.join(letters) # Make complement of reverse sequence

    return revComplementSeq

def baseComposition(seq):
    '''
    Compute mononucleotide frequencies for an input sequence

    Input:
        1. seq: DNA sequence
    Output:
        1. baseCounts: Dictionary containing counts per nucleotidic base for the input sequence 
    '''
    baseCounts = {}
    baseCounts['A'] = seq.count('A') 
    baseCounts['G'] = seq.count('G') 
    baseCounts['C'] = seq.count('C') 
    baseCounts['T'] = seq.count('T') 
    baseCounts['N'] = seq.count('N')
    baseCounts['total'] = baseCounts['A'] + baseCounts['G'] + baseCounts['C'] + baseCounts['T'] + baseCounts['N']

    return baseCounts


def find_monomers(seq, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity):
    '''
    Find nucleotide monomers in an input sequence

    Input:
        1. seq: DNA sequence
        2. targetMonomer: Type of monomer to search for
        3. windowSize: Window size used to scan the sequence
        4. maxWindowDist: Maximum number of windows without the monomer between two monomeric windows to group them in the same cluster
        5. minMonomerSize: Minimum size for calling a monomer
        6. minPurity: minimum % of bases corresponding to the target monomer in a window or monomer cluster

    Output:
        1. monomers: List containing the monomer instances identified 
    '''
    seqLen = len(seq)

    windowDist = 0
    monomers = []
    noMonomeric = []
    indexes = []

    ## 1. Split seq into consecutive slices 
    # ---ws--->---ws--->---ws--->---ws--->
    slices = [(pos, pos + windowSize, seq[pos:pos + windowSize]) for pos in range(0, seqLen, windowSize)]

    ## 2. Search for candidate monomers by iterating over the slices
    for index, value in enumerate(slices):
        
        indexes.append(index)

        # Collect slide info
        beg, end, sliceSeq = value
        sliceLen = len(sliceSeq)

        # Compute % of target monomer in the slice 
        baseCounts = baseComposition(sliceSeq)
        percMonomer = float(baseCounts[targetMonomer])/sliceLen*100
        
        ## A) Slice correspond to a monomer (---ws--- == AAA.../TTT.../...) if:
        # 1. sliceLen >= windowSize AND 
        # 2. purity >= minPurity
        if (sliceLen >= windowSize) and (percMonomer >= minPurity):

            ## a) CREATE new monomer object if:
            #  a.1) No monomer already identified OR
            #  a.2) Distance from new monomer to the previously reported monomer > maxWindowDist 
            if (not monomers) or (windowDist > maxWindowDist):

                ## Create monomer
                monomerObj = monomer(beg, sliceSeq, index)
                monomers.append(monomerObj)

            ## b) EXTENSION if monomer within maxWindowDist to previous identified monomer 
            else:
                previousMonomer = monomers[-1]
                
                ## b.1 Extend previously identified monomer with all the not monomeric slices till the new monomer 
                for noMonomericSlice in noMonomeric:
                    noMonomericSliceSeq, noMonomericIndex = noMonomericSlice[2:]
                    previousMonomer.add(noMonomericSliceSeq, noMonomericIndex)

                ## b.2 Extend previously identified monomer with the new monomer         
                previousMonomer.add(sliceSeq, index)

            ## Initiate distance counter and no monomeric list
            windowDist = 0
            noMonomeric = []

        ## B) Slide do not correspond to a monomer
        else:
            windowDist += 1
            noMonomeric.append([beg, end, sliceSeq, index])
    
    ## 3. Refine monomers by attempting to extend both ends 
    # [ACGAA] AAAAAAAAAAAAAAAAAAAA [AAAGG]
    # Look to the slices flanking the monomer and extend till no A match is found
    # In the example the monomer would be extended as follows:
    # [ACG|AA]AAAAAAAAAAAAAAAAAAAA[AAA|GG]
    # '|' indicates the end of the extension 
    
    # For each monomer
    for monomerObj in monomers:
    
        ## 3.1 Left extension
        # [ACG|AA]AAAAAAAAAAAAAAAAAAAA
        # '|' indicates the end of the extension 
        index = monomerObj.indexes[0] - 1 

        # flanking slice available
        if index in indexes:

            # Pick sequence
            sliceSeq = slices[index][2]

            # Parse sequence backward one nucleotide at a time
            # [ACG|AA] *Monomer* 
            # <-<-<-<-
            for nucleotide in reversed(sliceSeq):

                # a) Extend monomer if nucleotide match
                if nucleotide == targetMonomer:
                    
                    monomerObj.seq = nucleotide + monomerObj.seq
                    monomerObj.beg -= 1

                    # Add slice index to the list
                    if index not in monomerObj.indexes:
                        monomerObj.indexes.insert(0, index)
 
                
                # b) Stop extension once no nucleotide match is found
                else:
                    break

        ## 3.2 Right extension
        # AAAAAAAAAAAAAAAAAAAA[AAA|GG]
        # '|' indicates the end of the extension 
        index = monomerObj.indexes[-1] + 1 
        
        # flanking slice available
        if index in indexes:

            # Pick sequence
            sliceSeq = slices[index][2]

            # Parse sequence forward one nucleotide at a time
            # *Monomer* [AAA|GG]
            #           ->->->->
            for nucleotide in sliceSeq:

                # a) Extend monomer if nucleotide match
                if nucleotide == targetMonomer:
                    
                    monomerObj.seq = monomerObj.seq + nucleotide 
                    monomerObj.end += 1

                    # Add slice index to the list
                    if index not in monomerObj.indexes:
                        monomerObj.indexes.append(index) 

                # b) Stop extension once no nucleotide match is found
                else:
                    break            
            
    ## 4. Filter candidate monomers 
    filteredMonomers = []

    # For each monomer
    for monomerObj in monomers:
        monomerSize = len(monomerObj.seq)
        baseCounts = baseComposition(monomerObj.seq)
        percMonomer = float(baseCounts[targetMonomer])/monomerSize*100

        ## Select those monomers with:
        # 1. size >= minMonomerSize AND 
        # 2. purity >= minPurity
        if (monomerSize >= minMonomerSize) and (percMonomer >= minPurity):
            filteredMonomers.append(monomerObj)
        
    return filteredMonomers


def aligmentMaxNbMatches(FASTA_file, db, PAF_file, outDir):
    '''
    '''

    # DESILENCIAAAAAR!!!
    # TODO: append en el error!
    err = open(outDir + '/identifyMate.err', 'w')
    command = 'minimap2 ' + db + ' ' + FASTA_file + ' > ' + PAF_file
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'IDENTIFY MATE SEQ'
        msg = 'Identify mate sequence failed' 
        log.step(step, msg)
        


    # If PAF file is not empty
    if not os.stat(PAF_file).st_size == 0:
        PAFObj = formats.PAF()
        PAFObj.read(PAF_file)

        # Pick the identity of the aligment with highest number of matches
        aligmentMaxNbMatches = PAFObj.sortNbMatches()[0]

    else:
        aligmentMaxNbMatches = None

    return aligmentMaxNbMatches


## CLASSES ##
class monomer():
    '''
    '''

    def __init__(self, beg, seq, index):
        '''
        Initialize empty class instance
        '''
        self.beg = beg
        self.end = self.beg + len(seq)
        self.seq = seq
        self.indexes = [index]

    def add(self, seq, index):
        '''
        Extend monomer by adding a new piece of sequence at its end

        Input:
            1. seq: piece of sequence to add at monomer end
            2. index: slice index
        '''
        self.end = self.end + len(seq)
        self.seq = self.seq + seq
        self.indexes.append(index)
    
    def length(self):
        '''
        Return monomer length
        '''
        return len(self.seq)

## [SR CHANGE]
def getConsensusSeq(FASTA_file, outDir):

    ### 2. Make multiple sequence alignment
    msfPath = FASTA_file.replace("fa", "msf")
    command = 'muscle -in ' + FASTA_file + ' -out ' + msfPath + ' -msf' 
    status = subprocess.call(command, shell=True)

    ### 3. Generate consensus sequence (cons tool from EMBOSS packagge)
    consensusPath = FASTA_file.replace("_supportingReads", "_consensus")

    command = 'cons -sequence ' + msfPath + ' -outseq ' + consensusPath + ' -identity 0 -plurality 0'
    status = subprocess.call(command, shell=True)

    ### Read consensus sequence 
    consensusFastaObj = formats.FASTA()
    consensusFastaObj.read(consensusPath)
    consensusSeq = consensusFastaObj.seqDict["EMBOSS_001"].upper()

    # TODO
    ### Do cleanup
    #command = 'rm ' + fastaPath + ' ' + msfPath + ' ' + consensusPath             
    #os.system(command) # returns the exit status

    ## Replace '-' by 'N' for ambiguous bases:
    consensusSeq = consensusSeq.replace('-', 'N')

    ## Convert consensus sequence into upper case:
    consensusSeq = consensusSeq.upper()

    return consensusPath, consensusSeq

## [SR CHANGE]
def getPAFAlign(FASTA_file, indexDb, outDir):
    # Alineo el fasta consenso
    # TODO ponerlo bien
    PAF_file = FASTA_file.replace(".fa", "_alignments.paf")

    #err = open(logDir + '/align.err', 'w') 
    command = 'minimap2 ' + indexDb + ' ' + FASTA_file + ' > ' + PAF_file
    status = subprocess.call(command, shell=True)

    if status != 0:
        step = 'ALIGN-INSERT'
        msg = 'Insert alignment failed' 
        log.step(step, msg)
    
    return PAF_file