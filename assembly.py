'''
Module 'consensus' - Contains funtions to create a consensus sequence from a set of input sequences
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
import log
import formats
import unix
import alignment
import os

## FUNCTIONS ##

def racon(reads, technology, quality, outDir):
    '''
    Build a consensus sequence from a set of sequencing reads  

    The algorithm selects an arbitrary read sequence as template and uses the remaining sequences to correct the template
    with racon. 
    
    Input:
        1. reads: FASTQ (if qualities available) or FASTA (if qualities NOT available) object containing all the reads to be used to build a consensus 
        2. technology: Sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        3. quality: True (sequence qualities available) or False (not available).
        4. outDir: Output directory

    Output:
        1. FASTA: FASTA object containing consensus sequence or None if no consensus sequence was generated 
    '''

    ## 0. Create logs directory:
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Divide the FASTQ/FASTA in two: 
    # 1) template: containing template read to be polished (template arbitrarily selected)
    # 2) reads2polish: containing the remaining reads to polish the template 
    
    # a) Quality available -> Create FASTQ files
    if quality:
        template = formats.FASTQ()
        reads2polish = formats.FASTQ()

        # Iterate over FASTQ entries
        for entry in reads.seqDict.values():
            # a) Select read to be polished
            if not template.seqDict:
                template.add(entry)

            # b) Collect remaining reads to polish
            else:
                reads2polish.add(entry)

        ## Write FASTQ files 
        template_file = outDir + '/template.fastq'
        template.write(template_file)

        reads2polish_file = outDir + '/reads2polish.fastq'
        reads2polish.write(reads2polish_file)

    # b) Quality not available > Create FASTA files
    else:    
        template = formats.FASTA()
        reads2polish = formats.FASTA()

        ## Iterate over FASTA entries
        for readId, seq in reads.seqDict.items():

            # a) Select read to be polished
            if not template.seqDict:
                template.seqDict[readId] = seq

            # b) Collect remaining reads to polish
            else:
                reads2polish.seqDict[readId] = seq

        ## Write FASTA files 
        template_file = outDir + '/template.fasta'
        template.write(template_file)

        reads2polish_file = outDir + '/reads2polish.fasta'
        reads2polish.write(reads2polish_file)

    ## 3. Align reads against the template 
    ## Set preset according to the technology
    preset = alignment.minimap2_presets(technology)

    ## Do alignment 
    PAF = outDir + '/alignments.paf'
    err = open(logDir + '/minimap2.err', 'w') 
    command = 'minimap2 -x ' + preset + ' ' + template_file + ' ' + reads2polish_file + ' > ' + PAF
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'CONSENSUS'
        msg = 'Alignment of sequences against template failed' 
        log.step(step, msg)

    ## 4. Template polishing with racon
    POLISHED = outDir + '/polished.fasta'
    err = open(logDir + '/racon.err', 'w') 
    command = 'racon ' + reads2polish_file + ' ' + PAF + ' ' + template_file + ' > ' + POLISHED
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'CONSENSUS'
        msg = 'Template polishing failed' 
        log.step(step, msg)

    ## 5. Read polished sequence 
    FASTA = formats.FASTA()
    FASTA.read(POLISHED)

    ## 6. Set FASTA as None if no consensus sequence was generated
    if not FASTA.seqDict:
        FASTA = None

    return FASTA

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

    if not os.stat(consensusPath).st_size == 0:
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
    else:
        consensusSeq = None

    return consensusPath, consensusSeq