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

## FUNCTIONS ##

def racon(FASTQ_all, outDir):
    '''
    Build a consensus sequence from a set of long-reads (Nanopore or Pacbio). 

    The algorithm selects an arbitrary long-read sequence as template and uses the remaining sequences to correct the template
    with racon. 
    
    If long-read correction fails return the selected uncorrected long-read sequence as consensus. This should be infrequent if 
    enough number of reads provided

    Input:
        1. FASTQ_all: FASTQ object containing all the reads to be used to build a consensus 
        2. outDir: Output directory

    Output:
        1. FASTA: FASTA object containing consensus sequence or None if no consensus sequence was generated 
    '''

    ## 0. Create logs directory:
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Divide the FASTQ in two: 
    # 1) FASTQ1 containing template read to be polished
    # 2) FASTQ2 containing the remaining reads to polish the template 
    FASTQ1 = formats.FASTQ()
    FASTQ2 = formats.FASTQ()

    # Iterate over FASTQ entries
    for FASTQ_entry in FASTQ_all.fastqDict.values():

        # a) Select read to be polished
        if not FASTQ1.fastqDict:
            FASTQ1.add(FASTQ_entry)

        # b) Collect remaining reads to polish
        else:
            FASTQ2.add(FASTQ_entry)

    ## 2. Write FASTQ files 
    FASTQ1_file = outDir + '/template.fastq'
    FASTQ1.write(FASTQ1_file)

    FASTQ2_file = outDir + '/reads2polish.fastq'
    FASTQ2.write(FASTQ2_file)

    ## 3. Align reads against the read to be corrected
    PAF = outDir + '/alignments.paf'
    err = open(logDir + '/minimap2.err', 'w') 
    command = 'minimap2 -x ava-ont ' + FASTQ1_file + ' ' + FASTQ2_file + ' > ' + PAF
    status = subprocess.call(command, stderr=err, shell=True)

    ## 4. Read polishing with racon
    POLISHED = outDir + '/polished.fasta'
    err = open(logDir + '/racon.err', 'w') 
    command = 'racon ' + FASTQ2_file + ' ' + PAF + ' ' + FASTQ1_file + ' > ' + POLISHED
    status = subprocess.call(command, stderr=err, shell=True)

    ## 5. Read polished sequence 
    FASTA = formats.FASTA()
    FASTA.read(POLISHED)

    ## 6. Set FASTA as None if no consensus sequence was generated
    if not FASTA.fastaDict:
        FASTA = None

    ## 7. Do cleanup 
    #unix.rm(FASTQ1_file)
    #unix.rm(FASTQ2_file)
    #unix.rm(PAF)
    #unix.rm(POLISHED)

    return FASTA
    