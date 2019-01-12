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
    ...

    Input:
        1. FASTQ_all: ....
        2. outDir: Output directory
    '''

    ## 0. Create logs directory:
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Divide the FASTQ in two: 
    # 1) FASTQ1 containing one read to be polished
    # 2) FASTQ2 containing the remaining reads to polish the read in FASTQ1 
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
    FASTQ1_file = outDir + '/targetRead.fastq'
    FASTQ1.write(FASTQ1_file)

    FASTQ2_file = outDir + '/reads2polish.fastq'
    FASTQ2.write(FASTQ2_file)

    ## 3. Align reads against read to be corrected
    PAF = outDir + '/alignments.paf'
    err = open(logDir + '/minimap2.err', 'w') 
    command = 'minimap2 -x ava-ont ' + FASTQ1_file + ' ' + FASTQ2_file + ' > ' + PAF
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'MINIMAP2'
        msg = 'minimap2 alignment failed' 
        log.step(step, msg)

    ## 4. Read polishing with racon
    POLISHED = outDir + '/polished.fasta'
    err = open(logDir + '/racon.err', 'w') 
    command = 'racon ' + FASTQ2_file + ' ' + PAF + ' ' + FASTQ1_file + ' > ' + POLISHED
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'RACON'
        msg = 'racon failed' 
        log.step(step, msg)

    ## 5. Read polished sequence and return
    FASTA = formats.FASTA()
    FASTA.read(POLISHED)

    return FASTA
    