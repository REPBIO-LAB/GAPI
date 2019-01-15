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
        1. FASTA: FASTA object containing consensus sequence
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
    command = 'racon -u ' + FASTQ2_file + ' ' + PAF + ' ' + FASTQ1_file + ' > ' + POLISHED
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'RACON'
        msg = 'racon failed' 
        log.step(step, msg)

    ## 5. Read polished sequence 
    FASTA = formats.FASTA()
    FASTA.read(POLISHED)

    ## 6. Use uncorrected read selected as template if no polished sequence was generated
    if not FASTA.fastaDict:
        step = 'RACON'
        msg = 'No polished sequence was generated. Use uncorrected template read as consensus' 
        log.step(step, msg)

        seqId = list(FASTQ1.fastqDict)[0]
        seq = FASTQ1.fastqDict[seqId].seq
        FASTA.fastaDict[seqId] = seq

    ## 7. Do cleanup 
    #unix.rm(FASTQ1_file)
    #unix.rm(FASTQ2_file)
    #unix.rm(PAF)
    #unix.rm(POLISHED)

    return FASTA
    