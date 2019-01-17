'''
Module 'alignment' - Contains funtions align sequences
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
import log
import formats
import unix
import bamtools 

## FUNCTIONS ##

def targeted_alignment_minimap2(FASTA, targetInterval, reference, outDir):
    '''
    Align a set of sequences into a reference genome target region. 
    
    Useful for doing local realignment of reads around SV breakpoints. Much faster than whole genome realignment

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. targetInterval: Reference genome interval where sequences will be aligned. The interval must be provided as chr:beg-end.
        3. reference: Path to the reference genome in fasta format. An index of the reference generated with samtools faidx must be located in the same directory
        4. outDir: Output directory
        
    Output:
        1. BAM_sorted: Path to sorted BAM file containing input sequences alignments 
    '''
    ## 0. Create logs directory
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Extract the reference target region prior alignment 
    target = outDir + '/target.fasta'
    err = open(logDir + '/target.err', 'w') 
    command = 'samtools faidx ' + reference + ' ' + targetInterval + ' > ' + target
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'TARGET'
        msg = 'Extraction of reference target region failed' 
        log.step(step, msg)

    ## 2. Align the sequences into the target region 
    SAM = outDir + '/alignments.sam'
    err = open(logDir + '/align.err', 'w') 
    command = 'minimap2 -ax map-ont ' + target + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)

    ## 3. Convert SAM to sorted BAM
    BAM = bamtools.SAM2BAM(SAM, outDir)

    return BAM    
