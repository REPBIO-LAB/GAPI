'''
Module 'alignment' - Contains funtions align sequences
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
import log
import unix
import bamtools 

## FUNCTIONS ##

def minimap2_presets(technology):
    '''
    Set minimap2 preset according to the sequencing technology

    Input:
        1. technology: Sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        
    Output:
        1. preset: proper preset  
    '''

    if technology == 'PACBIO':
        preset = 'map-pb'
        
    elif technology == 'NANOPORE':
        preset = 'map-ont'

    elif technology == 'ILLUMINA':
        preset = 'sr'

    else:
        log.info('ERROR: ' + technology + ' technology not supported')
        sys.exit(1)

    return preset

def alignment_minimap2(FASTA, index, outDir):
    '''
    Align a set of sequence into a reference with minimap2

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. index: Path to the the index of the reference in .mmi format (generated with minimap2)
        3. outDir: Output directory

    Output:
        1. BAM: Path to sorted BAM file containing input sequences alignments or 'None' if alignment failed 
    '''

    ## 1. Align the sequences into the reference
    # Use -Y to get soft clippings for supplementary alignments
    SAM = outDir + '/alignments.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'minimap2 -Y -a -k15 -w10 ' + index + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)
        
    ## 2. Convert SAM to sorted BAM
    BAM = bamtools.SAM2BAM(SAM, outDir)

    ## 3. Do cleanup
    unix.rm([SAM])

    return BAM


def targeted_alignment_minimap2(FASTA, targetInterval, reference, outDir):
    '''
    Align a set of sequences into a reference target region. 
    
    Useful for doing local realignment of reads around SV breakpoints. Much faster than whole genome realignment

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. targetInterval: Reference genome interval where sequences will be aligned. The interval must be provided as chr:beg-end.
        3. reference: Path to the reference sequences in fasta format. An index of the reference generated with samtools faidx must be located in the same directory
        4. outDir: Output directory
        
    Output:
        1. BAM: Path to sorted BAM file containing input sequences alignments or 'None' if realignment failed 
    '''
    ## 0. Create logs directory
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Extract the reference target region prior alignment 
    target = outDir + '/target.fa'
    err = open(logDir + '/target.err', 'w') 
    command = 'samtools faidx ' + reference + ' ' + targetInterval + ' > ' + target
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'TARGET'
        msg = 'Extraction of reference target region failed' 
        log.step(step, msg)
        return None

    ## 2. Align the sequences into the target region 
    # Use -Y to get soft clippings for supplementary alignments
    SAM = outDir + '/alignments.sam'
    err = open(logDir + '/align.err', 'w') 
    command = 'minimap2 -Y -a -k15 -w10 ' + target + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)
        return None

    ## 3. Convert SAM to sorted BAM
    BAM = bamtools.SAM2BAM(SAM, outDir)

    ## 4. Do cleanup
    unix.rm([target, SAM])

    return BAM    

def targetered2genomic_coord(event, ref, offset):
    '''
    Convert event coordinates resulting from the realignment of a sequence into a target region into genomic coordinates

    Input:
        1. event: INS, DEL or CLIPPING event object
        2. ref: reference corresponding to the targetered seq
        3. offset: offset to be added to event coordinates

    Output:
        1. event: modified event object
    '''
    ## 1. Convert reference                
    event.ref = ref

    ## 2. Add offset to event begin and end coordinates
    event.beg = event.beg + offset
    event.end = event.end + offset 

    ## 3. Add offset to supplementary alignments 
    if event.supplAlignment != None:

        supplAlignments = ''

        # For each supplementary alignment
        for supplAlignment in event.supplAlignment.split(';'):

            # Stop after the last supplementary alignment
            if supplAlignment == '':
                break

            supplRef, supplBeg, strand, cigar, MAPQ, NM = supplAlignment.split(',')

            supplRef = ref
            beg = int(supplBeg) + offset

            supplAlignment = supplRef + ',' + str(beg) + ',' + strand + ',' + cigar + ',' + MAPQ + ',' + NM
            supplAlignments = supplAlignments + supplAlignment + ';' if supplAlignments != '' else supplAlignment + ';'

            event.supplAlignment = supplAlignments

    return event