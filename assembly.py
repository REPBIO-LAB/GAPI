'''
Module 'consensus' - Contains funtions to create a consensus sequence from a set of input sequences
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
import log
import unix
import alignment
import formats

## FUNCTIONS ##
def extract_contigs(GFA, outDir):
    '''
    Extract contigs from Graphical Fragment Assembly (GFA) file

    Input:
        1. GFA: GFA file
        2. outDir: output directory

    Output:
        1. contigs: FASTA file containing contigs
    '''

    # TO DO. Implement function 
    contigs = outDir + '/contigs.fa'
    # awk '$1 ~/S/ {print ">"$2"\n"$3}' GFA > reads.fasta


def assembly_pipeline_miniasm_racon(sequences, technology, outDir):
    '''
    Assembly pipeline based on miniasm for generating contings and racon for contigs polishing. 

    Input:
        1. sequences: FASTA file containing set of sequences to be assembled
        2. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        3. outDir: output directory

    COMMENT: Incorporate nbRounds argument to specify how many rounds of polishing will be performed

    Output:
        1. polished: FASTA file containing polished contigs
    '''
    ## 1. Do assembly with miniasm
    assemblyDir = outDir + '/Assembly/'
    unix.mkdir(assemblyDir)

    contigs = assemble_miniasm(sequences, technology, assemblyDir)

    ## 2. Polish contigs with racon
    polishDir = outDir + '/Polish/'
    unix.mkdir(polishDir)

    polished = polish_racon(contigs, sequences, technology, assemblyDir)

    return polished

def assemble_miniasm(sequences, technology, outDir):
    '''
    Assemble a set of sequences into uncorrected contigs with miniasm

    Input:
        1. sequences: FASTA file containing set of sequences to be assembled
        2. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        3. outDir: output directory

    Output:
        1. contigs: FASTA file containing generated contigs
    '''

    print('INPUT: ', sequences, technology, outDir)

    ## 0. Create directories:
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align sequences All-vs-all 
    ## Set preset according to the technology
    preset = alignment.minimap2_presets(technology)

    print('PRESET: ', preset)
    
    ## Do alignment 
    PAF = outDir + '/overlaps.paf'
    err = open(logDir + '/minimap2.err', 'w') 
    command = 'minimap2 -x ' + preset + ' ' + sequences + ' ' + sequences + ' > ' + PAF

    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ASSEMBLE-MINIASM'
        msg = 'All-vs-all alignment failed' 
        log.step(step, msg)

    ## 2. Assemble reads with miniasm
    ### Relevant options:
    ## Pre-selection
    # -s INT      min span [2000]. Drop mappings shorter than INT -bp. This option also affects the second round of read filtering and minimal overlap length.
    # -c INT      min coverage [3]. Minimal coverage by other reads [3]. In the first round of filtering, miniasm finds the longest region covered by INT or more reads. In the second round, it in addition requires each remaining base to be covered by INT bases at least minSpan /2 from the ends of other reads.
    ## Overlap
    # -o INT      min overlap [same as -s]. Minimal overlap length [same as minSpan ]
    # -h INT      max over hang length [1000]. Maximum overhang length [1000]. An overhang is an unmapped region that should be mapped given a true overlap or true containment. If the overhang is too long, the mapping is considered an internal match and will be ignored.
    # -I FLOAT    min end-to-end match ratio [0.8]. Minimal ratio of mapping length to mapping+overhang length for a mapping considered a containment or an overlap [0.8]. This option has a similar role to -h, except that it controls the ratio, not length.

    GFA = outDir + '/assembled.gfa'
    err = open(logDir + '/racon.err', 'w') 
    command = 'miniasm -f ' + sequences + ' ' + PAF + ' > ' + GFA
    print(command)    

    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ASSEMBLE-MINIASM'
        msg = 'Sequence assembly failed' 
        log.step(step, msg)

    ## 3. Extract contigs from GFA file
    contigs = extract_contigs(GFA, outDir)

    ## Do cleanup

    return contigs
    
def assemble_overlap(fastaA, fastaB, technology, outDir):
    '''
    Assemble two sets sequences via the identification of reciprocal overlap at the sequence ends. 
    (custom algorithm based on minimap2 alignments)

     -: aligned; _: clipped      sequence_bkp
    sequence A:     -------------------*_____________________
    sequence B:                                    |________|______________*-------------------
    assembled:      ---------------------------------------------------------------------------

    Input:
        1. fastaA: FASTA file containing set of sequences A
        2. fastaB: FASTA file containing set of sequences B
        3. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        4. outDir: output directory

    Output:
        1. contigFile: FASTA file containing the generated contig
    '''    
    ## 0. Create directories ##
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align sequences A vs B ##
    ## Set preset according to the technology
    preset = alignment.minimap2_presets(technology)
    
    ## Do alignment 
    PAF_file = outDir + '/overlaps.paf'
    err = open(logDir + '/minimap2.err', 'w') 
    command = 'minimap2 -x ' + preset + ' ' + fastaB + ' ' + fastaA + ' > ' + PAF_file

    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ASSEMBLE-OVERLAP'
        msg = 'Sequences A vs B alignment failed' 
        log.step(step, msg)

    ## 2. Read PAF alignments ##
    PAF = formats.PAF()
    PAF.read(PAF_file)

    # Exit function if no hit found
    if not PAF.lines:
        return

    ## 3. Filter overlaps 
    # Filtering criteria:
    # Pick alignments on the + strand 
    # Pick alignments starting within 250 from sequences ends. 
    filtered = []

    for overlap in PAF.lines:

        # Compute the distance between the corresponding sequence end and the overlap begin or end position
        distA = overlap.qLen - overlap.qEnd 
        distB = overlap.tBeg
    
        # Apply filters
        if (overlap.strand == '+') and (distA <= 250) and (distB <= 250):
            filtered.append(overlap)

    # Replace raw by filtered overlaps
    PAF.lines = filtered

    # Exit function if all the overlapping hits have been filtered
    if not PAF.lines:
        return

    ## 4. Pick longest overlap passing the filters
    # Sort PAF in decreasing overlap lengths
    PAF.lines = PAF.sortByLen()

    # Pick longest overlap
    overlapLongest = PAF.lines[0]

    ## 5. Concatenate overlapping sequences to generate a contig 
    ## Read input fasta files
    sequencesA = formats.FASTA()
    sequencesA.read(fastaA)
    
    sequencesB = formats.FASTA()
    sequencesB.read(fastaB)
    
    ## Contig generation  
    #  -: aligned; _: clipped                     sequence_bkp    qBeg       qEnd
    # sequence A (query):                ---------------*___________|_________|___      sequence_bkp
    # sequence B (template):                                    ____|_________|______________*-------------------
    #                                                             tBeg       tEnd
    # contig:                            -------------------------------------|----------------------------------
    #                                                  sequence_1        :qEnd tEnd:       sequence_2
    seqA = sequencesA.seqDict[overlapLongest.qName][:overlapLongest.qEnd]
    seqB = sequencesB.seqDict[overlapLongest.tName][overlapLongest.tEnd:]
    contigSeq = seqA + seqB 

    ## Write contig into fasta
    contig = formats.FASTA()
    contig.seqDict['CONTIG'] = contigSeq

    contigFile = outDir + '/contig.fa'
    contig.write(contigFile)

    return contigFile

def polish_racon(templates, sequences, technology, outDir):
    '''
    Use a collection of sequences to polish a set of target sequences (i.e. assembled contig) 
    
    Input:
        1. templates: FASTA file containing sequences to be polished
        2. sequences: FASTA file containing set of sequences used to polish the templates
        3. technology: sequencing technology (NANOPORE, PACBIO or ILLUMINA)
        4. outDir: output directory

    COMMENT: Incorporate nbRounds argument to specify how many rounds of polishing will be performed

    Output:
        1. polished: FASTA file containing polished sequences
    '''
    ## 0. Create logs directory:
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align reads against the template ##
    ## Set preset according to the technology
    preset = alignment.minimap2_presets(technology)

    ## Do alignment 
    PAF = outDir + '/alignments.paf'
    err = open(logDir + '/minimap2.err', 'w') 
    command = 'minimap2 -x ' + preset + ' ' + templates + ' ' + sequences + ' > ' + PAF
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'POLISH-RACON'
        msg = 'Alignment of sequences against template failed' 
        log.step(step, msg)

    ## 2. Template polishing with racon ##
    polished = outDir + '/polished.fasta'
    err = open(logDir + '/racon.err', 'w') 
    command = 'racon ' + sequences + ' ' + PAF + ' ' + templates + ' > ' + polished
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'POLISH-RACON'
        msg = 'Template polishing failed' 
        log.step(step, msg)

    ## Do cleanup

    return polished
    