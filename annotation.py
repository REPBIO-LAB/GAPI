'''
Module 'annotation' - Contains functions for the annotation of genomic intervals according to different annotation resources
'''

## DEPENDENCIES ##
# External
import os
import subprocess
from operator import itemgetter

# Internal
import unix
import formats
import databases
import log

def load_annotations(annotations2load, refLengths, annotationsDir, threads, outDir):
    '''
    Load a set of annotation files in bed formats into a bin database

    Input:
        1. annotations2load: list of annotations to load. Annotations available: REPEATS, TRANSDUCTIONS and EXONS
        2. refLengths: Dictionary containing reference ids as keys and as values the length for each reference  
        3. annotationsDir: Directory containing annotation files
        4. threads: number of threads used to parallelize the bin database creation
        5. outDir: Output directory
    
    Output:
        1. annotations: directory containing one key per type of annotation loaded and bin databases containing annotated features as values (None for those annotations not loaded)
    '''

    ## 0. Initialize dictionary
    annotations = {}
    annotations['REPEATS'] = None
    annotations['TRANSDUCTIONS'] = None
    annotations['EXONS'] = None

    ## 1. Create output directory
    unix.mkdir(outDir)

    ## 2. Load annotated repeats into a bin database
    if 'REPEATS' in annotations2load:
        log.info('2. Load annotated repeats into a bin database')
        #repeatsBed = annotationsDir + '/repeats_repeatMasker.bed'
        repeatsBed = annotationsDir + '/repeats_repeatMasker.L1.bed'
        annotations['REPEATS'] = formats.bed2binDb(repeatsBed, refLengths, threads)

    ## 3. Create transduced regions database
    if 'TRANSDUCTIONS' in annotations2load:
        log.info('3. Create transduced regions database')

        ## Create bed file containing transduced regions
        sourceBed = annotationsDir + '/srcElements.bed'
        transducedPath = databases.create_transduced_bed(sourceBed, 15000, outDir)

        ## Load transduced regions into a bin database
        annotations['TRANSDUCTIONS'] = formats.bed2binDb(transducedPath, refLengths, threads)

    ## 4. Create exons database
    if 'EXONS' in annotations2load:

        log.info('4. Create exons database')
        exonsBed = annotationsDir + '/exons.bed'
        annotations['EXONS'] = formats.bed2binDb(exonsBed, refLengths, threads)

    return annotations

def annotate_interval(ref, beg, end, annotDb):
    '''
    Intersect input interval (ref:beg-end) with a given annotation 

    Input: 
        1. ref: reference id
        2. beg: begin position
        3. end: end position
        4. annotDb: dictionary containing annotated features organized per chromosome (keys) into genomic bins (values)

    Output:
        1. sortedOverlaps. List of lists sorted in decreasing percentage of overlap. Each tuple corresponds to one overlapping event and is composed by 3 elements: 
            1. Overlapping event
            2. Number of overlapping base pairs
            3. Percentage of base pairs of the input interval that are overlapping  
            4. Tuple with input interval coordinates overlapping with the event
    '''

    # a) Annotated features available in the same ref 
    if ref in annotDb:
            
        ## Select features bin database for the corresponding reference 
        binDb = annotDb[ref]        

        ## Retrieve all the annotated features overlapping with the input interval
        overlaps = binDb.collect_interval(beg, end, 'ALL')    

        ## Order overlapping features in decreasing order of perc of overlap
        sortedOverlaps = sorted(overlaps, key=lambda x: x[2], reverse=True)
         
    # b) No feature in the same ref as the interval
    else:
        sortedOverlaps = []

    return sortedOverlaps
        

def repeats_annotation(events, repeatsDb, buffer):
    '''
    For each input event assess if overlaps with an annotated repeat in the reference genome

    Input: 
        1. events: list containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. repeatsDb: dictionary containing annotated repeats organized per chromosome (keys) into genomic bins (values)
        3. buffer: number of base pairs to extend begin and end coordinates for each event prior assessing overlap

    Output:
        New 'repeatAnnot' attribute set for each input event. 
        'repeatAnnot' is a tuple containing two lists of overlapping repeat's families and subfamilies 
    '''

    ## Assess for each input event if it overlaps with an annotated repeat
    for event in events:

        # A) Annotated repeat in the same ref where the event is located
        if event.ref in repeatsDb:
            
            ## Select repeats bin database for the corresponding reference 
            repeatsBinDb = repeatsDb[event.ref]        

            ## Retrieve all the annotated repeats overlapping with the event interval
            overlaps = repeatsBinDb.collect_interval(event.beg - buffer, event.end + buffer, 'ALL')    

            ## Make list of overlapping repeats
            families = [overlap[0].optional['family'] for overlap in overlaps]

            ## Collect as well subfamily info
            subfamilies = [overlap[0].optional['subfamily'] for overlap in overlaps]

        # B) No repeat in the same ref as the event
        else:
            families = []
            subfamilies = []
        
        ## Add repeat annotation as attribute 
        # A) Event overlapping repeat
        if families:
            event.repeatAnnot = (families, subfamilies)

        # B) Event NOT overlapping repeat
        else:    
            event.repeatAnnot = (None, None)


def gene_annotation(events, annovarDir, outDir):
    '''
    Perform gene-based annotation for a list of input events
 
    Input: 
        1. events: List containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. annovarDir: Directory containing the two files used by ANNOVAR to perform gene based-annotation:
                a) build_annot.txt     - Text file containing annotated transcript coordinates
                b) build_annotMrna.fa  - Fasta containing annotated transcript sequences
        3. outDir: Output directory

    Output:
    
        New 'geneAnnot' attribute set for each input event. 
        'geneAnnot' is a tuple(region,gene) 
    '''

    ## 1. Create output directory
    unix.mkdir(outDir)

    ## 2. Create input file containing events intervals for ANNOVAR 
    create_annovar_input(events, 'events.annovar', outDir)
    annovarInput = outDir + '/events.annovar'

    ## 3. Annotate events intervals with ANNOVAR 
    out1, out2 = run_annovar(annovarInput, annovarDir, outDir)

    ## 4. Add gene annotation info to the events
    # Read annovar output file into a dict
    out1Dict = read_annovar_out1(out1)

    # Add to each event gene annotation info
    for event in events:

        ## Add gene annotation info
        name = event.ref + ':' + str(event.beg) + '-' + str(event.end)
        event.geneAnnot = out1Dict[name]

    ## Do cleanup
    unix.rm([annovarInput, out1, out2])


def create_annovar_input(events, fileName, outDir):
    '''
    Write events intervals into a format compatible with annovar 

    Input:
        1. events: List containing input events. Events should be objects containing ref, beg and end attributes.
        2. fileName: Output file name 
        3. outDir: Output directory
    '''
    ## 1. Initialize BED 
    BED = formats.BED()

    ## 2. Add events intervals to the BED  
    BED.lines = []
    BED.structure = 'List'

    # For each event
    for event in events:

        # Create BED entry:
        header = ['ref', 'beg', 'end', 'name']
        name = event.ref + ':' + str(event.beg) + '-' + str(event.end)
        fields = [event.ref, event.beg, event.end, name]
        entry = formats.BED_line(fields, header)
        
        # Add entry to BED
        BED.lines.append(entry)

    ## 3. Write bed into a format compatible with annovar 
    outFile = outDir + '/' + fileName
    BED.write_annovar(outFile)


def read_annovar_out1(out1):
    '''
    Read annovar out1 and store info in a dictionary

    Input:
        1. out1: Annovar output file 1 (region annotation for all the variants) 

    Output:
        1. out1Dict: Dictionary with comment field as keys and tuples(region,gene) as values        
    '''

    ## Open annovar out file
    with open(out1, "r") as out1File:
    
        out1Dict = {}

        # Read line by line adding the relevant info to the dict in each iteration
        for line in out1File:
            fields = line.split()
            region = fields[0]
            gene = fields[1]
            name = fields[8] 
            
            out1Dict[name] = (region, gene)

    return out1Dict
            

def run_annovar(inputFile, annovarDir, outDir):
    '''
    For each input entry perform gene-based annotation with Annovar

    Input: 
        1. inputFile: Annovar-like file containing regions to be annotated
        2. annovarDir: Directory where annovar annotation files are located
        3. outDir: Output directory

    Output: 
        1. out1: Annovar output file 1 (region annotation for all the variants) 
        2. out2: Annovar output file 2 (amino acid changes as a result of the exonic variant) 

    NOTE: perl annotate_variation.pl call should be available as 'ANNOVAR' environmental variable
    NOTE: check annovar documentation for explanation about how to interpret output files:

    http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-1-refseq-gene-annotation
    http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-2-refseq-gene-annotation
    '''

    ## Run annovar in the input file
    ANNOVAR = os.environ['ANNOVAR']
    out = open(outDir + '/annovar.out', 'w') 
    err = open(outDir + '/annovar.err', 'w') 
    command = ANNOVAR + ' -buildver build -out ' + outDir + '/annovar -dbtype annot ' + inputFile + ' ' + annovarDir
    status = subprocess.call(command, stdout=out, stderr=err, shell=True)

    ## Return results
    out1 = outDir + '/annovar.variant_function'
    out2 = outDir + '/annovar.exonic_variant_function'

    return out1, out2


def intersect_mate_annotation(discordants, annotation):
    '''
    For each input read assess if the mate aligns in a retrotransposon located elsewhere in the reference genome

    Input: 
        1) discordants: list containing input discordant read pair events
        2) annotation: dictionary containing annotated features organized per chromosome (keys) into genomic bins (values)

    Output:
        1) matesIdentity: dictionary containing lists of discordant read pairs organized taking into account their orientation and if the mate aligns in an annotated feature 
                               This info is encoded in the dictionary keys as follows. Keys composed by 3 elements separated by '-':
                                
                                    - Orientation: read orientation (PLUS or MINUS)
                                    - Event type: DISCORDANT   
                                    - featureType: Feature type or 'None' if mate does not align in a retrotransposon
    '''
    matesIdentity = {}

    ##  For each input discordant intersect mate alignment coordinates with the provided annotation 
    for discordant in discordants:
        
        # A) Annotated repeat in the same ref where the mate aligns
        if discordant.mateRef in annotation:

            ## Select features bin database for the corresponding reference 
            featureBinDb = annotation[discordant.mateRef]        
 
            ## Retrieve all the annotated features overlapping with the mate alignment interval
            overlappingFeatures = featureBinDb.collect_interval(discordant.mateStart, discordant.mateStart + 100, 'ALL')    

            ## Determine mate status 
            # a) Mate does not align within an annotated repeat
            if len(overlappingFeatures) == 0:
                featureType = 'None'

            # b) Mate aligns within a single repeat
            elif len(overlappingFeatures) == 1:
                featureType = overlappingFeatures[0][0].optional['family']

            # c) Mate overlaps multiple repeats
            else:
                overlappingFeatures = sorted(overlappingFeatures,key=itemgetter(1), reverse=True)
                featureType = overlappingFeatures[0][0].optional['family']

        # B) No repeat in the same ref as mate
        else:
            featureType = 'None'

        ## Set discordant identity
        discordant.identity = featureType

        ## Add discordant read pair to the dictionary
        identity = discordant.side + '_DISCORDANT_' + featureType

        # a) There are already discordant read pairs with this identity
        if identity in matesIdentity:
            matesIdentity[identity].append(discordant)

        # b) First discordant read pair with this identity
        else:
            matesIdentity[identity] = [ discordant ] 
    
    return matesIdentity


