'''
Module 'virus' - for dealing with virus specific needs
'''
## External
import unix
import pysam
from cigar import Cigar
import numpy as np
import os
import Bio.SeqUtils
from Bio.SeqUtils import lcc
import multiprocessing as mp

## Internal
import formats
import log
import sequences
import alignment

# Multiprocessing lock as global variable. Useful for safely writting from different processes to the same output file.
def init(l):
    global lock
    lock = l

# TODO: This function are quite similar to retrotransposons ones, so maybe we can merge them in some way
def insertion_type(chain):
    '''
    Scan alignments chain to determine the type of insertion

    Input:
        1. chain: Sequence chain of alignments over viral consensus sequences
        
    Output:
        1. insType: Insertion type (viralSolo (only one viral sequence is inserted), viralFamilyNested (different viruses from the same family are inserted), viralNested (different viruses from different families are inserted) or None)
        2. family: List of viral families
        3. srcId: accession number of viral consensus sequence.
    ''' 

    # TODO: Add virus description to virus features
    ## Make list containing all the different families the sequence aligns into
    families = list(set([alignment.tName.split("|")[0] for alignment in chain.alignments]))
    nbFamilies = len(families)

    # TODO: this is a temp fix for unclassified viruses as it has an space in the name, but it should be fixed in another way
    try:
        ## Make list containing the accesion numbers the sequence aligns into 
        accNbs = list(set([alignment.tName.split("|")[1] for alignment in chain.alignments]))
        nbaccNbs = len(accNbs)
    except IndexError:
        accNbs = []
        nbaccNbs = len(accNbs)
    
    ## Make list containing the sequence description the sequence aligns into
    try:
        virusDescription = list(set([alignment.tName.split("|")[2] for alignment in chain.alignments]))
    except IndexError:
        virusDescription = []

    ## a) viralSolo insertion (only one viral sequence is inserted)
    if (nbFamilies == 1) and (nbaccNbs == 1):
        insType = 'viralSolo'
        family = families
        srcId = accNbs

    ## b) viralFamilyNested (different viruses from the same family are inserted)
    elif (nbFamilies == 1) and (nbaccNbs > 1):
        insType = 'viralFamilyNested'
        family = families
        srcId = accNbs
        
    ## c) viralNested (different viruses from different families are inserted)
    elif (nbFamilies > 1):
        insType = 'viralNested'
        family = families
        srcId = accNbs

    ## d) Unknown insertion type
    else:
        insType = 'unknown'
        family = []
        srcId = [] 

    return insType, family, srcId
    
def infer_lengths(insType, chain, strand):
    '''
    Determine the length of each type of sequence composing the insertion
    
    Input:
        1. insType: Insertion type (solo, nested, orphan, partnered or None)
        2. chain: Sequence chain of alignments over viral consensus sequences
        3. strand: Insertion strand (+ or -)

    Output:
        1. lengths: dictionary containing length information
    ''' 
    ### Initialize dictionary
    lengths = {}

    ### 1. Compute the length of each type of sequence composing the insertion

    ## Pick only those hits over viral consensus sequence
    viralHits = [hit for hit in chain.alignments]

    for hitNb, viralHit in enumerate(viralHits):

        for feature in ['VIRAL_COORD', 'VIRAL_LEN', 'IS_FULL', 'TRUNCATION_5_LEN', 'TRUNCATION_3_LEN']:
            lengths['viralHit_' + str(hitNb)] = {}
            lengths['viralHit_' + str(hitNb)][feature] = None
        
        ## Determine piece of consensus sequence that has been integrated
        #ref = viralHit.tName.split('|')[1]
        viralBeg = viralHit.tBeg
        viralEnd = viralHit.tEnd

        ## Compute length
        lengths['viralHit_' + str(hitNb)]['VIRAL_LEN'] = viralEnd - viralBeg

        ## Assess if full length viral insertion
        consensusLen = viralHit.tLen 
        percConsensus = float(lengths['viralHit_' + str(hitNb)]['VIRAL_LEN']) / consensusLen * 100
        lengths['viralHit_' + str(hitNb)]['IS_FULL'] = True if percConsensus >= 95 else False

        ## Compute truncation length at both ends
        lengths['viralHit_' + str(hitNb)]['TRUNCATION_5_LEN'] = viralBeg   
        lengths['viralHit_' + str(hitNb)]['TRUNCATION_3_LEN'] = consensusLen - viralEnd
    
    return lengths

# TODO: This function are quite similar to retrotransposons ones, so maybe we can merge them in some way
def virus_structure(FASTA_file, index, outDir):
    '''    
    Infer the insertion size, structure and strand of viral insertions

    Input:
        1. FASTA_file: Path to FASTA file containing the sequence
        2. index: Minimap2 index for consensus viral sequences database
        3. outDir: Output directory
        
    Output:
        1. structure: dictionary containing insertion structure information
    '''     
    structure = {}

    ## 0. Create logs directory ##
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Align the sequence into the viral sequences database ##
    PAF_file = alignment.alignment_minimap2(FASTA_file, index, 'alignment2consensusVirus', 1, outDir)

    ## 2. Read PAF alignments ##
    PAF = formats.PAF()
    PAF.read(PAF_file)

    # Exit function if no hit on the viral database
    if not PAF.alignments:
        return structure

    ## 3. Chain complementary alignments ##
    # TODO; De momento dejo 50 pero ajustarlo!
    chain = PAF.chain(100, 50)
    for ali in chain.alignments:
        print ('EVAH')
        print (ali.qBeg)
        print (ali.qEnd)
    print (chain.perc_query_covered())
    print (PAF_file)

    ## 4. Infer insertion features ##
    ## Retrieve inserted seq
    #FASTA = formats.FASTA()
    #FASTA.read(FASTA_file)
    #sequence = list(FASTA.seqDict.values())[0]

    # TODO: Maybe it's a good idea to change the junciton.consSeq for this one
    
    ## 4.1 Insertion type
    # TODO: Put VIRUSDSC in output
    structure['INS_TYPE'], structure['FAMILY'], structure['CYTOBAND'] = insertion_type(chain)

    ## 4.2 Insertion strand
    #structure['STRAND'] = retrotransposons.infer_strand_alignment(structure['INS_TYPE'], chain)
    structure['STRAND'] = None

    ## 4.3 Sequence lengths 
    lengths = infer_lengths(structure['INS_TYPE'], chain, structure['STRAND'])
    structure.update(lengths)


    ## 4.5 Target site duplication (TO DO LATER...)
    #search4tsd()
    
    ## 4.6 Percentage resolved
    structure['PERC_RESOLVED'] = chain.perc_query_covered()
    
    return structure


def is_virusSR(events, viralSeqs):
    '''
    TODO SR: NOW THIS IS SAME AS IN RT, SO I SHOULD DO ONLY 1!!
    '''
    matesIdentity = {}

    for discordant in events:

        if discordant.readName in viralSeqs.keys():
            identity = viralSeqs[discordant.readName]
        else:
            identity = None
        
        discordant.setIdentity(identity)

        if discordant.identity:
            featureType = discordant.identity
        else:
            featureType = 'None'

        ## Add discordant read pair to the dictionary
        #identity = discordant.orientation + '-DISCORDANT-' + featureType
        identity = 'DISCORDANT-' + featureType


        if featureType != 'None':
            discordant.element = 'VIRUS'

        # a) There are already discordant read pairs with this identity
        if identity in matesIdentity:
            matesIdentity[identity].append(discordant)

        # b) First discordant read pair with this identity
        else:
            matesIdentity[identity] = [ discordant ] 
    
    return matesIdentity
    '''
    ## 1. Collect mate sequence of discordant events ##
    msg = '[Start bamtools.collectMatesSeq]'
    log.subHeader(msg)
    start = time.time()
    bamtools.collectMatesSeq(events, tumourBam, normalBam, True, 20)
    end = time.time()
    print("TIEMPO DE collectMatesSeq" + str(end - start))
    msg = '[End bamtools.collectMatesSeq]'
    log.subHeader(msg)

    ## 2. Identify mate sequence of discordant events ##
    msg = '[Start identifySequence]'
    log.subHeader(msg)
    eventsIdentityDict = identifySequence(events, outDir, viralDb)
    msg = '[End identifySequence]'
    log.subHeader(msg)

    return eventsIdentityDict
    '''

def filterBAM2FastaDict(BAM, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='SR'):
    '''
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}

    # For each read alignment
    for alignmentObj in iterator:
        alignmentPass = False
        #numMatches = 0
        queryCoord = 0
        print ('alignmentObj.query_name ' + str(alignmentObj.query_name))
        #print ('alignmentObj.is_unmapped ' + str(alignmentObj.is_unmapped))

        if not alignmentObj.is_unmapped:
            ctuples = alignmentObj.cigartuples
            allMatches = [t[1] for t in ctuples if t[0] == 0]
            totalMatch = sum (allMatches)
            print ('totalMatch ' + str(totalMatch))
            if totalMatch >= minTotalMatchVirus:
                c = Cigar(alignmentObj.cigarstring)
                for citem  in list(c.items()):
                    # If cigar is query consuming, update query coordinates:
                    if citem[1] != 'M' and citem[1] != 'D' and citem[1] != '=':
                        queryCoord = queryCoord + int(citem[0])
                    elif citem[1] == 'M' or citem[1] == '=':
                        print ('citem[0] ' + str(citem[0]))
                        print ('maxMatchCheckMAPQVirus'  + str(maxMatchCheckMAPQVirus))
                        print ('alignmentObj.mapping_quality'  + str(alignmentObj.mapping_quality))
                        print ('minMAPQVirus'  + str(minMAPQVirus))
                        if citem[0] >= minParcialMatchVirus:
                            if (citem[0] <= maxMatchCheckMAPQVirus and alignmentObj.mapping_quality > minMAPQVirus) or citem[0] > maxMatchCheckMAPQVirus:
                                sequence = alignmentObj.query_sequence[queryCoord:(queryCoord + int(citem[0]))]
                                # Calculate base percentage
                                basePercs = sequences.baseComposition(sequence)[1]
                                # Delete total value of base percentage result
                                del basePercs['total']
                                print ('basePercs ' + str(basePercs))
                                # Only those sequences with base percentage lower than 85 are collected:
                                if all(perc < maxBasePercVirus for perc in basePercs.values()):
                                    # Calculate complexity
                                    complexity = Bio.SeqUtils.lcc.lcc_simp(sequence)
                                    print ('complexity ' + str(complexity))
                                    print ('minLccVirus ' + str(minLccVirus))
                                    if complexity > minLccVirus:
                                        #print ('PassAll')
                                        # TODO: put this in a good way!!!!!
                                        if mode=='LR':
                                            if len(sequence)>10:
                                                alignmentPass = True
                                                break # Cambio 28/08 11:53
                                        else:
                                            alignmentPass = True
                                            break # Cambio 28/08 11:53
                                    else:
                                        queryCoord = queryCoord + int(citem[0])
                                else:
                                    queryCoord = queryCoord + int(citem[0])
                        else:
                            queryCoord = queryCoord + int(citem[0])
        '''
        print (numMatches)
        if numMatches > 90: # TODO SR: put this as an option
            alignmentPass = True
        elif alignmentObj.mapping_quality >= viralBamMAPQ and numMatches >= int(selectedPartialMatches):
            alignmentPass = True
        # New condition:
        # TODO SR: Put as options!!!
        elif alignmentObj.mapping_quality >= 40 and numMatches >= 60:
            alignmentPass = True
        '''
        print ('alignmentPass ' + str(alignmentPass))
        if alignmentPass == True:
            # Add to fasta dict
            if alignmentObj.query_name in fastaDict.keys():
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
            else:
                fastaDict[alignmentObj.query_name] = []
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)

    return fastaDict