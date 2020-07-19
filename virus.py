'''
Module 'virus' - for dealing with virus specific needs
'''
## External
import subprocess
import unix
import time
import multiprocessing as mp
import itertools
from Bio import SeqIO

## Internal
import bamtools
import formats
import log
import sequences
import filters
import alignment
import retrotransposons
import callers

## [SR CHANGE]
import os
import re

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

# TODO: UNUSED FUNCTION
def identifySequence(events, outDir, viralDb):
    '''
    From a list of events, perform a search with their sequence in order to find out its identity.
    
    Input:
        1. events: list of events
        2. outDir
    Output:
        1. eventsIdentityDict -> keys: eventType + identity. values: list of events
    '''

    ## TODO: Make specific temp folder for these filer

    eventsIdentityDict = {}
    fastaObj = formats.FASTA()
    
    for event in events:
        if event.mateSeq != None:

            # Write fasta file containing the read sequence
            fastaDict={}
            fastaDict[event.readName]=event.mateSeq
            fastaObj.seqDict = fastaDict
            FASTA_file = outDir + '/' + str(event.id) + '.fasta'
            fastaObj.write(FASTA_file)

            PAF_file = outDir + '/' + str(event.id) + 'alignments.paf'

            #db = '/lustre/scratch117/casm/team154/jt14/3vi/data/databases/RVDBv12.2_MEIGA_HUMAN/consensusViralDb.mmi'

            # Perform aligment
            aligmentMaxNbMatches = sequences.aligmentMaxNbMatches(FASTA_file, viralDb, PAF_file, outDir)

            command = 'rm ' + FASTA_file               
            os.system(command) 

            if aligmentMaxNbMatches != None:

                # Require a minimum percentage of aligment
                lenAlig = abs(aligmentMaxNbMatches.qEnd - aligmentMaxNbMatches.qBeg)

                percAlig = filters.fraction(lenAlig, aligmentMaxNbMatches.qLen)

                # TODO: ponerlo de parametrooo???
                if percAlig:
                    if percAlig > 0.5: 
                    
                        identity = aligmentMaxNbMatches.tName.split('|')[0]
                        # In order to know more than simply the species.
                        #specificIdentity = aligmentMaxNbMatches.tName.split('|')[1]
                        specificIdentity = re.split('\|| ',aligmentMaxNbMatches.tName)[1]

                        # Add identities to event object
                        event.specificIdentity = specificIdentity

                        # Add identity to the eventType and make the output dictionary
                        eventTypeIdentity = event.orientation + '-' + event.type + '-' + identity

                        if eventTypeIdentity not in eventsIdentityDict.keys():
                            eventsIdentityDict[eventTypeIdentity] = []

                        eventsIdentityDict[eventTypeIdentity].append(event)

            ## DESILENCIAAAAR!
            #command = 'rm ' + PAF_file               
            #os.system(command) 
    
    ## TODO: Do cleanup

    return eventsIdentityDict

def init(l):
    global lock
    lock = l

#virus.find_viral_reads(self.bam, self.normalBam, self.confDict['viralDb'], self.confDict['komplexityThreshold'], self.confDict['minTotalMatchVirus'], self.confDict['minParcialMatchVirus'], self.confDict['maxMatchCheckMAPQVirus'], self.confDict['minMAPQVirus'], self.confDict['maxBasePercVirus'], self.confDict['minLccVirus'], self.confDict['processes'], self.outDir)
def find_virus_discordants(bam, normalBam, viralDb, komplexityThreshold, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, processes, outDir):
    
    # Create directory
    collectVirusDir = outDir + '/COLLECT_VIRUS'

    # Filter by complexity (with komplexity)
    allFastas = sequences.komplexityFilter(komplexityThreshold, 'allFastas_all.fasta', 'allFastas.fasta', collectVirusDir)
    
    # Align with bwa allFastas vs viralDb and filter resulting bam
    # TODO SR: bwa allFastas vs viralDb: check if bwa -T parameter does something that we need
    allSam = alignment.alignment_bwa(allFastas, viralDb, 'allSam', processes, collectVirusDir)

    # Index bam
    bamtools.samtools_index_bam(allSam, collectVirusDir)

    # WHEN THIS PART IS DESILENCED THIS IS NEEDED:
    #allSam = collectVirusDir + '/allSam.sam'
    viralSeqs = bamtools.filterBAM2FastaDict(allSam, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus)

    # Collect all identities that have hits in viralBam
    fastaIdentities = list(set(list(itertools.chain(*viralSeqs.values()))))

    # Write a fasta file containing only those sequences that are in identities:
    # self.confDict['viralDb'] -> Papillomaviridae|LC270039.1 02-AUG-2017
    # fastaIdentities = self.viralSeqs.values() -> Papillomaviridae|LC270039.1
    identDbPath = collectVirusDir + '/identDb.fasta'
    identDb = open(identDbPath, 'w')
    viralDb = open(viralDb, 'r')

    for record in SeqIO.parse(viralDb,'fasta'):
        for fastaIdentity in fastaIdentities:
            if fastaIdentity in record.id:
                identDb.write(">" + record.id + "\n")
                identDb.write(str(record.seq) + "\n")
                break

    viralDb.close()
    identDb.close()

    return viralSeqs, identDbPath

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
    chain = PAF.chain(100, 30)

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
    structure['STRAND'] = retrotransposons.infer_strand_alignment(structure['INS_TYPE'], chain)

    ## 4.3 Sequence lengths 
    lengths = infer_lengths(structure['INS_TYPE'], chain, structure['STRAND'])
    structure.update(lengths)


    ## 4.5 Target site duplication (TO DO LATER...)
    #search4tsd()
    
    ## 4.6 Percentage resolved
    structure['PERC_RESOLVED'] = chain.perc_query_covered()
    
    return structure

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
