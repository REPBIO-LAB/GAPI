'''
Module 'virusSR' - for dealing with virus specific needs
'''

## External
import time
import multiprocessing as mp
import itertools
from Bio import SeqIO

## Internal
import bamtools
import log
import sequences
import alignment
from VIGA_SR import bamtools_VIGASR


## [SR CHANGE]
import os
import re

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
    viralSeqs = bamtools_VIGASR.filterBAM2FastaDict(allSam, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus)

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