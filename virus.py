'''
Module 'virus' - for dealing with virus specific needs
'''
## External
import subprocess
import time

## Internal
import bamtools
import formats
import log
import sequences
import filters

## [SR CHANGE]
import os
import re

def is_virusSR(events):
    '''
    TODO: NOW THIS IS SAME AS IN RT, SO I SHOULD DO ONLY 1!!
    '''
    matesIdentity = {}

    for discordant in events:

        if discordant.identity:
            featureType = discordant.identity
        else:
            featureType = 'None'

        ## Add discordant read pair to the dictionary
        identity = discordant.orientation + '-DISCORDANT-' + featureType

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