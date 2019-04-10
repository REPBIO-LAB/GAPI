'''
Module 'virus' - for dealing with virus specific needs
'''
## External
import subprocess
import os

## Internal
import bamtools
import formats
import log

def is_virusSR(events, tumourBam, normalBam, eventType, outDir):
    '''
    '''

    ## 1. Collect mate sequence of discordant events ##
    bamtools.collectMatesSeq(events, tumourBam, normalBam, True, 20)

    ## 2. Identify mate sequence of discordant events ##
    eventsIdentityDict = identifySequence(events, eventType, outDir)

    return eventsIdentityDict


def identifySequence(events, eventType, outDir):
    '''
    From a list of events, perform a search with their sequence in order to find out its identity.
    
    Input:
        1. events: list of events
        2. eventType
        3. outDir
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
            fastaDict[event.readId]=event.mateSeq
            fastaObj.seqDict = fastaDict
            FASTA_file = outDir + '/' + str(event.id) + '.fasta'
            fastaObj.write(FASTA_file)

            # Perform aligment
            PAF_file = outDir + '/' + str(event.id) + 'alignments.paf'
            err = open(outDir + '/identifyMate.err', 'w')
            command = 'minimap2 /lustre/scratch117/casm/team154/jt14/3vi/data/databases/RVDBv12.2_MEIGA_HUMAN/consensusViralDb.mmi ' + FASTA_file + ' > ' + PAF_file
            status = subprocess.call(command, stderr=err, shell=True)

            if status != 0:
                step = 'IDENTIFY MATE SEQ'
                msg = 'Identify mate sequence failed' 
                log.step(step, msg)

            # If PAF file is not empty
            if not os.stat(PAF_file).st_size == 0:
                PAFObj = formats.PAF()
                PAFObj.read(PAF_file)

                # Pick the identity of the aligment with highest number of matches
                aligmentMaxNbMatches = PAFObj.sortNbMatches()[0]
                identity = aligmentMaxNbMatches.tName.split('|')[2]

                # Add identity to the eventType and make the output dictionary
                eventTypeIdentity = eventType + '-' + identity

                if eventTypeIdentity not in eventsIdentityDict.keys():
                    print (eventTypeIdentity)
                    eventsIdentityDict[eventTypeIdentity] = []
                eventsIdentityDict[eventTypeIdentity].append(event)
    
    ## TODO: Do cleanup

    return eventsIdentityDict