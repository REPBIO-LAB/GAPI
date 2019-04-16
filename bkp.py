'''
Module to solve the bkp INGUAL DESPUES SE PUEDE UNIR!
'''
import databases
import formats
import os
import subprocess
import log
import unix

def analizeBkp(clustersBinDb, identityDbFasta, reference, side, outDir): # 'RIGHT' OR 'LEFT'
    
    bkpDir = outDir + '/databases'
    dictMetaclusters = {}

    for metacluster in clustersBinDb.collect(['METACLUSTERS']):

        dictMetaclusters[metacluster] = {}

        bkpDir = outDir + '/' + metacluster.ref + '_' + str(metacluster.beg) + '_' + str(metacluster.end)
        unix.mkdir(bkpDir)

        consensusSeq = bkp(metacluster, side, bkpDir)

        # Build db for minimap2
        ref = metacluster.ref
        beg = metacluster.beg
        end = metacluster.end

        # From all specific identitities from all events in the metacluster, pick the specific identity that is repeated more times (si hay dos que coinciden devuelve una al azar.)
        # Coger la identity del primer evento DISCORDANT que aparezca (pq los clipping no tienen identity)
        identity = next(event.identity for event in metacluster.events if event.type == "DISCORDANT")
        

        specificIdentity = max(set([event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"]), key=[event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"].count)

        # Build database con identity + ref
        # TODO coger la identity que ya hemos reconocido, en lugar de toda la db!
        refIdentityDbIndex = databases.buildRefIdentityDb(ref, beg, end, identity, specificIdentity, identityDbFasta, reference, side, bkpDir)

        PAF_file = pafAligment(metacluster, consensusSeq, side, refIdentityDbIndex, bkpDir)

        # if PAF_file is not empty
        if os.path.getsize(PAF_file) > 0:

            PAFObj = formats.PAF()
            PAFObj.read(PAF_file)
            chain = PAFObj.chain()

            dictMetaclusters[metacluster]['refBkpRef'], dictMetaclusters[metacluster]['refBkpPos'], dictMetaclusters[metacluster]['intBkpPos'] = defineBkpPos(chain, side)
            dictMetaclusters[metacluster]['seq'] = defineBkpSeq(chain, consensusSeq, side)

    return dictMetaclusters


def pafAligment(metacluster, consensusSeq, side, refIdentityDbIndex, outDir):

        # Create consensus fasta
        consensusSeqName = side + 'ConsensusSeq'
        consensusFastaObj = formats.FASTA()
        consensusFastaDict = {}
        consensusFastaDict[consensusSeqName] = consensusSeq

        consensusFastaObj.seqDict = consensusFastaDict

        consensusFastaPath = outDir + '/' + str(metacluster.beg)+ '_' + side + '_consensusFasta.fa'

        consensusFastaObj.write(consensusFastaPath)

        # Aling with minimap2
        # TODO ponerlo bien
        PAF_file = outDir + '/' + str(metacluster.beg)+ '_' + side + '_alignments.paf'
        #err = open(logDir + '/align.err', 'w') 
        command = 'minimap2 ' + refIdentityDbIndex + ' ' + consensusFastaPath + ' > ' + PAF_file
        status = subprocess.call(command, shell=True)

        if status != 0:
            step = 'ALIGN-INSERT'
            msg = 'Insert alignment failed' 
            log.step(step, msg)

        return PAF_file



def defineBkpPos(chain, side):
    refBkpRef = None
    refBkpPos = None
    intBkpPos = None
    for alig in chain.alignments:
        qName = side + 'ConsensusSeq'
        if alig.qName == qName: 
            if not 'consensus' in alig.tName:
                refBkpRef = alig.tName.split(':')[0]
                refBkpPos = int(alig.tName.split('-')[0].split(':')[1]) + alig.tBeg
            else:
                intBkpPos = alig.tBeg
    
    return refBkpRef, refBkpPos, intBkpPos

def defineBkpSeq(chain, ConsensusSeq, side):
    seq = ''
    refBeg = 0
    refEnd = 0
    intBeg = 0
    intEnd = 0
    intName = ''


    # Pick begs and ends and int name
    for alig in chain.alignments:
        qName = side + 'ConsensusSeq'
        if alig.qName == qName:
            if not 'consensus' in alig.tName:
                refBeg = alig.qBeg
                refEnd = alig.qEnd
            else:
                intBeg = alig.qBeg
                intEnd = alig.qEnd
                intName = alig.tName

    # Construct sequence.
    if side == 'RIGHT':
        # In case the aligment is not continuos, print the number of bases that it lacks.
        if intBeg != 0 and refEnd != 0 and (intBeg - refEnd) > 1:
            arrow = '(' + str(intBeg - refEnd) + ')]>'
        else:
            arrow = ']>'
        seq = ConsensusSeq[refBeg:refEnd] + '[' + intName + arrow + ConsensusSeq[intBeg:intEnd]

    elif side == 'LEFT':
        if refBeg != 0 and intEnd != 0 and (refBeg - intEnd) > 1:
            arrow = '<[(' + str(refBeg - intEnd) + ')'
        else:
            arrow = '<['
        seq = ConsensusSeq[intBeg:intEnd] + arrow + intName + ']' + ConsensusSeq[refBeg:refEnd]
    
    return seq


def bkp(metacluster, side, outDir):

    consensusSeq = None
    fastaObj = formats.FASTA()
    fastaDict = {}

    # Hago un dictionario con los diferentes subclusters que hay en el metacluster
    subclustersDict = metacluster.subclusters
    # Cojo solo el subcluster CLIPPING que se corresponde con el lado elegido (left or right)
    # PQ LAS SECUENCIAS DE LOS DISCORDANT NO LAS QUIERO, NO?
    subclustersType = [key for key in subclustersDict.keys() if side in key]

    # Si hay clusters del tipo que buscamos
    if len(subclustersType) > 0:

        # Meto todos los eventos de los subclusters elegidos en un dictionario de fasta
        for event in subclustersDict[subclustersType[0]].events:
            fastaDict[event.readName] = event.readSeq


        fastaObj.seqDict = fastaDict

        fastaPath = outDir + '/' + str(metacluster.beg) + '_' + side +'_supportingReads.fa'
        fastaObj.write(fastaPath)


        ### 2. Make multiple sequence alignment
        msfPath = outDir + '/' + str(metacluster.beg) + '_' + side +'_supportingReads.msf'
        command = 'muscle -in ' + fastaPath + ' -out ' + msfPath + ' -msf' 
        status = subprocess.call(command, shell=True)

        ### 3. Generate consensus sequence (cons tool from EMBOSS packagge)
        consensusPath = outDir + '/' + str(metacluster.beg)  + '_' + side +'_consensus.fa'

        command = 'cons -sequence ' + msfPath + ' -outseq ' + consensusPath + ' -identity 0 -plurality 0'
        status = subprocess.call(command, shell=True)

        ### Read consensus sequence 
        consensusFastaObj = formats.FASTA()
        consensusFastaObj.read(consensusPath)
        consensusSeq = consensusFastaObj.seqDict["EMBOSS_001"].upper()

        # TODO
        ### Do cleanup
        #command = 'rm ' + fastaPath + ' ' + msfPath + ' ' + consensusPath             
        #os.system(command) # returns the exit status

        ## Replace '-' by 'N' for ambiguous bases:
        consensusSeq = consensusSeq.replace('-', 'N')

        ## Convert consensus sequence into upper case:
        consensusSeq = consensusSeq.upper()

    return consensusSeq