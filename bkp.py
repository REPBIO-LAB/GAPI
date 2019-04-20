'''
Module to solve the bkp INGUAL DESPUES SE PUEDE UNIR!
'''
import databases
import formats
import os
import subprocess
import log
import unix
## [SR CHANGE]
import sequences
import assembly


def analyzeMetaclusters(clustersBinDb, confDict, bam, normalBam, mode, db, indexDb, outDir):
    '''
    1. Por cada evento en la BinDB (en este caso especificamente por cada metacluster)
    a. Anadir supporting clipping al metacluster de discordant.
    b. Buscar el bkp (punto mas apoyado por los clippings que acabamos de anadir) y quitar aquellos clipping events que no lo soporten
    c. Hacer las cadenas de secuencias para ambos bkps.
    '''

    dictMetaclusters = {}

    for metacluster in clustersBinDb.collect(['METACLUSTERS']):

        leftIntConsensusSeq = None
        rightIntConsensusSeq = None

        # Set origin
        if mode == 'PAIRED':
            metacluster.setIntOrigin()

        dictMetaclusters[metacluster] = {}

        bkpDir = outDir + '/' + metacluster.ref + '_' + str(metacluster.beg) + '_' + str(metacluster.end)
        unix.mkdir(bkpDir)

        # a. Anadir supporting clipping al metacluster de discordant.
        CLIPPING_cluster = metacluster.supportingCLIPPING(100, confDict, bam, normalBam, mode)

        # b. Buscar el bkp (punto mas apoyado por los clippings que acabamos de anadir) y quitar aquellos clipping events que no lo soporten
        dictMetaclusters[metacluster]['refLeftBkp'], dictMetaclusters[metacluster]['refRightBkp'] = clippingBkp(CLIPPING_cluster)

        # c. Hacer las cadenas de secuencias para ambos bkps.
        #dictMetaclusters[metacluster]['leftSeq'], dictMetaclusters[metacluster]['rightSeq'] = makeConsSeqs(CLIPPING_cluster, 'REF', db, indexDb, bkpDir)
        leftRefConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'left', 'REF', db, indexDb, bkpDir)[1]
        rightRefConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'right', 'REF', db, indexDb, bkpDir)[1]

        leftIntConsensusPath, leftIntConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'left', 'INT', db, indexDb, bkpDir)
        rightIntConsensusPath, rightIntConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'right', 'INT', db, indexDb, bkpDir)

        if leftIntConsensusSeq != None:
            leftSeq = leftIntConsensusSeq + '<[INT]' + leftRefConsensusSeq
        else:
            leftSeq = None
        if rightIntConsensusSeq != None:
            rightSeq = rightRefConsensusSeq + '[INT]>' + rightIntConsensusSeq
        else:
            rightSeq = None

        dictMetaclusters[metacluster]['leftSeq'], dictMetaclusters[metacluster]['rightSeq'] = leftSeq, rightSeq

        if leftIntConsensusPath != None:
            dictMetaclusters[metacluster]['intLeftBkp'] =  bkpINT(metacluster, leftIntConsensusPath, db, bkpDir)
        if rightIntConsensusPath != None:
            dictMetaclusters[metacluster]['intRightBkp'] = bkpINT(metacluster, rightIntConsensusPath, db, bkpDir)

        ### Do cleanup
        unix.rm([bkpDir])

        return dictMetaclusters

def clippingBkp(CLIPPING_cluster):
    '''
    Buscar el bkp (punto mas apoyado por los clippings ) y quitar aquellos clipping events que no lo soporten
    '''
    leftBkp = None
    rightBkp = None
    leftBkps = []
    rightBkps = []
    for event in CLIPPING_cluster.events:
        if event.clippedSide == 'left':
            leftBkps.append(event.beg)
        elif event.clippedSide == 'right':
            rightBkps.append(event.beg)

    if len(leftBkps) > 0:
        leftBkp = max(set(leftBkps), key=leftBkps.count)
    if len(rightBkps) > 0:
        rightBkp = max(set(rightBkps), key=rightBkps.count)

    # TODO: AQUI ES MEJOR QUITARLO DE LA LISTA QUE HACER UNA LSITA NUEVA, PERO DE MOMENTO VA ASI!
    newEvents = []
    # Eliminar los clipping reads que no tengan ese bkp:
    for event in CLIPPING_cluster.events:
        if event.clippedSide == 'left':
            if event.beg == leftBkp:
                # TODO: AQUI IMAGINO QUE SE PUEDE ELIMINAR EL OBJETO DIRECTAMENTE, EN VEZ DE SACARLO DE LA LISTA!
                newEvents.append(event)
        if event.clippedSide == 'right':
            if event.beg == rightBkp:
                # TODO: AQUI ES MEJOR QUITARLO DE LA LISTA QUE HACER UNA LSITA NUEVA, PERO DE MOMENTO VA ASI!
                newEvents.append(event)

    CLIPPING_cluster.events = newEvents

    return leftBkp, rightBkp


def makeConsSeqs(CLIPPING_cluster, clippedSide, seqSide, db, indexDb, outDir):
    '''
    Hacer las cadenas de secuencias para ambos bkps.
    '''
    consensusPath = None
    consensusSeq = None

    clippingEvents = [event for event in CLIPPING_cluster.events if event.clippedSide == clippedSide]

    if len (clippingEvents) > 0:

        consensusPath, consensusSeq = clippingConsensusSeq(clippingEvents, CLIPPING_cluster.id, clippedSide, seqSide, db, indexDb, outDir)
    
    return consensusPath, consensusSeq


def clippingConsensusSeq(clippingEvents, CLIPPING_clusterID, clippedSide, seqSide, db, indexDb, outDir):

    # Retrieve fasta file with sequence from match or clipped side of clipping reads
    supportingReadsFasta = clippingSeq(clippingEvents, CLIPPING_clusterID, clippedSide, seqSide, outDir)

    # Consensus from the previous fasta
    consensusPath, consensusSeq = assembly.getConsensusSeq(supportingReadsFasta, outDir)

    return consensusPath, consensusSeq


def clippingSeq(clippingEvents, CLIPPING_clusterID, clippedSide, seqSide, outDir):
    '''
    Retrieve fasta file with sequence from match or clipped side of clipping reads
    '''

    fastaObj = formats.FASTA()
    fastaDict = {}
    # Determine bkp
    for event in clippingEvents:
        # Si queremos sacar la secuencia del lado de la integracion:
        if seqSide == 'INT':
            if clippedSide == 'left': 
                # cojo la seq desde el principio  hasta el bkp:
                # TODO MIRAR LO DE MAS UUNO
                fastaDict[event.readName] = event.readSeq[:event.readBkp]
            elif clippedSide == 'right': 
                # cojo seq desde el bkp hasta el final:
                fastaDict[event.readName] = event.readSeq[event.readBkp:]

        elif seqSide == 'REF':
            if clippedSide == 'left': 
                # TODO MIRAR LO DE MAS UUNO
                fastaDict[event.readName] = event.readSeq[event.readBkp:]
            elif clippedSide == 'right': 
                fastaDict[event.readName] = event.readSeq[:event.readBkp]
        
        fastaObj.seqDict = fastaDict

    fastaPath = outDir + '/' + str(CLIPPING_clusterID) +'_'+ str(clippedSide) +'_'+ str(seqSide) +'_supportingReads.fa'
    fastaObj.write(fastaPath)

    return fastaPath


def bkpINT(metacluster, consensusPath, db, outDir):

    indexDbSpecificIdentity = databases.buildIdentityDb(metacluster, db, outDir)

    ###### databases ######

    ###### databases ######    

    PAF_file = sequences.getPAFAlign(consensusPath, indexDbSpecificIdentity, outDir)
    PAFObj = formats.PAF()
    PAFObj.read(PAF_file)

    if not os.stat(PAF_file).st_size == 0:
        # DE AQUI SACAMOS LA INFO QUE QUERAMOS DEL VIRUS
        intBkp = [line.tBeg for line in PAFObj.lines][0]
    else:
        intBkp = None

    return intBkp