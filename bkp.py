'''
Module to solve the bkp INGUAL DESPUES SE PUEDE UNIR!
'''
import databases
import formats
import os
import subprocess
import log
import unix


def analyzeClipping(clustersBinDb, confDict, bam, normalBam, mode, db, indexDb, outDir):
    '''
    1. Por cada evento en la BinDB (en este caso especificamente por cada metacluster)
    a. Anadir supporting clipping al metacluster de discordant.
    b. Buscar el bkp (punto mas apoyado por los clippings que acabamos de anadir) y quitar aquellos clipping events que no lo soporten
    c. Hacer las cadenas de secuencias para ambos bkps.
    '''

    dictMetaclusters = {}

    for metacluster in clustersBinDb.collect(['METACLUSTERS']):
        dictMetaclusters[metacluster] = {}

        bkpDir = outDir + '/' + metacluster.ref + '_' + str(metacluster.beg) + '_' + str(metacluster.end)
        unix.mkdir(bkpDir)

        # a. Anadir supporting clipping al metacluster de discordant.
        CLIPPING_cluster = metacluster.supportingCLIPPING(1, confDict, bam, normalBam, mode)

        # b. Buscar el bkp (punto mas apoyado por los clippings que acabamos de anadir) y quitar aquellos clipping events que no lo soporten
        dictMetaclusters[metacluster]['refLeftBkp'], dictMetaclusters[metacluster]['refRightBkp'] = clippingBkp(CLIPPING_cluster)

        print (dictMetaclusters)
        # c. Hacer las cadenas de secuencias para ambos bkps.
        #dictMetaclusters[metacluster]['leftSeq'], dictMetaclusters[metacluster]['rightSeq'] = makeConsSeqs(CLIPPING_cluster, 'REF', db, indexDb, bkpDir)
        leftRefConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'left', 'REF', db, indexDb, bkpDir)[1]
        rightRefConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'right', 'REF', db, indexDb, bkpDir)[1]

        leftIntConsensusPath, leftIntConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'left', 'INT', db, indexDb, bkpDir)
        rightIntConsensusPath, rightIntConsensusSeq = makeConsSeqs(CLIPPING_cluster, 'right', 'INT', db, indexDb, bkpDir)

        leftSeq = leftIntConsensusSeq + '<[INT]' + leftRefConsensusSeq
        rightSeq = rightRefConsensusSeq + '[INT]>' + rightIntConsensusSeq

        dictMetaclusters[metacluster]['leftSeq'], dictMetaclusters[metacluster]['rightSeq'] = leftSeq, rightSeq


        dictMetaclusters[metacluster]['intLeftBkp'] =  bkpINT(metacluster, leftIntConsensusPath, db, bkpDir)
        dictMetaclusters[metacluster]['intRightBkp'] = bkpINT(metacluster, rightIntConsensusPath, db, bkpDir)

        print (dictMetaclusters)
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
            print (event.beg)
            leftBkps.append(event.beg)
        elif event.clippedSide == 'right':
            print (event.clippedSide)
            rightBkps.append(event.beg)

    leftBkp = max(set(leftBkps), key=leftBkps.count)
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


def makeConsSeqs(CLIPPING_cluster, clippedSide, seqSide,  db, indexDb, outDir):
    '''
    Hacer las cadenas de secuencias para ambos bkps.
    '''

    consensusSeq = None

    clippingEvents = [event for event in CLIPPING_cluster.events if event.clippedSide == clippedSide]

    if len (clippingEvents) > 0:

        consensusPath, consensusSeq = clippingConsensusSeq(clippingEvents, CLIPPING_cluster.id, clippedSide, seqSide, db, indexDb, outDir)
    
    return consensusPath, consensusSeq


def clippingConsensusSeq(clippingEvents, CLIPPING_clusterID, clippedSide, seqSide, db, indexDb, outDir):

    # Retrieve fasta file with sequence from match or clipped side of clipping reads
    supportingReadsFasta = clippingSeq(clippingEvents, CLIPPING_clusterID, clippedSide, seqSide, outDir)

    # Consensus from the previous fasta
    consensusPath, consensusSeq = getConsensus(supportingReadsFasta, outDir)

    # De aqui sacas el bkp en lo que esta insertado + ¿secuencia? AQUI PUEDES MIRAR SI FALTAN BASES COMO ANTES!
    # If seq side == INT (pq si es el ref no nos interesa hacer esto)
    print (seqSide)

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
                print (event.readName)
                print (event.readSeq[event.readBkp:])
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

def getConsensus(FASTA_file, outDir):

    ### 2. Make multiple sequence alignment
    # MIRAR SI ESTO ESTA BIEN
    msfPath = FASTA_file.replace("fa", "msf")
    command = 'muscle -in ' + FASTA_file + ' -out ' + msfPath + ' -msf' 
    status = subprocess.call(command, shell=True)

    ### 3. Generate consensus sequence (cons tool from EMBOSS packagge)
    consensusPath = FASTA_file.replace("_supportingReads", "_consensus")

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

    return consensusPath, consensusSeq

def getPAF(FASTA_file, indexDb, outDir):
    # Alineo el fasta consenso
    # TODO ponerlo bien
    PAF_file = FASTA_file.replace("_consensus.fa", "_alignments.paf")

    #err = open(logDir + '/align.err', 'w') 
    command = 'minimap2 ' + indexDb + ' ' + FASTA_file + ' > ' + PAF_file
    status = subprocess.call(command, shell=True)

    if status != 0:
        step = 'ALIGN-INSERT'
        msg = 'Insert alignment failed' 
        log.step(step, msg)
    
    return PAF_file


def bkpINT(metacluster, consensusPath, db, outDir):

    indexDbSpecificIdentity = databases.buildIdentityDb(metacluster, db, outDir)

    ###### databases ######

    ###### databases ######    

    PAF_file = getPAF(consensusPath, indexDbSpecificIdentity, outDir)
    PAFObj = formats.PAF()
    PAFObj.read(PAF_file)

    # DE AQUI SACAMOS LA INFO QUE QUERAMOS DEL VIRUS
    intBkp = [line.tBeg for line in PAFObj.lines][0]
    for line in PAFObj.lines:
        print (line.qName)
        print (line.qLen)
        print (line.qBeg)
        print (line.qEnd)
        print (line.strand)
        print (line.tName)
        print (line.tLen)
        print (line.tBeg)
        print (line.tEnd)
        print (line.nbMatches)

    return intBkp