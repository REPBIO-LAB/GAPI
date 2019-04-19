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


##################################

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

    ###### databases ######
    # Coger la identity del primer evento DISCORDANT que aparezca (pq los clipping no tienen identity)
    identity = next(event.identity for event in metacluster.events if event.type == "DISCORDANT")
    
    specificIdentity = max(set([event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"]), key=[event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"].count)

    specificHeader = '"consensus|' + specificIdentity + '|' + identity + '"'

    #dbSpecificIdentity = outDir + '/' + str(CLIPPING_clusterID) + '_specificIdentity.fa'
    dbSpecificIdentity = outDir + '/' + str(metacluster.id) + '_specificIdentity.fa' 

    # Build database con identity + ref
    # TODO coger la identity que ya hemos reconocido, en lugar de toda la db!
    # 2. Cojo de la db de identities la secuencia que fue asignada como identity:
    # TODO poner bien el status y todo eso
    #err = open(logDir + '/index.err', 'w') 
    command = 'samtools faidx  ' + db + ' ' + specificHeader + ' -o ' + dbSpecificIdentity
    status = subprocess.call(command, shell=True)
    # indexo
    indexDbSpecificIdentity = dbSpecificIdentity.replace("_specificIdentity.fa", "_refIdentityDb.mmi")

    ## DESILENCIAAAAR
    # TODO
    # ponerlo bien
    #err = open(logDir + '/index.err', 'w')
    command = 'minimap2 -k 10 -w 1 -d ' + indexDbSpecificIdentity + ' ' + dbSpecificIdentity 
    status = subprocess.call(command, shell=True)

    if status != 0:
        step = 'BUILD-VIRUS-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)

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