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
import pysam
from collections import Counter

## Internal
import bamtools
import formats
import log
import sequences
import filters
import alignment
import retrotransposons
import callers
from VIGA_LR import bamtools_VIGALR
import assembly
from VIGA_LR import formats_VIGALR
import mappy
import statistics


## [SR CHANGE]
import os
import re

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




# UNUSED!
def check_viral_clipping(events, FASTA_path, name, viralDb, processes, outDir):
    '''
    
    Input:
        1. events: List of clipping events of same clipping side.
        2. FASTA_path: Path to FASTA file with polished sequence.
        3. name: Uniq tag to name files (useful when running function in parallel).
        4. viralDb
        5. processes
        6. outDir: output directory.
    Output:
        1. hits_viral: PAF alignment object
    '''

    # 0. Initialize variables
    hits_viral = None
    name = 'alignments_viral_' + str(name)
    begsEnds = []

    # 1. Aling FASTA against viral db.
    SAM_viral = alignment.alignment_bwa(FASTA_path, viralDb, name, processes, outDir)

    BAM_sorted = bamtools.SAM2BAM(SAM_viral, outDir)

    # Convert SAM to PAF
    PAF_viral = alignment.sam2paf(SAM_viral, name, outDir)

    # 2. Filter match hits
    # NOTE: Vuelve a pasar los filtros pero esta vez de la consenso. Vuelvo a filtrar (pero ahora la consenso) por los mismos parametros que antes
    minParcialMatchVirus = 0
    minTotalMatchVirus = 125 # >= 125
    maxMatchCheckMAPQVirus = 190
    minMAPQVirus = 0
    maxBasePercVirus = 60 # <= 60 (no mucho pero algo filtra)
    minLccVirus = 1.7 # > 1.7

    # Dictionary containing those match hits that passed the filters
    eventsIdentity = bamtools_VIGALR.filterBAM2FastaDict_LR(BAM_sorted, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='LR')
    #print ('eventsIdentity ' + str(eventsIdentity))
    

    # 3. Ver que el alineamiento es CERCA del principio (o del final) del clipping, depende el lado

    # Organize hits according to their corresponding metacluster
    allHits_viral = alignment.organize_hits_paf(PAF_viral)
    # If one sequence of the FASTA has passed filters, store all begginings and ends of all alignments with same sequence (query) name.
    # NOTE: Esto podría ser un problema, ya que estoy evaluando alineamientos que igual no pasaron. De momento lo dejo asi pero tenerlo en cuenta!!
    for query, paf in allHits_viral.items():
        for ali in paf.alignments:
            # miro si paso los filtros:
            if ali.qName in eventsIdentity.keys():
                begsEnds.append(ali.qBeg)
                begsEnds.append(ali.qEnd)
                # Store also the query length
                ali_query_Len = ali.qLen
                #print ('ali_query_Len' + str(ali_query_Len))

    #print ('begsEnds '  + str(begsEnds))
    

    if begsEnds:
        # TODO: poner esto como parametro
        # Si el minimo del principio y final de todos los alignments esta cerca del principio o el maximo esta cerca de la longitud de la query:
        if min(begsEnds) < 250 or abs(ali_query_Len - max(begsEnds)) < 250:
        # Si si que hay return el PAF del alineamiento.
            # El 0 es el que se uso de template asi que es el que va a salir
            # NOTE: Creo que este if no hace falta, porque es solo una secuencia en el fasta, con lo cual si hay alguno sera ese. Asi tampoco haria falta el argumento events.
            if events[0].readName in allHits_viral:
                hits_viral = allHits_viral[events[0].readName]
    
    return hits_viral


def checkEventsIdentity_deprecated(fastaPath, viralDb, processes, outDir):
    ## 3.3 Align consensus inserted sequences into the viral database
    msg = '3.3 Align consensus inserted sequences into the viral database'
    log.info(msg)  
    
    #start_time = time.time()
    SAM_viral = alignment.alignment_bwa(fastaPath, viralDb, 'alignments_viral', processes, outDir)
    #print("--- %s seconds SAM_viral ---" % (time.time() - start_time))

    '''
    # NOTE: I think it's not neccesary
    SAMmapped = outDir + '/alignments_viral_onlyMapped.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'samtools view -b -F4 ' + SAM_viral + ' > ' + SAMmapped
    status = subprocess.call(command, stderr=err, shell=True)
    '''

    # index
    BAM_sorted = bamtools.SAM2BAM(SAM_viral, outDir)

    # Check bam has alignments:
    # TODO: Poner esto de otra forma!!!
    if pysam.idxstats(BAM_sorted) == '*\t0\t0\t0\n':
        emptyPAF = formats.PAF()
        return {}, emptyPAF
    # TODO: put as input options!!!
    minParcialMatchVirus = 0
    minTotalMatchVirus = 125 # >= 125
    maxMatchCheckMAPQVirus = 200
    #maxMatchCheckMAPQVirus = 1000
    minMAPQVirus = 0
    maxBasePercVirus = 60 # <= 60 (no mucho pero algo filtra)
    minLccVirus = 1.7 # > 1.7
    eventsIdentity = virus.filterBAM2FastaDict(BAM_sorted, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='LR')
    
    ## Convert SAM to PAF
    PAF_viral = alignment.sam2paf(SAM_viral, 'alignments_viral', outDir)

    ## Organize hits according to their corresponding metacluster
    allHits_viral = alignment.organize_hits_paf(PAF_viral)

    print ('allHits_viral ' + str(allHits_viral))

   
    return eventsIdentity, allHits_viral

def checkEventsIdentity(fastaPath, viralDb, processes, outDir):
    ## 3.3 Align consensus inserted sequences into the viral database
    msg = '3.3 Align consensus inserted sequences into the viral database'
    log.info(msg)  
    
    #start_time = time.time()
    SAM_viral = alignment.alignment_bwa(fastaPath, viralDb, 'alignments_viral', processes, outDir)
    #print("--- %s seconds SAM_viral ---" % (time.time() - start_time))

    '''
    # NOTE: I think it's not neccesary
    SAMmapped = outDir + '/alignments_viral_onlyMapped.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'samtools view -b -F4 ' + SAM_viral + ' > ' + SAMmapped
    status = subprocess.call(command, stderr=err, shell=True)
    '''

    # index
    BAM_sorted = bamtools.SAM2BAM(SAM_viral, outDir)

    # Check bam has alignments:
    # TODO: Poner esto de otra forma!!!
    if pysam.idxstats(BAM_sorted) == '*\t0\t0\t0\n':
        emptyPAF = formats.PAF()
        return {}, emptyPAF
    # TODO: put as input options!!!
    minParcialMatchVirus = 0
    minTotalMatchVirus = 125 # >= 125
    maxMatchCheckMAPQVirus = 200
    #maxMatchCheckMAPQVirus = 1000
    minMAPQVirus = 0
    maxBasePercVirus = 60 # <= 60 (no mucho pero algo filtra)
    minLccVirus = 1.7 # > 1.7
    eventsIdentity = bamtools_VIGALR.filterBAM2FastaDict_LR(BAM_sorted, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='LR')
    
    ## Convert SAM to PAF
    PAF_viral = alignment.sam2paf(SAM_viral, 'alignments_viral', outDir)

    ## Organize hits according to their corresponding metacluster
    allHits_viral = alignment.organize_hits_paf(PAF_viral)

    print ('allHits_viral ' + str(allHits_viral))

    # Minimap


   
    return eventsIdentity, allHits_viral

def checkBkpProximity_wSide(event, nombre, allHits_viral, buffer):
    '''
    '''
    bkpProximity = None
    begsEnds = []

    for ali in allHits_viral[nombre].alignments:

        begsEnds.append(ali.qBeg)
        begsEnds.append(ali.qEnd)
        # Store also the query length
        ali_query_Len = ali.qLen

        print ('ali_query_Len ' + str(ali_query_Len))
        print ('begsEnds ' + str(begsEnds))

        if begsEnds:
            # TODO: poner esto como parametro
            # Si el minimo del principio y final de todos los alignments esta cerca del principio o el maximo esta cerca de la longitud de la query:
            if event.clippedSide == 'right':
                if min(begsEnds) <= buffer:
                # Si si que hay return el PAF del alineamiento.
                    # El 0 es el que se uso de template asi que es el que va a salir
                    # NOTE: Creo que este if no hace falta, porque es solo una secuencia en el fasta, con lo cual si hay alguno sera ese. Asi tampoco haria falta el argumento events.
                    bkpProximity = min(begsEnds)
            elif event.clippedSide == 'left':
                if abs(ali_query_Len - max(begsEnds)) <= buffer:
                    bkpProximity = max(begsEnds)
    
    return bkpProximity


def checkBkpProximity_woSide(event, nombre, allHits_viral, buffer):
    '''
    '''
    bkpProximity = None
    begsEnds = []


    for ali in allHits_viral[nombre].alignments:

        begsEnds.append(ali.qBeg)
        begsEnds.append(ali.qEnd)
        # Store also the query length
        ali_query_Len = ali.qLen

        print ('ali_query_Len ' + str(ali_query_Len))
        print ('begsEnds ' + str(begsEnds))


        if begsEnds:
            # TODO: poner esto como parametro
            # Si el minimo del principio y final de todos los alignments esta cerca del principio o el maximo esta cerca de la longitud de la query:
            if min(begsEnds) < buffer:
            # Si si que hay return el PAF del alineamiento.
                # El 0 es el que se uso de template asi que es el que va a salir
                # NOTE: Creo que este if no hace falta, porque es solo una secuencia en el fasta, con lo cual si hay alguno sera ese. Asi tampoco haria falta el argumento events.
                bkpProximity = min(begsEnds)
            elif abs(ali_query_Len - max(begsEnds)) < buffer:
                bkpProximity = max(begsEnds)

    return bkpProximity


###### KIND OF GENERAL FUNCTIONS, THAT BELONGS TO VIRUS FOR THE MOMENT ######

#### Module alignments:

def alignment_minimap2_SAM(FASTA, index, fileName, processes, outDir):
    '''
    Align a set of sequence into a reference with minimap2

    Input:
        1. FASTA: Path to FASTA file with sequences to align
        2. index: Path to the the index of the reference in .mmi format (generated with minimap2)
        3. fileName: output file will be named accordingly
        4. processes: Number of processes used by minimap2
        5. outDir: Output directory

    Output:
        1. PAF: Path to PAF file containing input sequences alignments or 'None' if alignment failed 
    '''
    ## Align the sequences into the reference
    # Note, condider to use -Y to get soft clippings for supplementary alignments
    SAM = outDir + '/' + fileName + '.sam'
    err = open(outDir + '/align.err', 'w') 
    command = 'minimap2 -a -t ' + str(processes) + ' ' + index + ' ' + FASTA + ' > ' + SAM
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'ALIGN'
        msg = 'Local alignment failed' 
        log.step(step, msg)

    return SAM

#### Module clusters:

# UNUSED
def check_metaclusters_identity(SV_type, metacluster, eventsIdentity):
    
    newAllMetaclusters = {}
    newAllMetaclusters[SV_type] = []
    identity = []
    specificIdentity = []


    #for metacluster in metaclusters:
    for event in metacluster.events:
        if event.readName in eventsIdentity.keys(): 
            print ('eventIDEnTITY ' + str(event.readName) + ' ' + str(eventsIdentity[event.readName]))
            for ident in eventsIdentity[event.readName]:
                splitedtName = ident.split('|')
                identity.append(splitedtName[0])
                if len(splitedtName) > 1:
                    specificIdentity.append(splitedtName[1])
            event.identity = identity
            event.specificIdentity = specificIdentity
            #metacluster.identity = eventsIdentity[event.readName]

    newAllMetaclusters[SV_type].append(metacluster)
    
    return newAllMetaclusters

def check_metaclusters_identity_allHits(SV_type, metacluster, eventsIdentity, allHits_viral, outDir):
    '''
    fastaDict es un dictionario FASTA que puede contener lo siguiente:
        - Las secuencias consenso de los lados de las INS_noOverlap
        - La secuencia consenso de las INS que hacen overlap.
        - La secuencia consenso de todos aquellos metaclusters cuya proporcion de eventos con la misma identidad no es suficiente para determinar la identidad del metacluster directamente.
    '''
    
    newAllMetaclusters = {}
    newAllMetaclusters[SV_type] = []
    identity = []
    specificIdentity = []
    #begsEnds = []
    fastaDict = {}

    # Special for INS_noOverlap and INS overlap, because they already have a consensus.
    if SV_type == 'INS_noOverlap' and (metacluster.consRightSeq or metacluster.consLeftSeq):
        INS_noOverlapDict = {}
        if metacluster.consRightSeq:
            for event in metacluster.events:
                if event.clippedSide == 'right':
                    INS_noOverlapDict[event.readName + '_'+ str(metacluster.beg)] = ''.join(metacluster.consRightSeq)
                    break
        if metacluster.consLeftSeq:
            for event in metacluster.events:
                if event.clippedSide == 'left':
                    INS_noOverlapDict[event.readName + '_'+ str(metacluster.beg)] = ''.join(metacluster.consLeftSeq)
                    break
        fastaDict.update(INS_noOverlapDict)

        newAllMetaclusters[SV_type].append(metacluster)
        return newAllMetaclusters, fastaDict

    elif SV_type == 'INS' and not any([event.type == 'INS' for event in metacluster.events]):
        INS_OverlapDict = {}
        if metacluster.consensusEvent:
            # TODO: Fix this in another way:
            if metacluster.consensusEvent.readName != 'CONTIG':
                INS_OverlapDict[metacluster.consensusEvent.readName + '_'+ str(metacluster.beg)] = metacluster.consensusEvent.readSeq
            else:
                INS_OverlapDict[metacluster.events[0].readName + '_'+ str(metacluster.beg)] = metacluster.consensusEvent.readSeq


        fastaDict.update(INS_OverlapDict)

        newAllMetaclusters[SV_type].append(metacluster)
        return newAllMetaclusters, fastaDict


    # Check bkp proximity in CLIPPINGS
    for event in metacluster.events:
        begsEnds = []
        if event.type == 'CLIPPING':
            if event.readName in allHits_viral.keys() and event.readName in eventsIdentity.keys(): 

                # If one sequence of the FASTA has passed filters, store all begginings and ends of all alignments with same sequence (query) name.
                # NOTE: Esto podría ser un problema, ya que estoy evaluando alineamientos que igual no pasaron. De momento lo dejo asi pero tenerlo en cuenta!!

                for ali in allHits_viral[event.readName].alignments:

                    begsEnds.append(ali.qBeg)
                    begsEnds.append(ali.qEnd)
                    # Store also the query length
                    ali_query_Len = ali.qLen
                    print ('ali_query_Len' + str(ali_query_Len))

                    print ('begsEnds '  + str(begsEnds))
                    

                    if begsEnds:
                        # TODO: poner esto como parametro
                        # Si el minimo del principio y final de todos los alignments esta cerca del principio o el maximo esta cerca de la longitud de la query:
                        if min(begsEnds) < 250:
                        # Si si que hay return el PAF del alineamiento.
                            # El 0 es el que se uso de template asi que es el que va a salir
                            # NOTE: Creo que este if no hace falta, porque es solo una secuencia en el fasta, con lo cual si hay alguno sera ese. Asi tampoco haria falta el argumento events.
                            event.bkpProximity = min(begsEnds)
                        elif abs(ali_query_Len - max(begsEnds)) < 250:
                            event.bkpProximity = max(begsEnds)
                    # NOTE: event.bkpProximity es el principio o el final del alig que esta cerca del bkp


    # Check events identity

    for event in metacluster.events:
        print ('event.type != CLIP1 ' + str(metacluster.beg))
        if event.type != 'CLIPPING':
            print ('event.type != CLIP2 ' + str(metacluster.beg))
            if event.readName in allHits_viral.keys() and event.readName in eventsIdentity.keys():
                print ('event.type != CLIP3 ' + str(metacluster.beg))
                for alig in allHits_viral[event.readName].alignments:
                    splitedtName = alig.tName.split('|')
                    identity.append(splitedtName[0])
                    if len(splitedtName) > 1:
                        specificIdentity.append(splitedtName[1])
                print ('event.type != CLIP4 ' + str(metacluster.beg))
                event.identity = list(set(identity))
                event.specificIdentity = list(set(specificIdentity))
                #metacluster.identity = eventsIdentity[event.readName]
        
        # NOTE: asi los clipping tienen solo una identity que es la del alineamiento que esta mas cerca del bkp.
        elif event.type == 'CLIPPING':
            if (event.readName in allHits_viral.keys()):
                print ('EO1 ' + str(event.readName))
            print ('EO2 ' + str(event.bkpProximity) + ' ' + str(event.readName))
            if event.readName in eventsIdentity.keys():
                print ('EO3 ' + str(event.readName))
            if (event.readName in allHits_viral.keys()) and (event.bkpProximity != None) and (event.readName in eventsIdentity.keys()): 
                for alig in allHits_viral[event.readName].alignments:
                    if alig.qBeg == event.bkpProximity or alig.qEnd == event.bkpProximity:
                        splitedtName = alig.tName.split('|')
                        identity.append(splitedtName[0])
                        if len(splitedtName) > 1:
                            specificIdentity.append(splitedtName[1])
                event.identity = list(set(identity))
                event.specificIdentity = list(set(specificIdentity))
                #metacluster.identity = eventsIdentity[event.readName]

    # Check metacluster identity

    # Segun el procentaje de eventos con la misma ident le doy ident o no:
    identities = list(itertools.chain(*[event.identity for event in metacluster.events if event.identity != None]))
    print ('EGAAAAH')
    print (metacluster.beg)
    if identities:
        print (identities)
        countIdentities = Counter(identities)
        print (countIdentities)
        mostCommonIdentity = countIdentities.most_common(1)[0]
        print (mostCommonIdentity)
        identitiesProportion = mostCommonIdentity[1]/metacluster.nbEvents()[0]
        print (identitiesProportion)
        # La identitdad mas comun tiene que estar apoyada al mens por el 70% de los eventos
        if identitiesProportion >= 0.70:
            metacluster.SV_features['IDENTITY'] = countIdentities.most_common(3)
            specificIdentities = list(itertools.chain(*[event.specificIdentity for event in metacluster.events if event.identity != None]))
            if specificIdentities:
                countSpecificIdentities = Counter(specificIdentities)
                metacluster.SV_features['SPECIDENTITY'] = countSpecificIdentities.most_common(3)
                        
        # If identitiesProportion == 0.5 (for example, if tehre are 2 reads and only one supporting), do consensus:
        elif identitiesProportion < 0.70 and identitiesProportion >= 0.18:
            name = 'lowIdentitiesProportion_' + str(metacluster.beg)
            polishDir = outDir + '/lowIdentitiesProportion/' + name
            unix.mkdir(polishDir)
            if SV_type == 'INS':
                polished_FASTA = polish_INS_sequence(metacluster, metacluster.events, name, polishDir)
            else:
                polished_FASTA = polish_CLIP_sequence(metacluster, metacluster.events, name, polishDir)
            if polished_FASTA:
                # Read FASTA
                fastaObj = formats.FASTA()
                fastaObj.read(polished_FASTA)
                fastaDict.update(fastaObj.seqDict)
                '''
                # If there is consensus, put identity:
                if fastaObj.seqDict.values():
                    metacluster.SV_features['IDENTITY'] = countIdentities.most_common(3)
                    specificIdentities = list(itertools.chain(*[event.specificIdentity for event in metacluster.events if event.identity != None]))
                    if specificIdentities:
                        countSpecificIdentities = Counter(specificIdentities)
                        metacluster.SV_features['SPECIDENTITY'] = countSpecificIdentities.most_common(3)
                # Alineo
                '''

            unix.rm([polishDir])


    newAllMetaclusters[SV_type].append(metacluster)
    print ('fastaDict0' + str(fastaDict))
    return newAllMetaclusters, fastaDict

def supportingEventsFasta(metacluster, outDir):
    fastaEventsObj = formats_VIGALR.FASTA()
    for event in metacluster.events:
            if event.type == 'INS':
                fastaEventsObj.seqDict[event.readName]=event.pick_insert()
            elif event.type == 'CLIPPING':
                fastaEventsObj.seqDict[event.readName]=event.clipped_seq()

    # Set output FASTA file name
    soloBND_eventsFasta = outDir + "/events.fasta"
    # Write output FASTA
    fastaEventsObj.write(soloBND_eventsFasta, 'append', True)


def get_consensus_nonClassisSVTypes(metacluster, SV_type, processes, outDir):
    '''
    Make consensus of clipping sequences of metaclusters and determine their identity.
    '''
    # 0. Create output directory
    metaOutDir = outDir +'/'+ str(metacluster.id) + '_' + str(metacluster.bkpPos)
    unix.mkdir(metaOutDir)
    
    # A. Analize BND (break-ends):
    if SV_type == 'BND':
        # A1. Get polished clipping sequence
        polished_clippling_FASTA = polish_CLIP_sequence(metacluster, metacluster.events, metacluster.bkpPos, metaOutDir)
        if polished_clippling_FASTA:
            # Read FASTA
            fastaObj = formats.FASTA()
            fastaObj.read(polished_clippling_FASTA)
        
            if metacluster.events[0].clippedSide == 'right':
                metacluster.consRightSeq = list(fastaObj.seqDict.values())
        
            elif metacluster.events[0].clippedSide == 'left':
                metacluster.consLeftSeq = list(fastaObj.seqDict.values())
    
    # TODO: sustituir
    #elif SV_type == 'INS_noOverlap':
    # B. Analize Big_INS (big insertions -those whose ends do not overlap-):
    else:

        # B0. Initialize variables

        eventsClipRight=[]
        eventsClipLeft=[]
        #polished_clippling_FASTAs = []

        # B1. Split metacluster events based on clippedSide
        for event in metacluster.events:
            if event.clippedSide == 'right':
                eventsClipRight.append(event)
            elif event.clippedSide == 'left':
                eventsClipLeft.append(event)

        # B2. Analize right clipping events and left clipping events separately.
        for eventList in [eventsClipRight, eventsClipLeft]:
            if eventList:
                # B3. Get PAF alignment object with viral alignments
                polished_clippling_FASTA = polish_CLIP_sequence(metacluster, eventList, metacluster.bkpPos, metaOutDir)

                if polished_clippling_FASTA and all([event.clippedSide == 'right' for event in eventList]):
                    fastaObjR = formats.FASTA()
                    fastaObjR.read(polished_clippling_FASTA)
                    # Consensus fasta Right:
                    metacluster.consRightSeq = list(fastaObjR.seqDict.values())

                if polished_clippling_FASTA and all([event.clippedSide == 'left' for event in eventList]):
                    fastaObjL = formats.FASTA()
                    fastaObjL.read(polished_clippling_FASTA)
                    # Consensus fasta Left:
                    metacluster.consLeftSeq = list(fastaObjL.seqDict.values())
    print ('metacluster.consRightSeq ' + str(metacluster.beg) +' '+ str(metacluster.consRightSeq))
    print ('metacluster.consLeftSeq ' + str(metacluster.beg) +' '+ str(metacluster.consLeftSeq))
    return metacluster

# UNUSED

def determine_soloBND_type(metacluster, allHits_genome, groupedEntries, allHits_viral):

        ## 4.1 Collect consensus inserted sequence hits
        metaId = str(metacluster.ref) + ':' + str(metacluster.beg) + '-' + str(metacluster.end)

        ## Hits in the reference genome
        if metaId in allHits_genome:
            hits_genome = allHits_genome[metaId]
        
        else:
            hits_genome = formats.PAF()

        ## Hits in the reference genome (splice-aware alignment)
        if metaId in groupedEntries:
            hits_splicing = groupedEntries[metaId]

        else:
            hits_splicing = []

        ## Hits in the viral database
        if metaId in allHits_viral:
            hits_viral = allHits_viral[metaId]

        else:
            hits_viral = formats.PAF()

        ## 4.2 Insertion type inference
        #metacluster.determine_INS_type(hits_genome, hits_splicing, hits_viral, annotations['REPEATS'], annotations['TRANSDUCTIONS'], annotations['EXONS'])
        metacluster.determine_INS_type(hits_genome, hits_splicing, hits_viral, None, None, None, types2Search=['VIRUS'])

        return metacluster

def checkLowProportion(allMetaclustersIdentity, event, metacluster, allHits_viral, eventsIdentity, bufferBND, bufferNoOverlap, sideBND, sideNoOverlap):
    '''
    '''
    bkpProximity = None
    SV_TypeTested = None
    identity = []
    specificIdentity = []

    # TODO: Do this in another way to find the event that is the key
    # Find the event that was used as consensus name

    nombre =  event.readName + '_'+ str(metacluster.beg)
    
    for alig in allHits_viral[nombre].alignments:
        splitedtName = alig.tName.split('|')
        identity.append(splitedtName[0])
        if len(splitedtName) > 1:
            specificIdentity.append(splitedtName[1])
    
    # Check if it's near bkp
    if 'BND' in allMetaclustersIdentity.keys():
        
        if metacluster in allMetaclustersIdentity['BND']:
            SV_TypeTested = True

            if sideBND:
                bkpProximity = checkBkpProximity_wSide(event, nombre, allHits_viral, bufferBND)
            else:
                bkpProximity = checkBkpProximity_woSide(event, nombre, allHits_viral, bufferBND)

    # Check if it's near bkp
    if 'INS_noOverlap' in allMetaclustersIdentity.keys():
        if metacluster in allMetaclustersIdentity['INS_noOverlap']:
            SV_TypeTested = True

            if sideNoOverlap:
                bkpProximity = checkBkpProximity_wSide(event, nombre, allHits_viral, bufferNoOverlap)
            else:
                bkpProximity = checkBkpProximity_woSide(event, nombre, allHits_viral, bufferNoOverlap)

            event.bkpProximity = bkpProximity
                    
    return bkpProximity, SV_TypeTested, identity, specificIdentity


#### Module output:

#UNUSED (but maybe it will be used)
def checkSpecIdent_LR(metacluster):
    checkSpecIdent = {}
    for event in metacluster.events:
        if event.specificIdentity != None:

            if any(event.specificIdentity) in checkSpecIdent.keys():
                checkSpecIdent[event.specificIdentity] += 1
            else:
                checkSpecIdent[event.specificIdentity] = 1
    
    if not checkSpecIdent:
        SPECIDENT = None
    else: 
        SPECIDENT = ",".join("{}:{}".format(k, checkSpecIdent[k]) for k in sorted(checkSpecIdent))

    return SPECIDENT

def INS2VCF_junction(metaclusters, index, refLengths, source, build, species, outName, outDir):
    '''
    Write INS calls into a VCF file

    Input:
        1. metaclusters: list containing list of INS metaclusters
        2. index: minimap2 index for the reference genome 
        3. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        4. source: software version used to generate the insertion calls
        5. build: reference genome build
        6. species: specie
        7. outName: Output file name
        8. outDir: Output directory

    Output: vcf file containing identified metaclusters
    '''
    ## 1. Initialize VCF 
    VCF = formats.VCF()

    ## 2. Create header
    ## Define info 
    info = {'VTYPE': ['.', 'String', 'Type of variant'], \
            'ITYPE': ['.', 'String', 'Type of structural variant'], \
            'MECHANISM': ['.', 'String', 'Insertion mechanism'], \
            'FAM': ['.', 'String', 'Repeat family'], \
            'SUBFAM': ['.', 'String', 'Repeat subfamily'], \
            'GERMDB': ['.', 'String', 'List of germline variation databases where the variant is reported'], \
            'CIPOS': ['2', 'Integer', 'Confidence interval around POS for imprecise variants'], \
            'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
            'NBEXONS': ['1', 'Integer', 'Number of exons for a processed pseudogene insertion'], \
            'SRCGENE': ['.', 'String', 'Source gene for a processed psendogene insertion'], \
            'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
            'REGION': ['.', 'String', 'Genomic region where insertion occurs'], \
            'GENE': ['.', 'String', 'HUGO gene symbol'], \
            'REP': ['.', 'String', 'Families for annotated repeats at the insertion region'], \
            'REPSUB': ['.', 'String', 'Subfamilies for annotated repeats at the insertion region'], \
            'DIST': ['.', 'Integer', 'Distance between insertion breakpoint and annotated repeat'], \
            'NBTOTAL': ['1', 'Integer', 'Total number of insertion supporting reads'], \
            'NBTUMOR': ['1', 'Integer', 'Number of insertion supporting reads in the tumour'], \
            'NBNORMAL': ['1', 'Integer', 'Number of insertion supporting reads in the normal'], \
            'NBSPAN': ['1', 'Integer', 'Number of spanning supporting reads'], \
            'NBCLIP': ['1', 'Integer', 'Number of clipping supporting reads'], \
            'LEN': ['1', 'Integer', 'Insertion length'], \
            'CV': ['1', 'Float', 'Length coefficient of variation'], \
            'RTLEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
            'TRUN5LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
            'TRUN3LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
            'FULL': ['0', 'Flag', 'Full length mobile element'], \
            'TDLEN': ['1', 'Integer', 'Transduction length'], \
            'INVLEN': ['1', 'Integer', '5-inversion length'], \
            'PERCR': ['1', 'Float', 'Percentage of inserted sequence that has been resolved'], \
            'QHITS': ['.', 'String', 'Coordinates for inserted sequence hits on the reference'], \
            'THITS': ['.', 'String', 'Inserted sequence hits on the reference'], \
            'RTCOORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
            'POLYA': ['0', 'Flag', 'PolyA tail identified'], \
            'INSEQ': ['.', 'String', 'Inserted sequence'], \
            'IDENTITY': ['.', 'String', 'Inserted identity'], \
            'SPECIDENTITY': ['.', 'String', 'Inserted identity'], \
            'INS': ['.', 'String', 'Insertion supporting reads'], \
            'CLIP': ['.', 'String', 'Clipping supporting reads'], \
            'BKPCSEQR': ['.', 'String', 'Consensus sequence at Right BKP.'], \
            'BKPCSEQL': ['.', 'String', 'Consensus sequence at Left BKP.'], \
            'MAPQ': ['1', 'Integer', 'Average of MAPQ of discordant reads supporting the insertion'], \
            'SUPPOBKP': ['.', 'String', 'Breakpoints that support the integrations (RIGHT or/and LEFT)'], \
            }
            
    ## Create header
    VCF.create_header(source, build, species, refLengths, info, {}, [])

    ## 3. Add insertion calls to the VCF
    ## 3.1 Load reference index
    reference = mappy.Aligner(fn_idx_in=index) # comment as time consuming

    ## 3.2 Iterate over INS metaclusters
    for metacluster in metaclusters:

        print ('metacluster.bridge ' + str(metacluster.bridge))

        ## Collect insertion basic features
        CHROM = metacluster.ref
        POS, CIPOS = metacluster.mean_pos()
        ID = '.'
        REF = reference.seq(CHROM, POS, POS + 1)
        ALT = '<INS>'
        QUAL = '.'
        FILTER = 'PASS' if not metacluster.failedFilters else ','.join(metacluster.failedFilters)
        
        ## Collect extra insertion features to include at info field
        INFO = {}
        repeats = metacluster.repeatAnnot if hasattr(metacluster, 'repeatAnnot') else []   
        INS_readNames = [ INS.readName for INS in metacluster.events if INS.type == 'INS' ]
        CLIP_readNames = [ CLIP.readName for CLIP in metacluster.events if CLIP.type == 'CLIPPING' ]     
        eventsMAPQ = [ event.mapQual for event in metacluster.events if event.mapQual != None ]
        if eventsMAPQ:
            eventsMAPQMean = int(statistics.mean(eventsMAPQ)) if eventsMAPQ else None
        else:
            eventsMAPQMean = None
        
        # Supporting bkp:
        suppoBkps = []
        for event in metacluster.events:
            if event.bkpProximity:
                suppoBkps.append(event.clippedSide)
        suppoBkp = list(set(suppoBkps))

        INFO['VTYPE'] = metacluster.mutOrigin
        INFO['ITYPE'] = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        INFO['MECHANISM'] = metacluster.SV_features['MECHANISM'] if 'MECHANISM' in metacluster.SV_features else None        
        INFO['FAM'] = ','.join(metacluster.SV_features['FAMILY']) if ('FAMILY' in metacluster.SV_features and metacluster.SV_features['FAMILY']) else None
        INFO['SUBFAM'] = ','.join(metacluster.SV_features['SUBFAMILY']) if ('SUBFAMILY' in metacluster.SV_features and metacluster.SV_features['SUBFAMILY']) else None
        INFO['GERMDB'] = metacluster.germlineDb if hasattr(metacluster, 'germlineDb') else None       
        INFO['CIPOS'] = str(CIPOS[0]) + ',' + str(CIPOS[1]) 
        INFO['CYTOID'] = ','.join(metacluster.SV_features['CYTOBAND']) if ('CYTOBAND' in metacluster.SV_features and metacluster.SV_features['CYTOBAND']) else None
        INFO['NBEXONS'] = metacluster.SV_features['NB_EXONS'] if 'NB_EXONS' in metacluster.SV_features else None
        INFO['SRCGENE'] = ','.join(metacluster.SV_features['SOURCE_GENE']) if 'SOURCE_GENE' in metacluster.SV_features else None
        INFO['STRAND'] = metacluster.SV_features['STRAND'] if 'STRAND' in metacluster.SV_features else None
        INFO['REGION'], INFO['GENE'] = metacluster.geneAnnot if hasattr(metacluster, 'geneAnnot') else (None, None)
        INFO['REP'] = ','.join([repeat['family'] for repeat in repeats]) if repeats else None 
        INFO['REPSUB'] = ','.join([repeat['subfamily'] for repeat in repeats]) if repeats else None   
        INFO['DIST'] = ','.join([str(repeat['distance']) for repeat in repeats]) if repeats else None
        INFO['NBTOTAL'], INFO['NBTUMOR'], INFO['NBNORMAL'] = str(metacluster.nbTotal), str(metacluster.nbTumour), str(metacluster.nbNormal) 
        INFO['NBSPAN'], INFO['NBCLIP'] = str(metacluster.nbINS), str(metacluster.nbCLIPPING)
        #INFO['LEN'] = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        INFO['CV'] = metacluster.cv
        INFO['RTLEN'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['TRUN5LEN'] = metacluster.SV_features['TRUNCATION_5_LEN'] if 'TRUNCATION_5_LEN' in metacluster.SV_features else None
        INFO['TRUN3LEN'] = metacluster.SV_features['TRUNCATION_3_LEN'] if 'TRUNCATION_3_LEN' in metacluster.SV_features else None
        INFO['FULL'] = metacluster.SV_features['IS_FULL'] if 'IS_FULL' in metacluster.SV_features else None
        INFO['TDLEN'] = metacluster.SV_features['TRANSDUCTION_LEN'] if 'TRANSDUCTION_LEN' in metacluster.SV_features else None
        INFO['INVLEN'] = metacluster.SV_features['INVERSION_LEN'] if 'INVERSION_LEN' in metacluster.SV_features else None
        INFO['PERCR'] = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        #INFO['QHITS'] = None if metacluster.insertHits is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertHits.alignments ])
        #INFO['THITS'] = None if metacluster.insertHits is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['RTCOORD'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['POLYA'] = metacluster.SV_features['POLYA'] if 'POLYA' in metacluster.SV_features else None
        INFO['IDENTITY'] = metacluster.SV_features['IDENTITY'] if 'IDENTITY' in metacluster.SV_features else None
        INFO['SPECIDENTITY'] = metacluster.SV_features['SPECIDENTITY'] if 'SPECIDENTITY' in metacluster.SV_features else None
        #INFO['SPECIDENTITY'] = checkSpecIdent_LR(metacluster)
        INFO['SPECIDENTITY'] = None
        INFO['INS'] = ','.join(INS_readNames) if INS_readNames else None
        INFO['CLIP'] = ','.join(CLIP_readNames) if CLIP_readNames else None
        INFO['BKPCSEQR'] = metacluster.consLeftSeq if metacluster.consLeftSeq != None else None
        INFO['BKPCSEQL'] = metacluster.consLeftSeq if metacluster.consLeftSeq != None else None
        INFO['MAPQ'] = eventsMAPQMean
        #INFO['INSEQ'] = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None
        INFO['SUPPOBKP'] = suppoBkp

        ## Create VCF variant object
        fields = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, {}]

        ## Add variant to the VCF
        INS = formats.VCF_variant(fields)
        VCF.add(INS)
        
    ## 4. Sort VCF
    VCF.sort()

    ## 5. Write VCF in disk
    IDS = ['VTYPE', 'ITYPE', 'MECHANISM', 'FAM', 'SUBFAM', 'GERMDB', 'CIPOS', 'CYTOID', \
           'NBEXONS', 'SRCGENE', 'STRAND', 'REGION', 'GENE', 'REP', 'REPSUB', 'DIST', \
           'NBTOTAL', 'NBTUMOR', 'NBNORMAL', 'NBSPAN', 'NBCLIP', 'LEN', 'CV', 'RTLEN', \
           'TRUN5LEN', 'TRUN3LEN', 'FULL', 'TDLEN', 'INVLEN', 'PERCR', \
           'QHITS', 'THITS', 'RTCOORD', 'POLYA', \
           'IDENTITY', 'SPECIDENTITY', 'MAPQ', 'SUPPOBKP', 'INSEQ', 'INS', 'CLIP', 'BKPCSEQR', 'BKPCSEQL']

    VCF.write(IDS, [], outName, outDir)


def INS2VCF(metaclusters, index, refLengths, source, build, species, outName, outDir):
    '''
    Write INS calls into a VCF file

    Input:
        1. metaclusters: list containing list of INS metaclusters
        2. index: minimap2 index for the reference genome 
        3. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        4. source: software version used to generate the insertion calls
        5. build: reference genome build
        6. species: specie
        7. outName: Output file name
        8. outDir: Output directory

    Output: vcf file containing identified metaclusters
    '''
    ## 1. Initialize VCF 
    VCF = formats.VCF()

    ## 2. Create header
    ## Define info 
    info = {'VTYPE': ['.', 'String', 'Type of variant'], \
            'ITYPE': ['.', 'String', 'Type of structural variant'], \
            'MECHANISM': ['.', 'String', 'Insertion mechanism'], \
            'FAM': ['.', 'String', 'Repeat family'], \
            'SUBFAM': ['.', 'String', 'Repeat subfamily'], \
            'GERMDB': ['.', 'String', 'List of germline variation databases where the variant is reported'], \
            'CIPOS': ['2', 'Integer', 'Confidence interval around POS for imprecise variants'], \
            'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
            'NBEXONS': ['1', 'Integer', 'Number of exons for a processed pseudogene insertion'], \
            'SRCGENE': ['.', 'String', 'Source gene for a processed psendogene insertion'], \
            'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
            'REGION': ['.', 'String', 'Genomic region where insertion occurs'], \
            'GENE': ['.', 'String', 'HUGO gene symbol'], \
            'REP': ['.', 'String', 'Families for annotated repeats at the insertion region'], \
            'REPSUB': ['.', 'String', 'Subfamilies for annotated repeats at the insertion region'], \
            'DIST': ['.', 'Integer', 'Distance between insertion breakpoint and annotated repeat'], \
            'NBTOTAL': ['1', 'Integer', 'Total number of insertion supporting reads'], \
            'NBTUMOR': ['1', 'Integer', 'Number of insertion supporting reads in the tumour'], \
            'NBNORMAL': ['1', 'Integer', 'Number of insertion supporting reads in the normal'], \
            'NBSPAN': ['1', 'Integer', 'Number of spanning supporting reads'], \
            'NBCLIP': ['1', 'Integer', 'Number of clipping supporting reads'], \
            'LEN': ['1', 'Integer', 'Insertion length'], \
            'CV': ['1', 'Float', 'Length coefficient of variation'], \
            'RTLEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
            'TRUN5LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
            'TRUN3LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
            'FULL': ['0', 'Flag', 'Full length mobile element'], \
            'TDLEN': ['1', 'Integer', 'Transduction length'], \
            'INVLEN': ['1', 'Integer', '5-inversion length'], \
            'PERCR': ['1', 'Float', 'Percentage of inserted sequence that has been resolved'], \
            'QHITS': ['.', 'String', 'Coordinates for inserted sequence hits on the reference'], \
            'THITS': ['.', 'String', 'Inserted sequence hits on the reference'], \
            'RTCOORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
            'POLYA': ['0', 'Flag', 'PolyA tail identified'], \
            'INSEQ': ['.', 'String', 'Inserted sequence'], \
            'IDENTITY': ['.', 'String', 'Inserted identity'], \
            'SPECIDENTITY': ['.', 'String', 'Inserted identity'], \
            'INS': ['.', 'String', 'Insertion supporting reads'], \
            'CLIP': ['.', 'String', 'Clipping supporting reads'], \
            'BKPCSEQR': ['.', 'String', 'Consensus sequence at Right BKP.'], \
            'BKPCSEQL': ['.', 'String', 'Consensus sequence at Left BKP.'], \
            'MAPQ': ['1', 'Integer', 'Average of MAPQ of discordant reads supporting the insertion'], \
            }
            
    ## Create header
    VCF.create_header(source, build, species, refLengths, info, {}, None)

    ## 3. Add insertion calls to the VCF
    ## 3.1 Load reference index
    reference = mappy.Aligner(fn_idx_in=index) # comment as time consuming

    ## 3.2 Iterate over INS metaclusters
    for metacluster in metaclusters:

        ## Collect insertion basic features
        CHROM = metacluster.ref
        POS, CIPOS = metacluster.mean_pos()
        ID = '.'
        REF = reference.seq(CHROM, POS, POS + 1)
        ALT = '<INS>'
        QUAL = '.'
        FILTER = 'PASS' if not metacluster.failedFilters else ','.join(metacluster.failedFilters)
        
        ## Collect extra insertion features to include at info field
        INFO = {}
        repeats = metacluster.repeatAnnot if hasattr(metacluster, 'repeatAnnot') else []        
        INS_readNames = [ INS.readName for INS in metacluster.events if INS.type == 'INS' ]
        CLIP_readNames = [ CLIP.readName for CLIP in metacluster.events if CLIP.type == 'CLIPPING' ]
        eventsMAPQ = [ event.mapQual for event in metacluster.events if event.mapQual != None ]
        if eventsMAPQ:
            eventsMAPQMean = int(statistics.mean(eventsMAPQ)) if eventsMAPQ else None
        else:
            eventsMAPQMean = None
        
        INFO['VTYPE'] = metacluster.mutOrigin
        INFO['ITYPE'] = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        INFO['MECHANISM'] = metacluster.SV_features['MECHANISM'] if 'MECHANISM' in metacluster.SV_features else None        
        INFO['FAM'] = ','.join(metacluster.SV_features['FAMILY']) if ('FAMILY' in metacluster.SV_features and metacluster.SV_features['FAMILY']) else None
        INFO['SUBFAM'] = ','.join(metacluster.SV_features['SUBFAMILY']) if ('SUBFAMILY' in metacluster.SV_features and metacluster.SV_features['SUBFAMILY']) else None
        INFO['GERMDB'] = metacluster.germlineDb if hasattr(metacluster, 'germlineDb') else None       
        INFO['CIPOS'] = str(CIPOS[0]) + ',' + str(CIPOS[1]) 
        INFO['CYTOID'] = ','.join(metacluster.SV_features['CYTOBAND']) if ('CYTOBAND' in metacluster.SV_features and metacluster.SV_features['CYTOBAND']) else None
        INFO['NBEXONS'] = metacluster.SV_features['NB_EXONS'] if 'NB_EXONS' in metacluster.SV_features else None
        INFO['SRCGENE'] = ','.join(metacluster.SV_features['SOURCE_GENE']) if 'SOURCE_GENE' in metacluster.SV_features else None
        INFO['STRAND'] = metacluster.SV_features['STRAND'] if 'STRAND' in metacluster.SV_features else None
        INFO['REGION'], INFO['GENE'] = metacluster.geneAnnot if hasattr(metacluster, 'geneAnnot') else (None, None)
        INFO['REP'] = ','.join([repeat['family'] for repeat in repeats]) if repeats else None 
        INFO['REPSUB'] = ','.join([repeat['subfamily'] for repeat in repeats]) if repeats else None   
        INFO['DIST'] = ','.join([str(repeat['distance']) for repeat in repeats]) if repeats else None
        INFO['NBTOTAL'], INFO['NBTUMOR'], INFO['NBNORMAL'] = str(metacluster.nbTotal), str(metacluster.nbTumour), str(metacluster.nbNormal) 
        INFO['NBSPAN'], INFO['NBCLIP'] = str(metacluster.nbINS), str(metacluster.nbCLIPPING)
        INFO['LEN'] = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        INFO['CV'] = metacluster.cv
        INFO['RTLEN'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['TRUN5LEN'] = metacluster.SV_features['TRUNCATION_5_LEN'] if 'TRUNCATION_5_LEN' in metacluster.SV_features else None
        INFO['TRUN3LEN'] = metacluster.SV_features['TRUNCATION_3_LEN'] if 'TRUNCATION_3_LEN' in metacluster.SV_features else None
        INFO['FULL'] = metacluster.SV_features['IS_FULL'] if 'IS_FULL' in metacluster.SV_features else None
        INFO['TDLEN'] = metacluster.SV_features['TRANSDUCTION_LEN'] if 'TRANSDUCTION_LEN' in metacluster.SV_features else None
        INFO['INVLEN'] = metacluster.SV_features['INVERSION_LEN'] if 'INVERSION_LEN' in metacluster.SV_features else None
        INFO['PERCR'] = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        INFO['QHITS'] = None if metacluster.insertHits is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['THITS'] = None if metacluster.insertHits is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['RTCOORD'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['POLYA'] = metacluster.SV_features['POLYA'] if 'POLYA' in metacluster.SV_features else None
        INFO['INSEQ'] = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None
        INFO['IDENTITY'] = metacluster.SV_features['IDENTITY'] if 'IDENTITY' in metacluster.SV_features else None
        INFO['SPECIDENTITY'] = metacluster.SV_features['SPECIDENTITY'] if 'SPECIDENTITY' in metacluster.SV_features else None
        #INFO['SPECIDENTITY'] = checkSpecIdent_LR(metacluster)
        INFO['SPECIDENTITY'] = None
        INFO['INS'] = ','.join(INS_readNames) if INS_readNames else None
        INFO['CLIP'] = ','.join(CLIP_readNames) if CLIP_readNames else None
        INFO['BKPCSEQR'] = metacluster.consLeftSeq if metacluster.consLeftSeq != None else None
        INFO['BKPCSEQL'] = metacluster.consLeftSeq if metacluster.consLeftSeq != None else None
        INFO['MAPQ'] = eventsMAPQMean

        ## Create VCF variant object
        fields = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, {}]

        ## Add variant to the VCF
        INS = formats.VCF_variant(fields)
        VCF.add(INS)
        
    ## 4. Sort VCF
    VCF.sort()

    ## 5. Write VCF in disk
    infoIds = ['VTYPE', 'ITYPE', 'MECHANISM', 'FAM', 'SUBFAM', 'GERMDB', 'CIPOS', 'CYTOID', \
           'NBEXONS', 'SRCGENE', 'STRAND', 'REGION', 'GENE', 'REP', 'REPSUB', 'DIST', \
           'NBTOTAL', 'NBTUMOR', 'NBNORMAL', 'NBSPAN', 'NBCLIP', 'LEN', 'CV', 'RTLEN', \
           'TRUN5LEN', 'TRUN3LEN', 'FULL', 'TDLEN', 'INVLEN', 'PERCR', \
           'QHITS', 'THITS', 'RTCOORD', 'POLYA', \
           'IDENTITY', 'SPECIDENTITY', 'MAPQ', 'INSEQ', 'INS', 'CLIP', 'BKPCSEQR', 'BKPCSEQL']

    VCF.write(infoIds, [], outName, outDir)


















#### Module sequences:

# UNUSED
def polish_sequence(events, name, outDir):
    '''
    Given clipping events of same clipping side, return polished clipping sequence.
    Input:
        1. events: List of clipping events of same clipping side.
        2. name: Uniq tag to name files (useful when running function in parallel).
        3. outDir: output directory.
    Output:
        1. polished_clippling_FASTA: FASTA file containing polished sequence or 'None'
    '''

    polished_FASTA = None
    
    # group keys with same values python:
    # https://stackoverflow.com/questions/54249400/python-how-to-group-keys-that-have-the-same-values-in-a-dictionary

    empiezos = {}
    for event in events:
        if event.clippedSide == 'right':
            if event.clipped_seq()[0:2] != '':
                empiezos[event.readName] = event.clipped_seq()[0:2]
        elif event.clippedSide == 'left':
            if event.clipped_seq()[-3:-1] != '':
                print ('event.readName ' + str(event.readName) + ' ' + str(event.clipped_seq()) + ' ' + str(event.clipped_seq()[-3:]))
                empiezos[event.readName] = event.clipped_seq()[-3:]

    print ('empiezos ' + str(empiezos))
    
    empiezosKeys = {}
    empiezosKeys = set(empiezos.values())

    print ('empiezosKeys ' + str(empiezosKeys))

    d = {}
    for n in empiezosKeys:
        d[n] = [k for k in empiezos.keys() if empiezos[k] == n]

    if d:
        # from list of list choose max length
        # https://stackoverflow.com/questions/53406024/finding-the-longest-and-the-shortest-lists-within-a-list-in-python
        print ('d ' + str(d))
        print ('max ' + str(max(d.values(), key=len)))
        maximo = max(d.values(), key=len)

        # cojo el read con la maxima length
        length = 0
        for event in events:
            if len(event.clipped_seq()) >= length:
                maxi = event
                length=len(event.clipped_seq())

        '''
        for event in events:
            if event.readName == maximo[0]:
                maxi = event
        '''


        polishedFastaEntireSequence = None

        # 1. Make fasta with template sequence
        # Pick first one of the list as template.
        templateFastaObj = formats.FASTA()
        templateFastaObj.seqDict[maxi.readName] = maxi.readSeq
        templateFastaPath = outDir + '/template_sequence_' + str(name) + '.fa'
        templateFastaObj.write(templateFastaPath)

        # 2. Make fasta with supporting reads
        supportingReadsFastaObj = formats.FASTA()
        for even in events:
            if even.readName != maxi.readName:

                supportingReadsFastaObj.seqDict[even.readName] = even.readSeq
        
        supportingReadsFastaPath = outDir + '/supporting_sequences_' + str(name) + '.fa'
        supportingReadsFastaObj.write(supportingReadsFastaPath)

        # 3. Polish with racon
        # TODO: change the following lines: 
        #polishedFasta = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, confDict['technology'], confDict['rounds'], outDir)
        polished_FASTA = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, 'NANOPORE', 1, outDir)

    return polished_FASTA

def polish_CLIP_sequence(metacluster, events, name, outDir):
    '''
    Given clipping events of same clipping side, return polished clipping sequence.
    Input:
        1. events: List of clipping events of same clipping side.
        2. name: Uniq tag to name files (useful when running function in parallel).
        3. outDir: output directory.
    Output:
        1. polished_clippling_FASTA: FASTA file containing polished sequence or 'None'
    '''

    polished_FASTA = None
    
    # group keys with same values python:
    # https://stackoverflow.com/questions/54249400/python-how-to-group-keys-that-have-the-same-values-in-a-dictionary

    empiezos = {}
    for event in events:
        if event.clippedSide == 'right':
            if event.clipped_seq()[0:2] != '':
                empiezos[event.readName] = event.clipped_seq()[0:2]
        elif event.clippedSide == 'left':
            if event.clipped_seq()[-3:-1] != '':
                print ('event.readName ' + str(event.readName) + ' ' + str(event.clipped_seq()) + ' ' + str(event.clipped_seq()[-3:]))
                empiezos[event.readName] = event.clipped_seq()[-3:]

    print ('empiezos ' + str(empiezos))
    
    empiezosKeys = {}
    empiezosKeys = set(empiezos.values())

    print ('empiezosKeys ' + str(empiezosKeys))

    d = {}
    for n in empiezosKeys:
        d[n] = [k for k in empiezos.keys() if empiezos[k] == n]

    if d:
        # from list of list choose max length
        # https://stackoverflow.com/questions/53406024/finding-the-longest-and-the-shortest-lists-within-a-list-in-python
        print ('d ' + str(d))
        print ('max ' + str(max(d.values(), key=len)))
        maximo = max(d.values(), key=len)

        # cojo el read con la maxima length
        length = 0
        for event in events:
            if not all ([event.supplementary for event in events]):
                if len(event.clipped_seq()) >= length and not event.supplementary:
                    maxi = event
                    length=len(event.clipped_seq())
            else:
                if len(event.clipped_seq()) >= length:
                    maxi = event
                    length=len(event.clipped_seq())

        '''
        for event in events:
            if event.readName == maximo[0]:
                maxi = event
        '''


        polishedFastaEntireSequence = None

        # 1. Make fasta with template sequence
        # Pick first one of the list as template.
        templateFastaObj = formats.FASTA()
        templateFastaObj.seqDict[maxi.readName + '_'+ str(metacluster.beg)] = maxi.clipped_seq()
        print ('templateFastaObj.seqDict ' +str(templateFastaObj.seqDict))
        templateFastaPath = outDir + '/template_sequence_' + str(name) + '.fa'
        templateFastaObj.write(templateFastaPath)

        # 2. Make fasta with supporting reads
        supportingReadsFastaObj = formats.FASTA()
        for even in events:
            if even.readName != maxi.readName:
                if even.clipped_seq():
                    supportingReadsFastaObj.seqDict[even.readName + '_'+ str(metacluster.beg)] = even.clipped_seq()
        
        print ('supportingReadsFastaObj.seqDict ' + str(supportingReadsFastaObj.seqDict))
        supportingReadsFastaPath = outDir + '/supporting_sequences_' + str(name) + '.fa'
        supportingReadsFastaObj.write(supportingReadsFastaPath)

        # 3. Polish with racon
        # TODO: change the following lines: 
        #polishedFasta = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, confDict['technology'], confDict['rounds'], outDir)
        polished_FASTA = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, 'NANOPORE', 1, outDir)
        print ('polished_FASTA ' + str(polished_FASTA))
        print ('outDir ' + str(outDir))

    return polished_FASTA

def polish_INS_sequence(metacluster, events, name, outDir):
    '''
    Given clipping events of same clipping side, return polished clipping sequence.
    Input:
        1. events: List of clipping events of same clipping side.
        2. name: Uniq tag to name files (useful when running function in parallel).
        3. outDir: output directory.
    Output:
        1. polished_clippling_FASTA: FASTA file containing polished sequence or 'None'
    '''
    polishedFastaEntireSequence = None

    # 1. Make fasta with template sequence
    # Pick first one of the list as template.
    templateFastaObj = formats.FASTA()
    if all([event.type == 'INS' for event in events]) or all([event.type == 'CLIPPING' for event in events]):
        print ('INS1')
        if events[0].type == 'INS':
            print ('INS1a')
            templateFastaObj.seqDict[events[0].readName + '_'+ str(metacluster.beg)] = events[0].pick_insert()
        if events[0].type == 'CLIPPING':
            print ('INS1b')
            if all([event.clippedSide == 'right' for event in events]) or all([event.clippedSide == 'left' for event in events]):
                print ('INS1c')
                templateFastaObj.seqDict[events[0].readName + '_'+ str(metacluster.beg)] = events[0].clipped_seq()

    # TODO: Limitation, if it has a lot of clippings but only one INS. Trying to fix this in last step.
    else:
        print ('INS2')
        for event in events:
            if event.type == 'INS':
                print ('INS3')
                templateFastaObj.seqDict[event.readName + '_'+ str(metacluster.beg)] = event.pick_insert()
                templateEvent=event.readName
                break
    templateFastaPath = outDir + '/template_sequence_' + str(name) + '.fa'
    templateFastaObj.write(templateFastaPath)

    '''
    if all([event.type == 'INS' for event in events]):
        templateFastaObj.seqDict[events[0].readName] = events[0].pick_insert()
    elif all([event.type == 'CLIPPING' for event in events]):
        # USual BND
        if all([event.clippedSide == 'right' for event in events]) or all([event.clippedSide == 'left' for event in events]):
            templateFastaObj.seqDict[events[0].readName] = events[0].clipped_seq()
        # INS_noOverlap
        
        # TODO: Solve this in another way
        else:
            polished_FASTAObj = formats.FASTA()
            for event in events:
                if event.identity != None:
                    eventRepre1 = event.readName
                    break
            for event in events:
                if event.identity != None and event.readName != eventRepre1:
                    eventRepre2 = event.readName
                    break
            if metacluster.consRightSeq:
                polished_FASTAObj.seqDict[eventRepre1]= ''.join(metacluster.consRightSeq)
            if metacluster.consLeftSeq:
                try:
                    polished_FASTAObj.seqDict[eventRepre2]= ''.join(metacluster.consLeftSeq)
                except UnboundLocalError:
                    polished_FASTAObj.seqDict[eventRepre1]= ''.join(metacluster.consLeftSeq)
            if metacluster.consRightSeq or metacluster.consLeftSeq:
                polished_FASTA = outDir + '/polished_1.fa'
                polished_FASTAObj.write(polished_FASTA)
                return polished_FASTA
        '''


    # 2. Make fasta with supporting reads
    supportingReadsFastaObj = formats.FASTA()

    if all([event.type == 'INS' for event in events]) or all([event.type == 'CLIPPING' for event in events]):
        for even in events[1:]:
            if even.type == 'INS':
                supportingReadsFastaObj.seqDict[even.readName + '_'+ str(metacluster.beg)] = even.pick_insert()
            elif even.type == 'CLIPPING':
                if all([event.clippedSide == 'right' for event in events]) or all([event.clippedSide == 'left' for event in events]):
                    supportingReadsFastaObj.seqDict[even.readName + '_'+ str(metacluster.beg)] = even.clipped_seq()
    else:
        for event in events:
            if event.type == 'INS' and event.readName != templateEvent:
                supportingReadsFastaObj.seqDict[event.readName + '_'+ str(metacluster.beg)] = event.pick_insert()
    
    supportingReadsFastaPath = outDir + '/supporting_sequences_' + str(name) + '.fa'
    supportingReadsFastaObj.write(supportingReadsFastaPath)


    # 3. Polish with racon
    # TODO: change the following lines: 
    #polishedFasta = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, confDict['technology'], confDict['rounds'], outDir)
    polished_FASTA = assembly.polish_racon(templateFastaPath, supportingReadsFastaPath, 'NANOPORE', 1, outDir)

    # If there is template but not supporting, return template.
    if not supportingReadsFastaObj.seqDict and not polished_FASTA and templateFastaObj.seqDict:
        polished_FASTA = templateFastaPath

    return polished_FASTA
