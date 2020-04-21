'''
Module to solve the bkp INGUAL DESPUES SE PUEDE UNIR!
'''
import databases
import formats
import os
import subprocess
import log
import unix
import operator
## [SR CHANGE]
import sequences
import assembly
import alignment
import events
import itertools
import bamtools


def analyzeMetaclusters(metaclusters, confDict, bam, normalBam, mode, outDir, binId, blatDbPath):
    '''
    SPANISH:
    1. Por cada evento en la binDB (en este caso especificamente por cada metacluster)
        a. Anadir supporting clipping al metacluster de discordant.
        b. Buscar el bkp (punto mas apoyado por los clippings que acabamos de anadir) y quitar aquellos clipping events que no lo soporten
        c. Hacer las cadenas de secuencias para ambos bkps.

    ENGLISH:
    1. For each event in the binDb (in this case, for each metacluster)
        a. Add supporting clipping events to the discordant metacluster.
        b. Look for reference breakpoint (most supported coordinate by clipping events added above) and remove those clipping events that dont support the bkp
        c. Make sequences of integrations for each bkp.
    
    Input:
        1. clustersBinDb
        2. confDict
        3. bam
        4. normalBam
        5. mode
        6. db: Fasta file db
        7. indexDb: Db indexed with minimap2
        8. outDir

    Output:
        1. dictMetaclusters: Nested dictionary -> First key: metacluster object. Second keys and corresponding values: refLeftBkp (value = None if UNK), refRightBkp (value = None if UNK), leftSeq (value = None if UNK), rightSeq (value = None if UNK), intLeftBkp (not present if UNK), intRightBkp (not present if UNK)
    '''
    metaclustersWODiscClip = {}

    # a. Add supporting clipping to discordant metacluster.

    for metacluster in metaclusters:

        # Make bkp directory
        bkpDir = outDir + '/' + metacluster.ref + '_' + str(metacluster.beg) + '_' + str(metacluster.end)
        unix.mkdir(bkpDir)

        if metacluster.orientation == 'PLUS':
            metaWODiscClip, clippingRightEventsToAdd, discClip = supportingCLIPPING(metacluster, 100, confDict, bam, normalBam, mode, bkpDir, 'PLUS')
            if metaWODiscClip and not discClip:
                metaclustersWODiscClip[metaWODiscClip] = clippingRightEventsToAdd

        elif metacluster.orientation == 'MINUS':
            metaWODiscClip, clippingLeftEventsToAdd, discClip = supportingCLIPPING(metacluster, 100, confDict, bam, normalBam, mode, bkpDir, 'MINUS')
            if metaWODiscClip and not discClip:
                metaclustersWODiscClip[metaWODiscClip] = clippingLeftEventsToAdd

        elif metacluster.orientation == 'RECIPROCAL':
            metaWODiscClip, clippingRightEventsToAdd, discClip = supportingCLIPPING(metacluster, 100, confDict, bam, normalBam, mode, bkpDir, 'PLUS')
            if metaWODiscClip and not discClip:
                metaclustersWODiscClip[metaWODiscClip] = clippingRightEventsToAdd
            metaWODiscClip, clippingLeftEventsToAdd, discClip = supportingCLIPPING(metacluster, 100, confDict, bam, normalBam, mode, bkpDir, 'MINUS')
            # Haces todo lo de abajo hasta add
            if metaWODiscClip and not discClip:
                if metaWODiscClip in  metaclustersWODiscClip.keys():
                    metaclustersWODiscClip[metaWODiscClip].extend(clippingLeftEventsToAdd)
                else:
                    metaclustersWODiscClip[metaWODiscClip] = clippingLeftEventsToAdd
    
    whenNoDiscClip(metaclustersWODiscClip, blatDbPath, binId, outDir)

    # b. Determine bkps.
    for metacluster in metaclusters:
        leftBkps = []
        rightBkps = []

        # Choose bkp with highest number of clipping events
        for event in metacluster.events:
            if event.type == 'CLIPPING':
                if event.clippedSide == 'left':
                    leftBkps.append(event.beg)
                elif event.clippedSide == 'right':
                    rightBkps.append(event.beg)

        if len(leftBkps) > 0:
            leftBkp = max(set(leftBkps), key=leftBkps.count)
            metacluster.refLeftBkp = leftBkp
        if len(rightBkps) > 0:
            rightBkp = max(set(rightBkps), key=rightBkps.count)
            metacluster.refRightBkp = rightBkp
        
        # Reconstruct bkp sequence:
        clipped_seqPlus = None
        clipped_seqMinus = None
        if metacluster.orientation == 'PLUS':
            clipped_seqPlus, consFastaBool = reconstructSeq(metacluster, confDict['consBkpSeq'], 'PLUS', outDir)
        elif metacluster.orientation == 'MINUS':
            clipped_seqMinus, consFastaBool = reconstructSeq(metacluster, confDict['consBkpSeq'], 'MINUS', outDir)
        elif metacluster.orientation == 'RECIPROCAL':
            clipped_seqPlus, consFastaBool = reconstructSeq(metacluster, confDict['consBkpSeq'], 'PLUS', outDir)
            clipped_seqMinus, consFastaBool = reconstructSeq(metacluster, confDict['consBkpSeq'], 'MINUS', outDir)
        
        # TODO: INT BKP -> simplmente contra la blat db

        if confDict['viralDb'] and metacluster.identity:
            if clipped_seqPlus != None:
                if not confDict['consBkpSeq'] or not consFastaBool:
                    # Create FASTA
                    clipped_seqPlusFastaPath = outDir + '/' + binId + '_clipped_seqPlusFastaPath.fasta'
                    clipped_seqPlusFasta = formats.FASTA()
                    clipped_seqPlusFasta.seqDict['clipped_seqPlusFasta'] = clipped_seqPlus
                    clipped_seqPlusFasta.write(clipped_seqPlusFastaPath)
                elif confDict['consBkpSeq']:
                    clipped_seqPlusFastaPath = clipped_seqPlus
                metacluster.intRightBkp = bkpINT(metacluster, clipped_seqPlusFastaPath, blatDbPath, bkpDir, metacluster.identity)
            if clipped_seqMinus != None:
                if not confDict['consBkpSeq'] or not consFastaBool:
                    # Create FASTA
                    clipped_seqMinusFastaPath = outDir + '/' + binId + '_clipped_seqMinusFastaPath.fasta'
                    clipped_seqMinusFasta = formats.FASTA()
                    clipped_seqMinusFasta.seqDict['clipped_seqMinusFasta'] = clipped_seqMinus
                    clipped_seqMinusFasta.write(clipped_seqMinusFastaPath)
                elif confDict['consBkpSeq']:
                    clipped_seqMinusFastaPath = clipped_seqMinus
                metacluster.intLeftBkp = bkpINT(metacluster, clipped_seqMinusFastaPath, blatDbPath, bkpDir, metacluster.identity)


    #clippings = [event for event in metacluster.events if event.type == 'CLIPPING']
    #discordants = [event for event in metacluster.events if event.type == 'DISCORDANT']

    #discClip = any(check in clippings for check in discordants)
        '''
        if metacluster.orientation == 'PLUS':
            # Collect discordant clipping left events
            plusDiscClip = {}
            for eventPlus in metacluster.events:
                if eventPlus.type == 'DISCORDANT':
                    lastOperation, lastOperationLen = eventPlus.cigarTuples[-1]
                    # Collect all clipping bkp genomic positions
                    if ((lastOperation == 4) or (lastOperation == 5)):
                        plusDiscClip[eventPlus.readName] = lastOperationLen
            
            # If there are discordant clipping left events      
            if plusDiscClip:
                # Collect clipping objects that are the same alignemten as discordant object:
                clippingsPlus = {}
                for eventPlus in metacluster.events:
                    if eventPlus.type == 'CLIPPING' and eventPlus.readName in plusDiscClip.keys():
                        clippingsPlus[eventPlus] = plusDiscClip[eventPlus.readName]
                
                # Choose the one with maximum clipping lenght
                if clippingsPlus: # This can be empty if clipping of discordant is too small.
                    largestClipping = max(clippingsPlus.items(), key=operator.itemgetter(1))[0]
                    # Make the representative sequence
                    metacluster.rightSeq = largestClipping.ref_seq() + '[INT]>' + largestClipping.clipped_seq()
            
            # If there are no discordant clipping left events  
            else:
                # Collect those clippings with blat hits
                clippingsPlus = {}
                for clippingPlus in metacluster.events:
                    if clippingPlus.type == 'CLIPPING':
                        if clippingPlus.blatIdentity == True:
                            clippingsPlus[clippingPlus] = clippingPlus.cigarTuples[-1][1]
                # If there are clippings with BLAT hits
                if clippingsPlus:
                    # Choose the one with maximum clipping lenght
                    largestClipping = max(clippingsPlus.items(), key=operator.itemgetter(1))[0]
                    # Make the representative sequence
                    metacluster.rightSeq = largestClipping.ref_seq() + '[INT]>' + largestClipping.clipped_seq()
                # If there are no clippings with BLAT hits
                else:
                    # Collect those clippings
                    clippingsPlus = {}
                    for clippingPlus in metacluster.events:
                        if clippingPlus.type == 'CLIPPING':
                            clippingsPlus[clippingPlus] = clippingPlus.cigarTuples[-1][1]
                    # If there are clippings
                    if clippingsPlus:
                        # Choose the one with maximum clipping length
                        largestClipping = max(clippingsPlus.items(), key=operator.itemgetter(1))[0]
                        # Make the representative sequence
                        # NOTE SR: this sequence will be less relayable
                        metacluster.rightSeq = largestClipping.ref_seq() + '[INT]>' + largestClipping.clipped_seq()
                    # If there are no clippings, there are not representative sequence.
                    else:
                        metacluster.rightSeq = None
        '''


def reconstructSeq(metacluster, consSeq, orientation, outDir):
    '''
    return clipped seq!
    '''
    # Collect discordant clipping left events
    discClip = {}
    clippingsDisc = {}
    for event in metacluster.events:
        if event.type == 'DISCORDANT':
            if orientation == 'PLUS':
                lastOperation, lastOperationLen = event.cigarTuples[-1]
                # Collect all clipping bkp genomic positions
                if ((lastOperation == 4) or (lastOperation == 5)):
                    discClip[event.readName] = lastOperationLen
            elif orientation == 'MINUS':
                firstOperation, firstOperationLen = event.cigarTuples[0]
                # Collect all clipping bkp genomic positions
                if ((firstOperation == 4) or (firstOperationLen == 5)):
                    discClip[event.readName] = firstOperationLen
    # If there are discordant clipping left events      
    if discClip:
        # Collect clipping objects that are the same alignemten as discordant object:
        clippingsDisc = {}
        for eventC in metacluster.events:
            if eventC.type == 'CLIPPING' and eventC.readName in discClip.keys():
                clippingsDisc[eventC] = discClip[eventC.readName]

    # Make representative sequence
    if clippingsDisc: # This can be empty if clipping of discordant is too small.
        if consSeq:
            clipped_seq, consFastaBool = conSeq(metacluster, clippingsDisc, orientation, outDir)
            return clipped_seq, consFastaBool
        elif not consSeq:
            clipped_seq = repreSeq(metacluster, orientation, clippingsDisc)
            return clipped_seq, False


    # If there are no discordant clipping left events  
    else:
        # Collect those clippings with blat hits
        clippingsBlat = {}
        for clippingB in metacluster.events:
            if clippingB.type == 'CLIPPING':
                if clippingB.blatIdentity == True:
                    clippingsBlat[clippingB] = clippingB.cigarTuples[-1][1]
        # If there are clippings with BLAT hits
        if clippingsBlat:
            if consSeq:
                clipped_seq, consFastaBool = conSeq(metacluster, clippingsBlat, orientation, outDir)
                return clipped_seq, consFastaBool

            elif not consSeq:
                clipped_seq = repreSeq(metacluster, orientation, clippingsBlat)
                return clipped_seq, False

        # If there are no clippings with BLAT hits
        else:
            # Collect metacluster clippings
            clippings = {}
            for clipping in metacluster.events:
                if clipping.type == 'CLIPPING':
                    clippings[clipping] = clipping.cigarTuples[-1][1]
            # If there are clippings
            if clippings:
                # NOTE SR: this sequence will be less relayable
                if consSeq:
                    clipped_seq, consFastaBool = conSeq(metacluster, clippings, orientation, outDir)
                    return clipped_seq, consFastaBool

                elif not consSeq:
                    clipped_seq = repreSeq(metacluster, orientation, clippingsBlat)
                    return clipped_seq, False

            # If there are no clippings, there are not representative sequence.
            else:
                if orientation == 'PLUS':
                    metacluster.rightSeq = None
                    return None, False
                elif orientation == 'MINUS':
                    metacluster.leftSeq = None
                    return None, False
    
    ### Do cleanup
    #unix.rm([bkpDir])

def repreSeq(metacluster, orientation, clippings):
    # Choose the one with maximum clipping lenght
    largestClipping = max(clippings.items(), key=operator.itemgetter(1))[0]
    # Make the representative sequence
    if orientation == 'PLUS':
        metacluster.rightSeq = largestClipping.ref_seq() + '[INT]>' + largestClipping.clipped_seq()
        return largestClipping.clipped_seq()
    elif orientation == 'MINUS':
        metacluster.leftSeq = largestClipping.clipped_seq() + '<[INT]' + largestClipping.ref_seq()
        return largestClipping.clipped_seq()

def conSeq(metacluster, clippings, orientation, outDir):
    # TODO SR: Luego poner el return de otra manera!!
    if len(clippings) > 1: 
        consFastaBool = False
        # Probar haciendo consenso COMPROBAR TIEMPO ETC list(itertools.chain(*clippingsDisc.keys()))
        if orientation == 'PLUS':
            metacluster.rightSeq, intConsensusPath = newConsensusSeq([*clippings], orientation, outDir)
            return intConsensusPath, consFastaBool
        elif orientation == 'MINUS':
            metacluster.leftSeq, intConsensusPath = newConsensusSeq([*clippings], orientation, outDir)
            return intConsensusPath, consFastaBool
    # If there are only one clipping, this will be the consensus.
    elif len(clippings) == 1:
        consFastaBool = True
        consClipping = [*clippings][0]
        if orientation == 'PLUS':
            metacluster.rightSeq = consClipping.ref_seq() + '[INT]>' + consClipping.clipped_seq()
            return consClipping.clipped_seq(), consFastaBool
        elif orientation == 'MINUS':
            metacluster.leftSeq = consClipping.clipped_seq() + '<[INT]' + consClipping.ref_seq()
            return consClipping.clipped_seq(), consFastaBool

def newConsensusSeq(clippings, orientation, outDir):
    refConsensusSeq = makeConsSeqs(clippings, 'REF', outDir)[1]
    intConsensusPath, intConsensusSeq = makeConsSeqs(clippings, 'INT', outDir)
    consSeq = None

    if intConsensusSeq != None:
        if orientation == 'PLUS':
            consSeq = refConsensusSeq + '[INT]>' + intConsensusSeq
        elif orientation == 'MINUS':
            consSeq = intConsensusSeq + '<[INT]' + refConsensusSeq

    #return consSeq, intConsensusPath
    # TEMP
    return consSeq, intConsensusSeq

def whenNoDiscClip(metaclustersWODiscClip, db, binId, outDir):
    clippingEventsToAdd = list(itertools.chain(*metaclustersWODiscClip.values()))
    clippingsFasta = writeClippingsFasta(clippingEventsToAdd, binId, outDir)

    # BLAT
    outName = binId + '_clippingsBlat'
    ## 3. Align clipped sequences with BLAT againsst db
    pslPath = alignment.alignment_blat(clippingsFasta, db, outName, outDir)
    pslDict = makePslDict(pslPath)
    matchClippings = collectMatchClippings(metaclustersWODiscClip, pslDict)
    if matchClippings != {}:
        clippings2Add = collectClipBkpMatch(matchClippings, clippingEventsToAdd)
        if clippings2Add:
            for metacluster, clippings in clippings2Add:
                metacluster.addEvents(clippings)



def makePslDict(pslPath):
    # Read PSL
    pslClipping = formats.PSL()
    pslClipping.read(pslPath)
    ## TODO SR: mirar filtros
    pslDict = {}
    for pslAlign in pslClipping.alignments:
        if pslAlign.qName in pslDict.keys():
            pslDict[pslAlign.qName].append(pslAlign.tName)
        else:
            pslDict[pslAlign.qName] = []
            pslDict[pslAlign.qName].append(pslAlign.tName)
    return pslDict
    

def collectMatchClippings(metaclustersWODiscClip, pslDict):
    matchClippings = {}
    for metaclusterWODiscClip, clippings in metaclustersWODiscClip.items():
        for clip in clippings:
            if clip.readName in pslDict.keys():
                if metaclusterWODiscClip.identity not in pslDict[clip.readName]:
                    if metaclusterWODiscClip in matchClippings.keys():
                        matchClippings[metaclusterWODiscClip].append(clip)
                    else:
                        matchClippings[metaclusterWODiscClip] = []
                        matchClippings[metaclusterWODiscClip].append(clip)

    return matchClippings
                
def collectClipBkpMatch(matchClippings, clippingEventsToAdd):
    
    clippings2Add = {}
    for metacluster, clippings in matchClippings.items():
        # Collect bkp position of BLAT hits clippping events
        bkp = []
        for matchClipping in clippings:
            bkp.append(matchClipping.readBkp)

        # Make coordinates of the region to collect those clipping events:
        binBeg = min(bkp) - 5 if min(bkp) >= 5 else 0
        binEnd = max(bkp) + 5

        # Collect clipping events whose bkp is in desired region
        for clipping in clippingEventsToAdd:
            if clipping.readBkp >= binBeg and clipping.readBkp <= binEnd:
                if metacluster in clippings2Add.keys():
                    clippings2Add[metacluster].append(clipping)
                    clipping.blatIdentity = True
                else:
                    clippings2Add[metacluster] = []
                    clippings2Add[metacluster].append(clipping)
                    clipping.blatIdentity = True
    
    return clippings2Add
                    


'''
# c. Make sequences of integrations for each bkp.
#dictMetaclusters[metacluster]['leftSeq'], dictMetaclusters[metacluster]['rightSeq'] = makeConsSeqs(CLIPPING_cluster, 'REF', db, indexDb, bkpDir)
    
    leftIntConsensusSeq = None
    rightIntConsensusSeq = None

    leftRefConsensusSeq = makeConsSeqs(allEvents, 'left', 'REF', bkpDir)[1]
    rightRefConsensusSeq = makeConsSeqs(allEvents, 'right', 'REF', bkpDir)[1]

    leftIntConsensusPath, leftIntConsensusSeq = makeConsSeqs(allEvents, 'left', 'INT', bkpDir)
    rightIntConsensusPath, rightIntConsensusSeq = makeConsSeqs(allEvents, 'right', 'INT', bkpDir)

    if leftIntConsensusSeq != None:
        leftSeq = leftIntConsensusSeq + '<[INT]' + leftRefConsensusSeq
    else:
        leftSeq = None
    if rightIntConsensusSeq != None:
        rightSeq = rightRefConsensusSeq + '[INT]>' + rightIntConsensusSeq
    else:
        rightSeq = None

    metacluster.leftSeq, metacluster.rightSeq = leftSeq, rightSeq

    if confDict['viralDb']:
        if leftIntConsensusPath != None:
            metacluster.intLeftBkp = bkpINT(metacluster, leftIntConsensusPath, confDict['viralDb'], bkpDir, metacluster.identity)
        if rightIntConsensusPath != None:
            metacluster.intRightBkp = bkpINT(metacluster, rightIntConsensusPath, confDict['viralDb'], bkpDir, metacluster.identity)
'''

        #else:
            #dictMetaclusters[metacluster]['refLeftBkp'], dictMetaclusters[metacluster]['refRightBkp'], dictMetaclusters[metacluster]['leftSeq'], dictMetaclusters[metacluster]['rightSeq'], dictMetaclusters[metacluster]['intLeftBkp'], dictMetaclusters[metacluster]['intRightBkp'] = None, None, None, None, None, None


def supportingCLIPPING(metacluster, buffer, confDict, bam, normalBam, mode, outDir, side):
    # Note: This function works but you have to allow duplicates in the clipping 

    clippingRightEventsToAdd = {}
    clippingLeftEventsToAdd = {}
    clippingConfDict = dict(confDict)
    clippingRightEventsToAddBlat = {}
    clippingLeftEventsToAddBlat = {}
    clippingConfDict['targetSV'] = ['CLIPPING']
    clippingConfDict['minMAPQ'] = 10
    # NOTE SR: Think and check if this is neccessary. I think it is not
    #confDict['minCLIPPINGlen'] = 2
    clippingEventsDict = {}
    CLIPPING_clusters = []

    ## Define region
    if metacluster.orientation != 'RECIPROCAL':
        if metacluster.orientation == 'PLUS':
            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determinePlusBkpArea(metacluster.beg, metacluster.end, metacluster.events, buffer)

        elif metacluster.orientation == 'MINUS':
            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determineMinusBkpArea(metacluster.beg, metacluster.end, metacluster.events, buffer)

        # Collect clippings
        ref = metacluster.ref

        if mode == "SINGLE":
            clippingEventsDict = bamtools.collectSV(ref, binBeg, binEnd, bam, clippingConfDict, None, False)

        elif mode == "PAIRED":
            clippingEventsDict = bamtools.collectSV_paired(ref, binBeg, binEnd, bam, normalBam, clippingConfDict)

        ## When the discordant cluster is RIGHT, add the biggest right clipping cluster if any:
        #if all (event.orientation == 'PLUS' for event in self.events):
        # When the metacluster is RIGHT:
        if metacluster.orientation == 'PLUS':
            if 'RIGHT-CLIPPING' in clippingEventsDict.keys():
                # Keep only those clipping event which have their clipping bkp is in the desired area
                clippingRightEventsToAdd = chooseBkpClippings(clippingEventsDict, 'RIGHT-CLIPPING', binBeg, binEnd)
                clippings = list(itertools.chain(*clippingRightEventsToAdd.values()))
                

                # If there are discordants with clipping, this is enough because the area of adding clippings is narrow:
                # TODO SR: For the moment, no BLAT id performed for MEs insertions.
                if discClip == True or not confDict['viralDb'] or not metacluster.identity:
                    # Make clipping clusters and add events to metacluster
                    # TODO SR: Think if making clusters is neccessary. Poner lo siguiente en vez de lo que ya hya!!
                    metacluster.addEvents(clippings)
                    #CLIPPING_clusters = metacluster.add_clippingEvents(ref, binBeg, binEnd, clippingRightEventsToAdd, ['RIGHT-CLIPPING'], confDict)
                    return metacluster, clippings, discClip
                elif metacluster.identity:
                    return metacluster, clippings, discClip
                else:
                    return None, None, None
                '''
                # If there are No discordants with clipping, the area of adding clippings is wide, so we have to select them performing BLAT aligment:
                elif discClip == False and confDict['viralDb'] and metacluster.identity:
                    clippingsFasta = bkp.writeClippingsFasta(clippingRightEventsToAdd, 'RIGHT-CLIPPING')
                    metacluster.clippingEvents = clippingRightEventsToAdd['RIGHT-CLIPPING']
                    #clippingRightEventsToAddBlat = bkp.blatClippings(clippingRightEventsToAdd, 'RIGHT-CLIPPING', self.identity, self.id, confDict['viralDb'], outDir)
                    #CLIPPING_clusters = self.add_clippingEvents(ref, binBeg, binEnd, clippingRightEventsToAddBlat, ['RIGHT-CLIPPING'], confDict)
                '''

        ## When the discordant cluster is LEFT, add the biggest left clipping cluster if any:
        #elif all (event.orientation == 'MINUS' for event in self.events):
        # When the metacluster is MINUS:
        elif metacluster.orientation == 'MINUS':
            if 'LEFT-CLIPPING' in clippingEventsDict.keys():
                # Keep only those clipping event which have their clipping bkp is in the desired area
                clippingLeftEventsToAdd = chooseBkpClippings(clippingEventsDict, 'LEFT-CLIPPING', binBeg, binEnd)
                clippings = list(itertools.chain(*clippingLeftEventsToAdd.values()))

                # If there are discordants with clipping:
                # TODO SR: For the moment, no BLAT id performed for MEs insertions.
                if discClip == True or not confDict['viralDb'] or not metacluster.identity:
                    # Make clipping clusters and add events to metacluster
                    # TODO SR: Think if making clusters is neccessary
                    metacluster.addEvents(clippings) 
                    #CLIPPING_clusters = metacluster.add_clippingEvents(ref, binBeg, binEnd, clippingLeftEventsToAdd, ['LEFT-CLIPPING'], confDict)
                    return metacluster, clippings, discClip
                elif metacluster.identity:
                    return metacluster, clippings, discClip
                else:
                    return None, None, None
                '''
                # If there are No discordants with clipping -> BLAT:
                elif discClip == False and confDict['viralDb'] and metacluster.identity:
                    clippingLeftEventsToAddBlat = bkp.blatClippings(clippingLeftEventsToAdd, 'LEFT-CLIPPING', metacluster.identity, metacluster.id, confDict['viralDb'], outDir)
                    CLIPPING_clusters = metacluster.add_clippingEvents(ref, binBeg, binEnd, clippingLeftEventsToAddBlat, ['LEFT-CLIPPING'], confDict)
                '''

    # When the metacluster is RECIPROCAL:
    # Fisrt separate PLUS AND MINUS events and them do the same as above 
    elif metacluster.orientation == 'RECIPROCAL':
        if side == 'PLUS':
            clippingPlusEventsDict = {}
            clippingMinusEventsDict = {}

            # Collect PLUS discordant events:
            reciprocalPlusEvents = [eventP for eventP in metacluster.events if eventP.orientation == 'PLUS']

            # HAcer el beg y el end en funcion solo de los PLUS!
            begPlus = min([eventPlus.beg for eventPlus in reciprocalPlusEvents])
            endPlus = max([eventPlus.end for eventPlus in reciprocalPlusEvents])
            #print ('reciprocalPlusEvents ' + str(reciprocalPlusEvents))
            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBegP, binEndP, discClipPlus = determinePlusBkpArea(begPlus, endPlus, reciprocalPlusEvents, buffer)

            # Collect right clippings
            ref = metacluster.ref

            if mode == "SINGLE":
                clippingPlusEventsDict = bamtools.collectSV(ref, binBegP, binEndP, bam, clippingConfDict, None, False)
            elif mode == "PAIRED":
                clippingPlusEventsDict = bamtools.collectSV_paired(ref, binBegP, binEndP, bam, normalBam, clippingConfDict)
            
            if 'RIGHT-CLIPPING' in clippingPlusEventsDict.keys():

                #print ('clippingPlusEventsDict ' + str(clippingPlusEventsDict))
                # Keep only those clipping event which have their clipping bkp is in the desired area
                clippingRightEventsToAdd = chooseBkpClippings(clippingPlusEventsDict, 'RIGHT-CLIPPING', binBegP, binEndP)
                clippings = list(itertools.chain(*clippingRightEventsToAdd.values()))

                #print ('clippingRightEventsToAdd ' + str(clippingRightEventsToAdd))

                # If there are No discordants with clipping, the area of adding clippings is wide, so we have to select them performing BLAT aligment.
                # Also, delete prrevious dictionary to avoid add the clippings of the "narrow" region to final dictionary
                # TODO SR: For the moment, no BLAT id performed for MEs insertions.

                if discClipPlus == True or not confDict['viralDb'] or not metacluster.identity:
                    metacluster.addEvents(clippings) 
                    #clippingRightEventsToAddBlat = bkp.blatClippings(clippingRightEventsToAdd, 'RIGHT-CLIPPING', metacluster.identity, metacluster.id, confDict['viralDb'], outDir)
                    #clippingRightEventsToAdd = {}
                    return metacluster, clippings, discClipPlus
                elif metacluster.identity:
                    return metacluster, clippings, discClipPlus
                else:
                    return None, None, None

        if side == 'MINUS':
            # Collect MINUS discordant events:
            reciprocalMinusEvents = []
            for eventM in metacluster.events:
                if eventM.type == 'DISCORDANT':
                    if eventM.orientation == 'MINUS':
                        reciprocalMinusEvents.append(eventM)


            begMinus = min([eventMinus.beg for eventMinus in reciprocalMinusEvents])
            endPlus = max([eventMinus.end for eventMinus in reciprocalMinusEvents])

            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBegM, binEndM, discClipMinus = determineMinusBkpArea(begMinus, endPlus, reciprocalMinusEvents, buffer)

            # Collect left clippings
            ref = metacluster.ref

            if mode == "SINGLE":
                clippingMinusEventsDict = bamtools.collectSV(ref, binBegM, binEndM, bam, clippingConfDict, None, False)
            elif mode == "PAIRED":
                clippingMinusEventsDict = bamtools.collectSV_paired(ref, binBegM, binEndM, bam, normalBam, clippingConfDict)
            
            if 'LEFT-CLIPPING' in clippingMinusEventsDict.keys():

                #print ('clippingMinusEventsDict ' + str(clippingMinusEventsDict))
                # Keep only those clipping event which have their clipping bkp is in the desired area
                clippingLeftEventsToAdd = chooseBkpClippings(clippingMinusEventsDict, 'LEFT-CLIPPING', binBegM, binEndM)
                clippings = list(itertools.chain(*clippingLeftEventsToAdd.values()))

                # If there are No discordants with clipping, the area of adding clippings is wide, so we have to select them performing BLAT aligment.
                # Also, delete prrevious dictionary to avoid add the clippings of the "narrow" region to final dictionary
                # TODO SR: For the moment, no BLAT id performed for MEs insertions.
                if discClipMinus == True and confDict['viralDb'] or not metacluster.identity:
                    metacluster.addEvents(clippings) 
                    #clippingLeftEventsToAddBlat = bkp.blatClippings(clippingLeftEventsToAdd, 'LEFT-CLIPPING', metacluster.identity, metacluster.id, confDict['viralDb'], outDir)
                    #clippingLeftEventsToAdd = {}
                    return metacluster, clippings, discClipMinus
                elif metacluster.identity:
                    return metacluster, clippings, discClipMinus
                else:
                    return None, None, None

                #print ('clippingLeftEventsToAdd ' + str(clippingLeftEventsToAdd))
                    '''
                    # Merge all clipping events to add in only one dictionary.
                    clippingRightEventsToAdd.update(clippingLeftEventsToAdd)
                    clippingRightEventsToAdd.update(clippingRightEventsToAddBlat)
                    clippingRightEventsToAdd.update(clippingLeftEventsToAddBlat)
                    clippingEventsDict = clippingRightEventsToAdd

                    #print ('clippingEventsDict ' + str(clippingEventsDict))

                    # Make clipping clusters and add them to metacluster:
                    if clippingEventsDict:
                        metacluster.addEvents(clippingEventsDict.values())
                        #CLIPPING_clusters = metacluster.add_clippingEvents(ref, binBegP, binEndM, clippingEventsDict, ['RIGHT-CLIPPING', 'LEFT-CLIPPING'], confDict)
                        return None, None
                    elif metacluster.identity:
                        return metacluster, clippingRightEventsToAdd
                    else:
                        return None, None
                    '''
# TODO SR: UNUSED FUNCTION clippingBkp
def clippingBkp(CLIPPING_clusters):
    '''
    Look for reference breakpoint (most supported coordinate by clipping events added above) and remove those clipping events that dont support the bkp

    Input:
        1. CLIPPING_cluster
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

    # TODO: AQUI ES MEJOR QUITARLO DE LA LISTA QUE HACER UNA LISTA NUEVA, PERO DE MOMENTO VA ASI!
    newEvents = []
    # Remove those clipping events that dont support the bkp
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


def makeConsSeqs(clippingEvents, seqSide, outDir):
    '''
    Make consesus sequence of one of the sides of the clipping cluster.

    Input:
        1. CLIPPING_cluster
        2. clippedSide: 'left' or 'right'
        3. seqSide: 'INT' if the consensus of the integrated sequence is wanted or 'REF' if the consensus of the reference is wanted
    
    Output:
        1. consensusPath: consensus file
        2. consensusSeq: consensus sequence
    '''
    consensusPath = None
    consensusSeq = None

    #clippingEvents = [event for event in CLIPPING_cluster.events if event.clippedSide == clippedSide]

    if len (clippingEvents) > 0:

        consensusPath, consensusSeq = clippingConsensusSeq(clippingEvents, clippingEvents[0].id, seqSide, outDir)
    
    return consensusPath, consensusSeq


def clippingConsensusSeq(clippingEvents, CLIPPING_clusterID, seqSide, outDir):

    # Retrieve fasta file with sequence from match or clipped side of clipping reads
    supportingReadsFasta = clippingSeq(clippingEvents, CLIPPING_clusterID, seqSide, outDir)

    # Consensus from the previous fasta
    consensusPath, consensusSeq = assembly.getConsensusSeq(supportingReadsFasta, outDir)

    return consensusPath, consensusSeq


def clippingSeq(clippingEvents, CLIPPING_clusterID, seqSide, outDir):
    '''
    Retrieve fasta file with sequence from match or clipped side of clipping reads
    '''

    fastaObj = formats.FASTA()
    fastaDict = {}
    # Determine bkp
    for event in clippingEvents:
        # Si queremos sacar la secuencia del lado de la integracion:
        if seqSide == 'INT':
            fastaDict[event.readName] = event.clipped_seq()

        elif seqSide == 'REF':
            fastaDict[event.readName] = event.ref_seq()
        
        fastaObj.seqDict = fastaDict

    fastaPath = outDir + '/' + str(CLIPPING_clusterID) +'_'+ str(seqSide) +'_supportingReads.fa'
    fastaObj.write(fastaPath)

    return fastaPath


def bkpINT(metacluster, consensusPath, db, outDir, identity):

    #indexDbSpecificIdentity = databases.buildIdentityDb(metacluster, db, outDir)   

    PAF_file = sequences.getPAFAlign(consensusPath, db, outDir)
    PAFObj = formats.PAF()
    PAFObj.read(PAF_file)

    intBkp = {}
    if not os.stat(PAF_file).st_size == 0:
        # DE AQUI SACAMOS LA INFO QUE QUERAMOS DEL VIRUS
        for alig in PAFObj.alignments:
            # Check identity
            if alig.tName.split('|')[0] in identity:
                intBkp[alig.tName.split('|')[1]] = alig.tBeg
        #intBkp = [line.tBeg for line in PAFObj.alignments][0]
    else:
        intBkp = {}

    return intBkp

def determinePlusBkpArea(beg, end, events, buffer):
    '''
    This function determines the coordinates from where clipping events should be added to a PLUS metacluster. Also, it determines if there are discordant clipping or not. 
    If there are discordant clipping, the region to look for clippings will be [the most left bkp -2 : the most left bkp +2].
    Otherwise, the the region to look for clipping events is from metacluster.beg - buffer to metacluster.end + buffer.

    Input:
        1. beg: beg position
        2. end: end position
        3. events: list of events
        4. buffer: bp of region extension when there are no discordant clippings.
    Ouput:
        1. binBeg: beg position of the region to look for clipping events.
        2. binEnd: end position of the region to look for clipping events.
        3. discClip: bool. True is there are discordatn events, False otherwise.
    '''

    plusBkp = []
    # If there are clippings in discordant events: Collect all clipping bkp genomic positions
    for discordantPlus in events:
        lastOperation, lastOperationLen = discordantPlus.cigarTuples[-1]
        # Collect all clipping bkp genomic positions
        if ((lastOperation == 4) or (lastOperation == 5)):
            plusBkp.append(discordantPlus.end)

    # The region to look for clippings will be [the most left bkp -2 : the most left bkp +2]
    if len(plusBkp) > 0:
        binBeg = min(plusBkp) - 5 if min(plusBkp) >= 5 else 0
        binEnd = max(plusBkp) + 5
        discClip = True

    
    # If there are NOT clippings in discordant events
    # The region to look for clippings will be metacluster positions + buffer
    else:
        binBeg = beg
        binEnd = end + buffer
        discClip = False

    return binBeg, binEnd, discClip

def determineMinusBkpArea(beg, end, events, buffer):
    '''
    This function determines the coordinates from where clipping events should be added to a MINUS metacluster. Also, it determines if there are discordant clipping or not. 
    If there are discordant clipping, the region to look for clippings will be [the most left bkp -2 : the most left bkp +2].
    Otherwise, the the region to look for clipping events is from metacluster.beg - buffer to metacluster.end + buffer.

    Input:
        1. beg: beg position
        2. end: end position
        3. events: list of events
        4. buffer: bp of region extension when there are no discordant clippings.

    Ouput:
        1. binBeg: beg position of the region to look for clipping events.
        2. binEnd: end position of the region to look for clipping events.
        3. discClip: bool. True is there are discordatn events, False otherwise.
    '''
    minusBkp = []
    # If there are clippings in discordant events: Collect all clipping bkp genomic positions
    for discordantMinus in events:
        firstOperation, firstOperationLen = discordantMinus.cigarTuples[0]
        # Collect all clipping bkp genomic positions
        if ((firstOperation == 4) or (firstOperationLen == 5)):
            minusBkp.append(discordantMinus.beg)

    # The region to look for clippings will be [the most left bkp -2 : the most left bkp +2]
    if len(minusBkp) > 0:
        binBeg = min(minusBkp) - 5 if min(minusBkp) >= 5 else 0
        binEnd = max(minusBkp) + 5
        discClip = True
    
    # If there are NOT clippings in discordant events
    # The region to look for clippings will be metacluster positions + buffer
    else:
        binBeg = beg - buffer if beg > buffer else 0
        # TODO SR: check as in determinePlusBkpArea
        binEnd = end
        discClip = False
    
    return binBeg, binEnd, discClip

def chooseBkpClippings(clippingEventsDict, eventType, binBeg, binEnd):
    '''
    Keep only those clipping events which clipping bkp is between input coordinates.

    Input:
        1. clippingEventsDict: clippingEventsDict[eventType] = []
        2. eventType: 'RIGHT-CLIPPING' or 'LEFT-CLIPPING'
        3. binBeg: beg position
        4. binEnd: end position
    Output:
        1. clippingEventsDict: same structure as input dictionary, containing only those clipping events which clipping bkp is between input coordinates.
    '''
    # eventType = 'RIGHT-CLIPPING' or 'LEFT-CLIPPING'
    clippingEventsToAdd = {}
    clippingEventsToAdd[eventType] = []
    
    ## Get clipping clusters:
    #clippingEventsDict = dict((key,value) for key, value in clippingEventsDict.items() if key == eventType)
    for clippingEvent in clippingEventsDict[eventType]:
        if eventType == 'RIGHT-CLIPPING':
            if binBeg <= clippingEvent.end <= binEnd:
                clippingEventsToAdd[eventType].append(clippingEvent)
        elif eventType == 'LEFT-CLIPPING':
            if binBeg <= clippingEvent.beg <= binEnd:
                clippingEventsToAdd[eventType].append(clippingEvent)

    return clippingEventsDict

def writeClippingsFasta(clippings, binId, outDir):
    # TODO SR: POner binid!!!!!!!
    '''
    Perform a BLAT search with clipping events, select bkp area of those that match and pick all clipping event which bkp is in this area.

    Input:
        1. clippingEventsDict: clippingEventsDict[eventType] = []
        2. eventType: 'RIGHT-CLIPPING' or 'LEFT-CLIPPING'
        3. identity: metacluster identity
        4. ID: metacluster ID
        5. db: reference DB
        6. outDir: output directory
    '''

    ## 1. Generate fasta containing soft clipped sequences
    clippedFasta = events.collect_clipped_seqs(clippings)

    ## 2. Write clipped sequences into fasta
    filePath = outDir + '/' + binId + '_clippings.fasta'
    clippedFasta.write(filePath, 'append', False)

    return filePath


def blatClippings(clippingsDict, eventType, identity, ID, db, outDir):

    '''
    Perform a BLAT search with clipping events, select bkp area of those that match and pick all clipping event which bkp is in this area.

    Input:
        1. clippingEventsDict: clippingEventsDict[eventType] = []
        2. eventType: 'RIGHT-CLIPPING' or 'LEFT-CLIPPING'
        3. identity: metacluster identity
        4. ID: metacluster ID
        5. db: reference DB
        6. outDir: output directory
    
    Output:
        1. clippings2Add: dictionary with clipping events to add to the metacluster.
    '''

    # Select those clippings that match with same sequences as identity metacluster.
    matchClippings = []
    for pslAlign in pslClipping.alignments:
        if identity in pslAlign.tName:
            matchClippings.append([clipping for clipping in clippingsDict[eventType] if clipping.fullReadName() == pslAlign.qName])
    
    matchClippingsUniq = list(set(matchClippings))

    clippings2Add = {}
    clippings2Add[eventType] = []

    if matchClippingsUniq != []:
        # Collect bkp position of BLAT hits clippping events
        bkp = []
        for matchClipping in matchClippingsUniq:
            bkp.append(matchClipping.readBkp)

        # Make coordinates of the region to collect those clipping events:
        binBeg = min(bkp) - 5 if min(bkp) >= 5 else 0
        binEnd = max(bkp) + 5

        # Collect clipping events whose bkp is in desired region
        for clipping in clippingsDict[eventType]:
            if clipping.readBkp >= binBeg and clipping.readBkp <= binEnd:
                clippings2Add[eventType].append(clipping)

    return clippings2Add