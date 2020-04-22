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
    Four main steps
        a. Add supporting clipping events to discordant metacluster.
        b. Determine metacluster bkps.
        c. Reconstruct bkp sequence.
        d. Determine inserted sequence bkp.
    
    Input
        1. metaclusters: List of metaclusters
        2. confDict
        3. bam
        4. normalBam
        5. mode
        6. outDir
        7. binId
        8. blatDbPath
    
    Output
        Doesn't return anything, just fill some metaclusters attributes:
        1. metacluster.refRightBkp
        2. metacluster.refLeftBkp
        3. metacluster.rightSeq
        4. metacluster.leftSeq
        5. metacluster.intRightBkp
        6. metacluster.intLeftBkp
    '''

    # Initiliaze dictionary: metaclustersWODiscClip[metacluster] = clippings -> key: metacluster without discordant clipping events; value: list of candidate clippings events to add.
    metaclustersWODiscClip = {}

    # a. Add supporting clipping to discordant metacluster.
    
    for metacluster in metaclusters:

        # Make bkp directory
        bkpDir = outDir + '/' + metacluster.ref + '_' + str(metacluster.beg) + '_' + str(metacluster.end)
        unix.mkdir(bkpDir)

        if metacluster.orientation != 'RECIPROCAL':
            metaWODiscClip, clippingEventsToAdd, discClip = supportingCLIPPING(metacluster, 100, confDict, bam, normalBam, mode, bkpDir, metacluster.orientation)
            # TODO SR: For the moment, no BLAT id performed for MEs insertions.
            # Add clippings if there are discodant clippings.
            if discClip == True or not confDict['viralDb'] or not metacluster.identity:
                metacluster.addEvents(clippingEventsToAdd)
            # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
            elif metaWODiscClip and not discClip and metacluster.identity:
                metaclustersWODiscClip[metaWODiscClip] = clippingEventsToAdd

        elif metacluster.orientation == 'RECIPROCAL':
            metaWODiscClip, clippingRightEventsToAdd, discClip = supportingCLIPPING(metacluster, 100, confDict, bam, normalBam, mode, bkpDir, 'PLUS')
            # Add clippings if there are discodant clippings.
            if discClip == True or not confDict['viralDb'] or not metacluster.identity:
                metacluster.addEvents(clippingRightEventsToAdd)
            # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
            elif metaWODiscClip and not discClip and metacluster.identity:
                metaclustersWODiscClip[metaWODiscClip] = clippingRightEventsToAdd
            
            metaWODiscClip, clippingLeftEventsToAdd, discClip = supportingCLIPPING(metacluster, 100, confDict, bam, normalBam, mode, bkpDir, 'MINUS')
            # Add clippings if there are discodant clippings.
            if discClip == True or not confDict['viralDb'] or not metacluster.identity:
                metacluster.addEvents(clippingLeftEventsToAdd)
            # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
            elif metaWODiscClip and not discClip and metacluster.identity:
                if metaWODiscClip in  metaclustersWODiscClip.keys():
                    metaclustersWODiscClip[metaWODiscClip].extend(clippingLeftEventsToAdd)
                else:
                    metaclustersWODiscClip[metaWODiscClip] = clippingLeftEventsToAdd
    
    # Add clippings when there are no discordant clippings, but they have BLAT matches and clippings without blat hits but same bkp as the ones that match
    addBlatClippings(metaclustersWODiscClip, blatDbPath, binId, outDir)

    
    for metacluster in metaclusters:
        # b. Determine bkps.

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
        
        # c. Reconstruct bkp sequence:

        clipped_seqPlus = None
        clipped_seqMinus = None
        if metacluster.orientation == 'PLUS':
            clipped_seqPlus, consFastaBoolPlus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'PLUS', outDir)
        elif metacluster.orientation == 'MINUS':
            clipped_seqMinus, consFastaBoolMinus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'MINUS', outDir)
        elif metacluster.orientation == 'RECIPROCAL':
            clipped_seqPlus, consFastaBoolPlus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'PLUS', outDir)
            clipped_seqMinus, consFastaBoolMinus = reconstructSeq(metacluster, confDict['consBkpSeq'], 'MINUS', outDir)
        
        # d. Determine inserted sequence bkp

        if confDict['viralDb'] and metacluster.identity:
            if clipped_seqPlus != None:
                if not confDict['consBkpSeq'] or not consFastaBoolPlus:
                    # Create FASTA
                    clipped_seqPlusFastaPath = outDir + '/' + binId + '_clipped_seqPlusFastaPath.fasta'
                    clipped_seqPlusFasta = formats.FASTA()
                    clipped_seqPlusFasta.seqDict['clipped_seqPlusFasta'] = clipped_seqPlus
                    clipped_seqPlusFasta.write(clipped_seqPlusFastaPath)
                elif confDict['consBkpSeq']:
                    clipped_seqPlusFastaPath = clipped_seqPlus
                metacluster.intRightBkp = bkpINT(metacluster, clipped_seqPlusFastaPath, blatDbPath, bkpDir, metacluster.identity)
            if clipped_seqMinus != None:
                if not confDict['consBkpSeq'] or not consFastaBoolMinus:
                    # Create FASTA
                    clipped_seqMinusFastaPath = outDir + '/' + binId + '_clipped_seqMinusFastaPath.fasta'
                    clipped_seqMinusFasta = formats.FASTA()
                    clipped_seqMinusFasta.seqDict['clipped_seqMinusFasta'] = clipped_seqMinus
                    clipped_seqMinusFasta.write(clipped_seqMinusFastaPath)
                elif confDict['consBkpSeq']:
                    clipped_seqMinusFastaPath = clipped_seqMinus
                metacluster.intLeftBkp = bkpINT(metacluster, clipped_seqMinusFastaPath, blatDbPath, bkpDir, metacluster.identity)


def supportingCLIPPING(metacluster, buffer, confDict, bam, normalBam, mode, outDir, side):
    '''
    # Determine the area where clipping events must be searched (narrow if there are discordant clippings, wide otherwise).
    # Collect clipping events in previous area.
    # Add clippings if there are discodant clippings. Otherwise, return metacluster and clippings list in order to perform BLAT search.

    Input:
        1. metacluster
        2. buffer: buffer to add in wide clippings search
        3. confDict
        4. bam
        5. normalBam
        6. mode
        7. outDir
        8. side: orientation to look for clippings: PLUS or MINUS.
    
    Output:
        1. metacluster
        2. candidate clippings list
        3. discClip: boolean. True if there are discordant clippings, False otherwise.
    '''
    # Note: This function works but you have to allow duplicates in the clipping

    # New dictionary for performing collecting collecting clippings
    clippingConfDict = dict(confDict)
    clippingConfDict['targetSV'] = ['CLIPPING']
    clippingConfDict['minMAPQ'] = 10
    # NOTE SR: Think and check if this is neccessary. I think it is not
    #confDict['minCLIPPINGlen'] = 2
    clippingEventsDict = {}

    ## Define region
    if side == 'PLUS':
        if metacluster.orientation != 'RECIPROCAL':
            # Determine the area where clipping events must be searched (narrow if there are discordant clippings, wide otherwise).
            binBeg, binEnd, discClip = determinePlusBkpArea(metacluster.beg, metacluster.end, metacluster.events, buffer)
        elif metacluster.orientation == 'RECIPROCAL':
            # Collect PLUS discordant events:
            reciprocalPlusEvents = [eventP for eventP in metacluster.events if eventP.orientation == 'PLUS']

            # Determine plus cluster beggining and end
            begPlus = min([eventPlus.beg for eventPlus in reciprocalPlusEvents])
            endPlus = max([eventPlus.end for eventPlus in reciprocalPlusEvents])

            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determinePlusBkpArea(begPlus, endPlus, reciprocalPlusEvents, buffer)

    elif side == 'MINUS':
        if metacluster.orientation != 'RECIPROCAL':
            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determineMinusBkpArea(metacluster.beg, metacluster.end, metacluster.events, buffer)
        elif metacluster.orientation == 'RECIPROCAL':
            # Collect MINUS discordant events
            # As some CLIPPING PLUS events could be added in previous step, now they have to be collected this way
            reciprocalMinusEvents = []
            for eventM in metacluster.events:
                if eventM.type == 'DISCORDANT':
                    if eventM.orientation == 'MINUS':
                        reciprocalMinusEvents.append(eventM)

            # Determine minus cluster beggining and end
            begMinus = min([eventMinus.beg for eventMinus in reciprocalMinusEvents])
            endPlus = max([eventMinus.end for eventMinus in reciprocalMinusEvents])

            # Determine the area where clipping events must be searched (narrow if there are discordatn clippings, wide otherwise).
            binBeg, binEnd, discClip = determineMinusBkpArea(begMinus, endPlus, reciprocalMinusEvents, buffer)

    # Collect clippings
    ref = metacluster.ref

    if mode == "SINGLE":
        clippingEventsDict = bamtools.collectSV(ref, binBeg, binEnd, bam, clippingConfDict, None, False)

    elif mode == "PAIRED":
        clippingEventsDict = bamtools.collectSV_paired(ref, binBeg, binEnd, bam, normalBam, clippingConfDict)

    # When the metacluster is RIGHT and there are some right clippings:
    if side == 'PLUS' and 'RIGHT-CLIPPING' in clippingEventsDict.keys():
        # Keep only those clipping event which have their clipping bkp is in the desired area
        clippingEventsToAdd = chooseBkpClippings(clippingEventsDict, 'RIGHT-CLIPPING', binBeg, binEnd)
        # From dictionary to list
        clippings = list(itertools.chain(*clippingEventsToAdd.values()))
        #metacluster, clippings = addClippings(metacluster, clippingEventsDict, 'RIGHT-CLIPPING', binBeg, binEnd, discClip, confDict['viralDb'])
        return metacluster, clippings, discClip

    # When the metacluster is LEFT and there are some minus clippings:
    elif side == 'MINUS' and 'LEFT-CLIPPING' in clippingEventsDict.keys():
        # Add clippings if there are discodant clippings. 
        # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
        # If there are not discordant clippings and the metacluster has no identity, dont add clippings.
        # Keep only those clipping event which have their clipping bkp is in the desired area
        clippingEventsToAdd = chooseBkpClippings(clippingEventsDict, 'LEFT-CLIPPING', binBeg, binEnd)
        # From dictionary to list
        clippings = list(itertools.chain(*clippingEventsToAdd.values()))
        #metacluster, clippings = addClippings(metacluster, clippingEventsDict, 'LEFT-CLIPPING', binBeg, binEnd, discClip, confDict['viralDb'])
        return metacluster, clippings, discClip
    else:
        return None, None, False


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

#def addDiscClippings(metacluster, clippings, discClip, viralDb):
    '''
    # Add clippings if there are discodant clippings. 
    # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
    # If there are not discordant clippings and the metacluster has no identity, dont add clippings.
    


    # TODO SR: For the moment, no BLAT id performed for MEs insertions.
    # Add clippings if there are discodant clippings.
    if discClip == True or not viralDb or not metacluster.identity:
        metacluster.addEvents(clippings) 
        return None, clippings
    # If there are not discordant clippings, but metacluster has identity, return metacluster and clippings list in order to perform BLAT clippings search
    elif metacluster.identity:
        return metacluster, clippings
    # If there are not discordant clippings and the metacluster has no identity, dont add clippings (return None)
    else:
        return None, None
    '''

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

def addBlatClippings(metaclustersWODiscClip, db, binId, outDir):
    # Write fasta with events collected in a wide region
    clippingEventsToAdd = list(itertools.chain(*metaclustersWODiscClip.values()))
    clippingsFasta = writeClippingsFasta(clippingEventsToAdd, binId, outDir)

    # Align clipped sequences with BLAT against db
    outName = binId + '_clippingsBlat'
    pslPath = alignment.alignment_blat(clippingsFasta, db, outName, outDir)
    pslDict = makePslDict(pslPath)
    matchClippings = collectMatchClippings(metaclustersWODiscClip, pslDict)
    if matchClippings != {}:
        # Collect those clippings that have their bkp in same position as ones in BLAT.
        clippings2Add = collectClipBkpMatch(matchClippings, clippingEventsToAdd)
        if clippings2Add:
            for metacluster, clippings in clippings2Add:
                metacluster.addEvents(clippings)

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
    '''
    Collect those clippings with hits in BLAT search whose hit match with metacluster identities.
    '''
    matchClippings = {}
    for metaclusterWODiscClip, clippings in metaclustersWODiscClip.items():
        for clip in clippings:
            if clip.readName in pslDict.keys():
                # TODO SR: Aqui no s esi habrá que mirar todos los de una lista en otra lista:
                if metaclusterWODiscClip.identity in pslDict[clip.readName]:
                    if metaclusterWODiscClip in matchClippings.keys():
                        matchClippings[metaclusterWODiscClip].append(clip)
                    else:
                        matchClippings[metaclusterWODiscClip] = []
                        matchClippings[metaclusterWODiscClip].append(clip)

    return matchClippings
                
def collectClipBkpMatch(matchClippings, clippingEventsToAdd):
    '''
    '''
    
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


def reconstructSeq(metacluster, consSeq, orientation, outDir):
    '''
    '''
    # Collect discordant clipping events
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
    # If there are discordant clipping events      
    if discClip:
        # Collect clipping objects that are the same alignment as discordant object:
        clippingsDisc = {}
        for eventC in metacluster.events:
            if eventC.type == 'CLIPPING' and eventC.readName in discClip.keys():
                clippingsDisc[eventC] = discClip[eventC.readName]

    
    if clippingsDisc: # This can be empty if clipping of discordant is too small.
        if consSeq: # Make consensus sequence with discordant clipping events
            clipped_seq, consFastaBool = conSeq(metacluster, clippingsDisc, orientation, outDir)
            return clipped_seq, consFastaBool
        elif not consSeq: # Make representative sequence with discordant clipping events
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
            if consSeq: # Make consensus sequence with BLAT clipping events
                clipped_seq, consFastaBool = conSeq(metacluster, clippingsBlat, orientation, outDir)
                return clipped_seq, consFastaBool

            elif not consSeq: # Make representative sequence with BLAT clipping events
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
                if consSeq: # Make consensus sequence with clipping events
                    clipped_seq, consFastaBool = conSeq(metacluster, clippings, orientation, outDir)
                    return clipped_seq, consFastaBool

                elif not consSeq: # Make representative sequence with clipping events
                    clipped_seq = repreSeq(metacluster, orientation, clippings)
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
    if len(clippings) > 1: 
        consFastaBool = True
        if orientation == 'PLUS':
            metacluster.rightSeq, intConsensusPath = newConsensusSeq([*clippings], orientation, outDir)
            return intConsensusPath, consFastaBool
        elif orientation == 'MINUS':
            metacluster.leftSeq, intConsensusPath = newConsensusSeq([*clippings], orientation, outDir)
            return intConsensusPath, consFastaBool
    # If there are only one clipping, this will be the consensus.
    elif len(clippings) == 1:
        consFastaBool = False
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

    return consSeq, intConsensusPath



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

##############################################################################################################################################

def blatClippings_Nouse(clippingsDict, eventType, identity, ID, db, outDir):

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