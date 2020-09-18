'''
Module 'bamtools VIGA SR' - Contains functions for extracting data from bam files specific for virus short reads
'''

## DEPENDENCIES ##
# External
import pysam
from cigar import Cigar
import numpy as np
import os
import Bio.SeqUtils
from Bio.SeqUtils import lcc

# Internal
import events
import formats
import sequences

def collectDiscodantsLowMAPQSeq(ref, binBeg, binEnd, bam, discordantMatesMaxMAPQ, discordantMatesCheckUnmapped, discordantMatesSupplementary, discordantMatesMaxBasePerc, discordantMatesMinLcc, outDir):
    '''
    Collecting read name and sequence of discordant low quality reads from all bam refs

    Input:
        1. ref: target referenge
        2. binBeg: bin begin
        3. binEnd: bin end
        4. bam: indexed BAM file
        5. discordantMatesMaxMAPQ: Maximum mapping quality used for collecting dicordant read mates.
        6. discordantMatesCheckUnmapped: Boolean. If True, when a dicordant read mate is unmapped, collect it no matter its MAPQ
        7. discordantMatesSupplementary: Boolean. When False, avoid collecting dicordant read mates that are supplementary alignments.
        8. discordantMatesMaxBasePerc: Maximum base percentage of discordant read mates sequences.
        9. discordantMatesMinLcc: Minimum local complexity of discordant read mates sequences.
        10. outDir: Output directory
    
    Output:
        1. Doesn't return anything. It creates a fasta file with discordant low quality reads from all bam refs.
    '''
    # TODO SR: Think if filterDuplicates step is neccesary in collectSeq method and implement it if so.
    #filterDuplicates = True

    ## Initialize dictionary to store SV events
    eventsSeqDict = {}

    ## Open BAM file for reading
    bamFile = pysam.AlignmentFile(bam, "rb")

    ## Extract alignments
    iterator = bamFile.fetch(ref, binBeg, binEnd)
    
    # For each read alignment
    for alignmentObj in iterator:

        ### Filter out alignments based on different criteria:
        MAPQ = int(alignmentObj.mapping_quality) # Mapping quality

        ## No query sequence available
        if alignmentObj.query_sequence == None:
            continue

        ## Aligments with MAPQ < threshold
        if (MAPQ > discordantMatesMaxMAPQ):
            continue

        # TODO SR: Think if filterDuplicates step is neccesary in collectSeq method and implement it if so.
        ## Duplicates filtering enabled and duplicate alignment
        #if (confDict['filterDuplicates'] == True) and (alignmentObj.is_duplicate == True):
            #continue

        # Filter supplementary alignments if TRUE. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)
        if discordantMatesSupplementary == False and alignmentObj.is_supplementary == True:
            continue
        
        ## Collect DISCORDANT low quality reads in dictionary -> eventsSeqDict[readName] = readSequence

        if not alignmentObj.is_proper_pair:

            # Pick sequences that are unmmapped or with mapping quality < discordantMatesMaxMAPQ
            if discordantMatesCheckUnmapped == True:
                # TODO: Put decoy as an option
                if (alignmentObj.is_unmapped == True) or (MAPQ < discordantMatesMaxMAPQ) or (alignmentObj.reference_name == 'NC_007605') or (alignmentObj.reference_name == 'EBV'):
                    # Calculate base percentage
                    basePercs = sequences.baseComposition(alignmentObj.query_sequence)[1]
                    # Delete total value of base percentage result
                    del basePercs['total']
                    # Only those sequences with base percentage lower than 85 are collected:
                    if all(perc < discordantMatesMaxBasePerc for perc in basePercs.values()):
                        # Check local complexity of sequences:
                        complexity = Bio.SeqUtils.lcc.lcc_simp(alignmentObj.query_sequence)
                        # Only those sequences with local complexity lower than 1.49 are collected:
                        if complexity > discordantMatesMinLcc:
                            eventsSeqDict[alignmentObj.query_name]=alignmentObj.query_sequence

            # Pick sequences with mapping quality < discordantMatesMaxMAPQ
            else:
                # TODO; Put decoy as an option
                if (MAPQ < discordantMatesMaxMAPQ) or (alignmentObj.reference_name == 'NC_007605') or (alignmentObj.reference_name == 'EBV'):
                    # Calculate base percentage
                    basePercs = sequences.baseComposition(alignmentObj.query_sequence)[1]
                    # Delete total value of base percentage result
                    del basePercs['total']
                    # Only those sequences with base percentage lower than 85 are collected:
                    if all(perc < discordantMatesMaxBasePerc for perc in basePercs.values()):
                        # Check local complexity of sequences:
                        complexity = Bio.SeqUtils.lcc.lcc_simp(alignmentObj.query_sequence)
                        # Only those sequences with local complexity lower than 1.49 are collected:
                        if complexity > discordantMatesMinLcc:
                            eventsSeqDict[alignmentObj.query_name]=alignmentObj.query_sequence
        
    ## Close 
    bamFile.close()
    collectVirusDir = outDir + '/COLLECT_VIRUS'
    # Set output FASTA file name
    allFastas_all = collectVirusDir + "/allFastas_all.fasta"
    # Create FASTA object
    seqsFastaObj= formats.FASTA()
    # Create FASTA dictionary
    seqsFastaObj.seqDict = eventsSeqDict
    # Write output FASTA
    seqsFastaObj.write(allFastas_all, 'append', True)

    # Delete useless variables
    del eventsSeqDict
    del seqsFastaObj

    return

def filterBAM2FastaDict(BAM, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='SR'):
    '''
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}

    # For each read alignment
    for alignmentObj in iterator:
        alignmentPass = False
        #numMatches = 0
        queryCoord = 0
        print ('alignmentObj.query_name ' + str(alignmentObj.query_name))
        #print ('alignmentObj.is_unmapped ' + str(alignmentObj.is_unmapped))

        if not alignmentObj.is_unmapped:
            ctuples = alignmentObj.cigartuples
            allMatches = [t[1] for t in ctuples if t[0] == 0]
            totalMatch = sum (allMatches)
            print ('totalMatch ' + str(totalMatch))
            if totalMatch >= minTotalMatchVirus:
                c = Cigar(alignmentObj.cigarstring)
                for citem  in list(c.items()):
                    # If cigar is query consuming, update query coordinates:
                    if citem[1] != 'M' and citem[1] != 'D' and citem[1] != '=':
                        queryCoord = queryCoord + int(citem[0])
                    elif citem[1] == 'M' or citem[1] == '=':
                        print ('citem[0] ' + str(citem[0]))
                        print ('maxMatchCheckMAPQVirus'  + str(maxMatchCheckMAPQVirus))
                        print ('alignmentObj.mapping_quality'  + str(alignmentObj.mapping_quality))
                        print ('minMAPQVirus'  + str(minMAPQVirus))
                        if citem[0] >= minParcialMatchVirus:
                            if (citem[0] <= maxMatchCheckMAPQVirus and alignmentObj.mapping_quality > minMAPQVirus) or citem[0] > maxMatchCheckMAPQVirus:
                                sequence = alignmentObj.query_sequence[queryCoord:(queryCoord + int(citem[0]))]
                                # Calculate base percentage
                                basePercs = sequences.baseComposition(sequence)[1]
                                # Delete total value of base percentage result
                                del basePercs['total']
                                print ('basePercs ' + str(basePercs))
                                # Only those sequences with base percentage lower than 85 are collected:
                                if all(perc < maxBasePercVirus for perc in basePercs.values()):
                                    # Calculate complexity
                                    complexity = Bio.SeqUtils.lcc.lcc_simp(sequence)
                                    print ('complexity ' + str(complexity))
                                    print ('minLccVirus ' + str(minLccVirus))
                                    if complexity > minLccVirus:
                                        #print ('PassAll')
                                        # TODO: put this in a good way!!!!!
                                        if mode=='LR':
                                            if len(sequence)>10:
                                                alignmentPass = True
                                                break # Cambio 28/08 11:53
                                        else:
                                            alignmentPass = True
                                            break # Cambio 28/08 11:53
                                    else:
                                        queryCoord = queryCoord + int(citem[0])
                                else:
                                    queryCoord = queryCoord + int(citem[0])
                        else:
                            queryCoord = queryCoord + int(citem[0])
        '''
        print (numMatches)
        if numMatches > 90: # TODO SR: put this as an option
            alignmentPass = True
        elif alignmentObj.mapping_quality >= viralBamMAPQ and numMatches >= int(selectedPartialMatches):
            alignmentPass = True
        # New condition:
        # TODO SR: Put as options!!!
        elif alignmentObj.mapping_quality >= 40 and numMatches >= 60:
            alignmentPass = True
        '''
        print ('alignmentPass ' + str(alignmentPass))
        if alignmentPass == True:
            # Add to fasta dict
            if alignmentObj.query_name in fastaDict.keys():
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
            else:
                fastaDict[alignmentObj.query_name] = []
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)

    return fastaDict

def filterBAM2FastaDict_LRGood(BAM, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='SR'):
    '''
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}
    totalMatch = None

    # For each read alignment
    for alignmentObj in iterator:
        alignmentPass = False
        #numMatches = 0
        queryCoord = 0
        print ('alignmentObj.query_name ' + str(alignmentObj.query_name))
        #print ('alignmentObj.is_unmapped ' + str(alignmentObj.is_unmapped))
        citems = []

        if not alignmentObj.is_unmapped:

            print ('HOLAalignmentObj.query_sequence ' + str(alignmentObj.query_sequence))
            TOTALcomplexity = Bio.SeqUtils.lcc.lcc_simp(alignmentObj.query_sequence)
            print ('TOTALcomplexity ' + str(TOTALcomplexity))
            # TODO put option!!!:
            if TOTALcomplexity > minLccVirus:

                ctuples = alignmentObj.cigartuples
                allMatches = [t[1] for t in ctuples if t[0] == 0]
                totalMatch = sum (allMatches)
                print ('totalMatch ' + str(totalMatch))
                if totalMatch >= minTotalMatchVirus:
                    c = Cigar(alignmentObj.cigarstring)
                    for citem  in list(c.items()):
                        # If cigar is query consuming, update query coordinates:
                        if citem[1] != 'M' and citem[1] != 'D' and citem[1] != '=':
                            queryCoord = queryCoord + int(citem[0])
                        elif citem[1] == 'M' or citem[1] == '=':
                            print ('citem[0] ' + str(citem[0]))
                            print ('maxMatchCheckMAPQVirus'  + str(maxMatchCheckMAPQVirus))
                            print ('alignmentObj.mapping_quality'  + str(alignmentObj.mapping_quality))
                            print ('minMAPQVirus'  + str(minMAPQVirus))
                            minMAPQVirus2 = 1
                            if citem[0] >= minParcialMatchVirus:
                                if (citem[0] <= maxMatchCheckMAPQVirus and alignmentObj.mapping_quality > minMAPQVirus) or citem[0] > maxMatchCheckMAPQVirus:
                                    #if ((citem[0]/alignmentObj.query_length) <= .05 and alignmentObj.mapping_quality > minMAPQVirus2) or (citem[0]/alignmentObj.query_length) > .05:
                                    if ((totalMatch/alignmentObj.query_length) <= .035 and alignmentObj.mapping_quality > minMAPQVirus2) or (totalMatch/alignmentObj.query_length) > .035: # cambio 28/08 12:57, 15:34, 15:47


                                        sequence = alignmentObj.query_sequence[queryCoord:(queryCoord + int(citem[0]))]

                                        #if len(sequence)/totalMatch >= 0.015:
                                        print ('len(sequence)/totalMatch ' + str(len(sequence)/totalMatch))
                                    
                                        # Calculate base percentage
                                        basePercs = sequences.baseComposition(sequence)[1]
                                        # Delete total value of base percentage result
                                        del basePercs['total']
                                        print ('basePercs ' + str(basePercs))
                                        # Only those sequences with base percentage lower than 85 are collected:
                                        if all(perc < maxBasePercVirus for perc in basePercs.values()):
                                            # Calculate complexity
                                            complexity = Bio.SeqUtils.lcc.lcc_simp(sequence)
                                            print ('complexity ' + str(complexity))
                                            print ('minLccVirus ' + str(minLccVirus))
                                            if complexity > minLccVirus:
                                                #print ('PassAll')
                                                # TODO: put this in a good way!!!!!
                                                if mode=='LR':
                                                    if len(sequence)>10:
                                                        alignmentPass = True
                                                        citems.append(citem[0])
                                                        if (sum(citems)/totalMatch) > 0.14:
                                                            break # Cambio 28/08 11:53
                                                else:
                                                    alignmentPass = True
                                                    citems.append(citem[0])
                                                    if (sum(citems)/totalMatch) > 0.14:
                                                        break # Cambio 28/08 11:53
                                            else:
                                                queryCoord = queryCoord + int(citem[0])
                                        else:
                                            queryCoord = queryCoord + int(citem[0])
                                    else:
                                        queryCoord = queryCoord + int(citem[0])
                                else:
                                    queryCoord = queryCoord + int(citem[0])
                            else:
                                queryCoord = queryCoord + int(citem[0])
        '''
        print (numMatches)
        if numMatches > 90: # TODO SR: put this as an option
            alignmentPass = True
        elif alignmentObj.mapping_quality >= viralBamMAPQ and numMatches >= int(selectedPartialMatches):
            alignmentPass = True
        # New condition:
        # TODO SR: Put as options!!!
        elif alignmentObj.mapping_quality >= 40 and numMatches >= 60:
            alignmentPass = True
        '''
        print ('citems ' + str(citems))
        print ('sum(citems) ' + str(sum(citems)))
        print ('alignmentPass ' + str(alignmentPass))
        if totalMatch:
            print ('(sum(citems)/totalMatch) ' + str((sum(citems)/totalMatch)))
        if totalMatch:
            print ('totalMatch ' + str(totalMatch))
            if (sum(citems)/totalMatch) > 0.14:
                # Add to fasta dict
                if alignmentObj.query_name in fastaDict.keys():
                    fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
                else:
                    fastaDict[alignmentObj.query_name] = []
                    fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)

    return fastaDict

def filterBAM2FastaDict_minimap(BAM, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='SR'):
    '''
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}
    totalMatch = None

    # For each read alignment
    for alignmentObj in iterator:
        alignmentPass = False
        #numMatches = 0
        queryCoord = 0
        print ('alignmentObj.query_name ' + str(alignmentObj.query_name))
        #print ('alignmentObj.is_unmapped ' + str(alignmentObj.is_unmapped))
        citems = []

        if not alignmentObj.is_unmapped and alignmentObj.query_sequence:
            print ('HOLAalignmentObj.query_sequence ' + str(alignmentObj.query_sequence))
            TOTALcomplexity = Bio.SeqUtils.lcc.lcc_simp(alignmentObj.query_sequence)
            print ('TOTALcomplexity ' + str(TOTALcomplexity))
            # TODO put option!!!:
            if TOTALcomplexity > minLccVirus:
                ctuples = alignmentObj.cigartuples
                allMatches = [t[1] for t in ctuples if t[0] == 0]
                totalMatch = sum (allMatches)
                print ('totalMatch ' + str(totalMatch))
                if totalMatch >= minTotalMatchVirus:
                    c = Cigar(alignmentObj.cigarstring)
                    for citem  in list(c.items()):
                        # If cigar is query consuming, update query coordinates:
                        print ('citam ' + str(citem))
                        if citem[1] != 'M' and citem[1] != 'D' and citem[1] != '=' and citem[1] != 'H': # CAMBIO 27/08
                            print ('citem[1] ' + str(citem[1]))
                            print ('queryCoord0 ' + str(queryCoord))
                            queryCoord = queryCoord + int(citem[0])
                            print ('queryCoord1 ' + str(queryCoord))
                        elif citem[1] == 'M' or citem[1] == '=':
                            print ('citem[0] ' + str(citem[0]))
                            print ('maxMatchCheckMAPQVirus'  + str(maxMatchCheckMAPQVirus))
                            print ('alignmentObj.mapping_quality'  + str(alignmentObj.mapping_quality))
                            print ('minMAPQVirus'  + str(minMAPQVirus))
                            print ('queryCoord ' + str(queryCoord))
                            print ('alignmentObj.query_length ' + str(alignmentObj.query_length))
                            minMAPQVirus2 = 1
                            if citem[0] >= minParcialMatchVirus:
                                #if citem[0] > maxMatchCheckMAPQVirus:
                                if (citem[0] <= maxMatchCheckMAPQVirus and alignmentObj.mapping_quality > minMAPQVirus) or citem[0] >= maxMatchCheckMAPQVirus: # cambio 28/08 15:32
                                    #if ((totalMatch/alignmentObj.query_length) <= .30 and alignmentObj.mapping_quality > minMAPQVirus) or (totalMatch/alignmentObj.query_length) > .30:
                                    if ((totalMatch/alignmentObj.query_length) <= .035 and alignmentObj.mapping_quality > minMAPQVirus2) or (totalMatch/alignmentObj.query_length) > .035: # cambio 28/08 12:57, 15:34, 15:47
                                        sequence = alignmentObj.query_sequence[queryCoord:(queryCoord + int(citem[0]))]
                                        print ('sequence ' + str(sequence))

                                        #if len(sequence)/totalMatch >= 0.022:
                                        print ('len(sequence)/totalMatch ' + str(len(sequence)/totalMatch))

                                        # Calculate base percentage
                                        basePercs = sequences.baseComposition(sequence)[1]
                                        # Delete total value of base percentage result
                                        del basePercs['total']
                                        print ('basePercs ' + str(basePercs))
                                        # Only those sequences with base percentage lower than 85 are collected:
                                        if all(perc < maxBasePercVirus for perc in basePercs.values()):
                                            # Calculate complexity
                                            complexity = Bio.SeqUtils.lcc.lcc_simp(sequence)
                                            print ('complexity ' + str(complexity))
                                            print ('minLccVirus ' + str(minLccVirus))
                                            if complexity > minLccVirus:
                                                #print ('PassAll')
                                                # TODO: put this in a good way!!!!!
                                                if mode=='LR':
                                                    if len(sequence)>10:
                                                        alignmentPass = True
                                                        citems.append(citem[0])
                                                        if (sum(citems)/totalMatch) > 0.14:
                                                            break # Cambio 28/08 11:53
                                                else:
                                                    alignmentPass = True
                                                    citems.append(citem[0])
                                                    if (sum(citems)/totalMatch) > 0.14:
                                                        break # Cambio 28/08 11:53
                                            else:
                                                queryCoord = queryCoord + int(citem[0])
                                        else:
                                            queryCoord = queryCoord + int(citem[0])                   
                                    else:
                                        queryCoord = queryCoord + int(citem[0])
                                else:
                                    queryCoord = queryCoord + int(citem[0])
                            else:
                                queryCoord = queryCoord + int(citem[0])
        '''
        print (numMatches)
        if numMatches > 90: # TODO SR: put this as an option
            alignmentPass = True
        elif alignmentObj.mapping_quality >= viralBamMAPQ and numMatches >= int(selectedPartialMatches):
            alignmentPass = True
        # New condition:
        # TODO SR: Put as options!!!
        elif alignmentObj.mapping_quality >= 40 and numMatches >= 60:
            alignmentPass = True
        '''
        print ('citems ' + str(citems))
        print ('sum(citems) ' + str(sum(citems)))
        print ('alignmentPass ' + str(alignmentPass))
        if totalMatch:
            print ('(sum(citems)/totalMatch) ' + str((sum(citems)/totalMatch)))
        if totalMatch:
            print ('totalMatch ' + str(totalMatch))
            if (sum(citems)/totalMatch) > 0.14:
                # Add to fasta dict
                if alignmentObj.query_name in fastaDict.keys():
                    fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
                else:
                    fastaDict[alignmentObj.query_name] = []
                    fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)

    return fastaDict

def filterBAM2FastaDict_LR_deprecated(BAM, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='SR'):
    '''
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}

    # For each read alignment
    for alignmentObj in iterator:
        alignmentPass = False
        #numMatches = 0
        queryCoord = 0
        print ('alignmentObj.query_name ' + str(alignmentObj.query_name))
        print ('query_length ' + str(alignmentObj.query_length))
        #print ('alignmentObj.is_unmapped ' + str(alignmentObj.is_unmapped))

        if not alignmentObj.is_unmapped:
            ctuples = alignmentObj.cigartuples
            allMatches = [t[1] for t in ctuples if t[0] == 0]
            totalMatch = sum (allMatches)
            print ('totalMatch ' + str(totalMatch))
            if totalMatch >= minTotalMatchVirus:
                c = Cigar(alignmentObj.cigarstring)
                for citem  in list(c.items()):
                    # If cigar is query consuming, update query coordinates:
                    if citem[1] != 'M' and citem[1] != 'D' and citem[1] != '=':
                        queryCoord = queryCoord + int(citem[0])
                    elif citem[1] == 'M' or citem[1] == '=':
                        print ('citem[0] ' + str(citem[0]))
                        print ('maxMatchCheckMAPQVirus'  + str(maxMatchCheckMAPQVirus))
                        print ('alignmentObj.mapping_quality'  + str(alignmentObj.mapping_quality))
                        print ('minMAPQVirus'  + str(minMAPQVirus))
                        if citem[0] >= minParcialMatchVirus:
                            if (citem[0] <= maxMatchCheckMAPQVirus and alignmentObj.mapping_quality > minMAPQVirus) or citem[0] > maxMatchCheckMAPQVirus:
                            #if ((totalMatch/alignmentObj.query_length) <= .30 and alignmentObj.mapping_quality > minMAPQVirus) or (totalMatch/alignmentObj.query_length) > .30:
                                sequence = alignmentObj.query_sequence[queryCoord:(queryCoord + int(citem[0]))]
                                # Calculate base percentage
                                basePercs = sequences.baseComposition(sequence)[1]
                                # Delete total value of base percentage result
                                del basePercs['total']
                                print ('basePercs ' + str(basePercs))
                                # Only those sequences with base percentage lower than 85 are collected:
                                if all(perc < maxBasePercVirus for perc in basePercs.values()):
                                    # Calculate complexity
                                    complexity = Bio.SeqUtils.lcc.lcc_simp(sequence)
                                    print ('complexity ' + str(complexity))
                                    print ('minLccVirus ' + str(minLccVirus))
                                    if complexity > minLccVirus:
                                        #print ('PassAll')
                                        # TODO: put this in a good way!!!!!
                                        if mode=='LR':
                                            if len(sequence)>10:
                                                alignmentPass = True
                                        else:
                                            alignmentPass = True
                                        break
                                    else:
                                        queryCoord = queryCoord + int(citem[0])
                                else:
                                    queryCoord = queryCoord + int(citem[0])
                        else:
                            queryCoord = queryCoord + int(citem[0])
        '''
        print (numMatches)
        if numMatches > 90: # TODO SR: put this as an option
            alignmentPass = True
        elif alignmentObj.mapping_quality >= viralBamMAPQ and numMatches >= int(selectedPartialMatches):
            alignmentPass = True
        # New condition:
        # TODO SR: Put as options!!!
        elif alignmentObj.mapping_quality >= 40 and numMatches >= 60:
            alignmentPass = True
        '''
        print ('alignmentPass ' + str(alignmentPass))
        if alignmentPass == True:
            # Add to fasta dict
            if alignmentObj.query_name in fastaDict.keys():
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
            else:
                fastaDict[alignmentObj.query_name] = []
                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)

    return fastaDict


def filterBAM2FastaDict_LR(BAM, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='SR', allHits_viral=None):
    '''
    Returns dictionary with [alignment query name]=alignment identity for those alignments that PASS the filters. 
    '''

    # Read bam and store in a dictionary
    bamFile = pysam.AlignmentFile(BAM, 'rb')

    iterator = bamFile.fetch()
    
    fastaDict= {}
    totalMatch = None

    # For each read alignment
    for alignmentObj in iterator:
        alignmentPass = False
        #numMatches = 0
        queryCoord = 0
        print ('alignmentObj.query_name ' + str(alignmentObj.query_name))
        #print ('alignmentObj.is_unmapped ' + str(alignmentObj.is_unmapped))
        citems = []

        if not alignmentObj.is_unmapped and alignmentObj.query_sequence:
            print ('HOLAalignmentObj.query_sequence ' + str(alignmentObj.query_sequence))
            TOTALcomplexity = Bio.SeqUtils.lcc.lcc_simp(alignmentObj.query_sequence)
            print ('TOTALcomplexity ' + str(TOTALcomplexity))
            # TODO put option!!!:
            if TOTALcomplexity > minLccVirus:
                ctuples = alignmentObj.cigartuples
                allMatches = [t[1] for t in ctuples if t[0] == 0]
                totalMatch = sum (allMatches)
                print ('totalMatch ' + str(totalMatch))
                if totalMatch >= minTotalMatchVirus:
                    c = Cigar(alignmentObj.cigarstring)
                    for citem  in list(c.items()):
                        # If cigar is query consuming, update query coordinates:
                        print ('citam ' + str(citem))
                        if citem[1] != 'M' and citem[1] != 'D' and citem[1] != '=' and citem[1] != 'H': # CAMBIO 27/08
                            print ('citem[1] ' + str(citem[1]))
                            print ('queryCoord0 ' + str(queryCoord))
                            queryCoord = queryCoord + int(citem[0])
                            print ('queryCoord1 ' + str(queryCoord))
                        elif citem[1] == 'M' or citem[1] == '=':
                            print ('citem[0] ' + str(citem[0]))
                            print ('maxMatchCheckMAPQVirus'  + str(maxMatchCheckMAPQVirus))
                            print ('alignmentObj.mapping_quality'  + str(alignmentObj.mapping_quality))
                            print ('minMAPQVirus'  + str(minMAPQVirus))
                            print ('queryCoord ' + str(queryCoord))
                            print ('alignmentObj.query_length ' + str(alignmentObj.query_length))
                            minMAPQVirus2 = 1
                            if citem[0] >= minParcialMatchVirus:
                                #if citem[0] > maxMatchCheckMAPQVirus:
                                if (citem[0] <= maxMatchCheckMAPQVirus and alignmentObj.mapping_quality > minMAPQVirus) or citem[0] >= maxMatchCheckMAPQVirus: # cambio 28/08 15:32
                                    #if ((totalMatch/alignmentObj.query_length) <= .30 and alignmentObj.mapping_quality > minMAPQVirus) or (totalMatch/alignmentObj.query_length) > .30:
                                    if ((totalMatch/alignmentObj.query_length) <= .035 and alignmentObj.mapping_quality > minMAPQVirus2) or (totalMatch/alignmentObj.query_length) > .035: # cambio 28/08 12:57, 15:34, 15:47
                                        sequence = alignmentObj.query_sequence[queryCoord:(queryCoord + int(citem[0]))]
                                        print ('sequence ' + str(sequence))

                                        #if len(sequence)/totalMatch >= 0.022:
                                        print ('len(sequence)/totalMatch ' + str(len(sequence)/totalMatch))

                                        # Calculate base percentage
                                        basePercs = sequences.baseComposition(sequence)[1]
                                        # Delete total value of base percentage result
                                        del basePercs['total']
                                        print ('basePercs ' + str(basePercs))
                                        # Only those sequences with base percentage lower than 85 are collected:
                                        if all(perc < maxBasePercVirus for perc in basePercs.values()):
                                            # Calculate complexity
                                            complexity = Bio.SeqUtils.lcc.lcc_simp(sequence)
                                            print ('complexity ' + str(complexity))
                                            print ('minLccVirus ' + str(minLccVirus))
                                            if complexity > minLccVirus:
                                                #print ('PassAll')
                                                # TODO: put this in a good way!!!!!
                                                if mode=='LR':
                                                    if len(sequence)>10:
                                                        alignmentPass = True
                                                        citems.append(citem[0])
                                                        if (sum(citems)/totalMatch) > 0.14:
                                                            break # Cambio 28/08 11:53
                                                else:
                                                    alignmentPass = True
                                                    citems.append(citem[0])
                                                    if (sum(citems)/totalMatch) > 0.14:
                                                        break # Cambio 28/08 11:53
                                            else:
                                                queryCoord = queryCoord + int(citem[0])
                                        else:
                                            queryCoord = queryCoord + int(citem[0])                   
                                    else:
                                        queryCoord = queryCoord + int(citem[0])
                                else:
                                    queryCoord = queryCoord + int(citem[0])
                            else:
                                queryCoord = queryCoord + int(citem[0])
        '''
        print (numMatches)
        if numMatches > 90: # TODO SR: put this as an option
            alignmentPass = True
        elif alignmentObj.mapping_quality >= viralBamMAPQ and numMatches >= int(selectedPartialMatches):
            alignmentPass = True
        # New condition:
        # TODO SR: Put as options!!!
        elif alignmentObj.mapping_quality >= 40 and numMatches >= 60:
            alignmentPass = True
        '''
        print ('citems ' + str(citems))
        print ('sum(citems) ' + str(sum(citems)))
        print ('alignmentPass ' + str(alignmentPass))
        if totalMatch:
            print ('(sum(citems)/totalMatch) ' + str((sum(citems)/totalMatch)))
        if totalMatch:
            print ('totalMatch ' + str(totalMatch))
            if allHits_viral != None and (sum(citems)/totalMatch) > 0.14 and (sum(citems)/totalMatch) < 0.144:
                # si va muy justo check la mapq de minimap
                #if (sum(citems)/totalMatch) > 0.14 and (sum(citems)/totalMatch) < 0.142:
                print ('ENTRA')
                if alignmentObj.query_name in allHits_viral.keys():
                    for alig in allHits_viral[alignmentObj.query_name].alignments:
                        if alig.MAPQ > 1 or alig.nbMatches > 1000:
                            # Add to fasta dict
                            if alignmentObj.query_name in fastaDict.keys():
                                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
                            else:
                                fastaDict[alignmentObj.query_name] = []
                                fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
            else:  
                if (sum(citems)/totalMatch) > 0.14:
                    print ('NO ENTRA')
                    # Add to fasta dict
                    if alignmentObj.query_name in fastaDict.keys():
                        fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)
                    else:
                        fastaDict[alignmentObj.query_name] = []
                        fastaDict[alignmentObj.query_name].append(alignmentObj.reference_name)

    return fastaDict