import structures


def write_INS(INS_metaclusters, outFileName, outDir):
    '''
    Write INS calls into a tsv file

    Input:
        1. INS_metaclusters: list containing list of INS metaclusters
        2. outFileName: Output file name
        3. outDir: Output directory
    Output: tsv file containing identified metaclusters
    '''
    ## 1. Open output file 
    outFilePath = outDir + '/' + outFileName
    outFile = open(outFilePath, 'w')

    ## 2. Write header 
    row = "#ref \t beg \t end \t filters \t mutOrigin \t insType \t family \t subfamily \t cytobandId \t gnName \t biotype \t nbTotal \t nbTumour \t nbNormal \t nbINS \t nbDEL \t nbCLIPPING \t length \t cv \t percResolved \t qHits \t tHits \t insertSeq \n"
    outFile.write(row)

    ## 3. Write INS metaclusters 
    # For each metacluster
    for metacluster in INS_metaclusters:
            
        ## General features 
        filters = 'PASS' if not metacluster.failedFilters else ','.join(metacluster.failedFilters)
        insType = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        nbTotal, nbTumour, nbNormal, nbINS, nbDEL, nbCLIPPING = metacluster.nbEvents()        
        meanLen, cv = metacluster.subclusters['INS'].cv_len() if 'INS' in metacluster.subclusters else (None, None)
        length = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        percResolved = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        qHits = None if metacluster.insertAnnot is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertAnnot.alignments ])
        tHits = None if metacluster.insertAnnot is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertAnnot.alignments ])
        insert = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None

        ## Repeat specific features
        family = metacluster.SV_features['FAMILY'] if 'FAMILY' in metacluster.SV_features else None
        subfamily = metacluster.SV_features['SUBFAMILY'] if 'SUBFAMILY' in metacluster.SV_features else None

        ## Transduction specific features
        cytobandId = metacluster.SV_features['CYTOBAND'] if 'CYTOBAND' in metacluster.SV_features else None

        ## Exon specific features
        gnName = metacluster.SV_features['GENE_NAME'] if 'GENE_NAME' in metacluster.SV_features else None
        biotype = metacluster.SV_features['BIOTYPE'] if 'BIOTYPE' in metacluster.SV_features else None
        
        # Write INS call into output file
        row = "\t".join([metacluster.ref, str(metacluster.beg), str(metacluster.end), str(filters), str(metacluster.mutOrigin), str(insType), str(family), str(subfamily), str(cytobandId), str(gnName), str(biotype), str(nbTotal), str(nbTumour), str(nbNormal), str(nbINS), str(nbDEL), str(nbCLIPPING), str(length), str(cv), str(percResolved), str(qHits), str(tHits), str(insert), "\n"])
        outFile.write(row)

    ## Close output file ##
    outFile.close()

def write_DISCORDANT(discordantClusters, outDir):
    '''
    Write DISCORDANT read pair calls into a tsv file
    '''

    ## 1. Open output file 
    fileName = "DISCORDANT_MEIGA.tsv"
    outFilePath = outDir + '/' + fileName
    outFile = open(outFilePath, 'w')

    ## 2. Write header 
    row = "#ref \t beg \t end \t clusterType \t family \t nbTotal \t nbTumour \t nbNormal \t repeats \t region \t gene \n"
    outFile.write(row)

    ## 3. Write clusters 
    # Iterate over the bins
    for discordantDict in discordantClusters:

        # Iterate over the cluster types
        for key, clusterList in discordantDict.items():
            orientation, clusterType, family = key.split('_')
            clusterType = orientation + '_' + clusterType

            # For each cluster from a given cluster type
            for DISCORDANT in clusterList:

                ## Collect info
                nbTotal, nbTumour, nbNormal = DISCORDANT.nbEvents()
                region, gene = DISCORDANT.geneAnnot if hasattr(DISCORDANT, 'geneAnnot') else ("None", "None")
            
                ## TEMPORARY: Only report discordant clusters that:
                # - Mate do not aligns on a retrotransposon
                if (family != 'None'):

                    # Write DISCORDANT cluster into output file
                    row = "\t".join([DISCORDANT.ref, str(DISCORDANT.beg), str(DISCORDANT.end), clusterType, family, str(nbTotal), str(nbTumour), str(nbNormal), str(DISCORDANT.repeatAnnot), region, gene, "\n"])
                    outFile.write(row)

def writeMetaclusters(metaclustersList, outDir):
    '''
    Write structural variation clusters into a tsv file
    '''

    ## Open output file ##
    fileName = "metaclusters.tsv"
    outFilePath = outDir + '/' + fileName
    outFile = open(outFilePath, 'w')

    for dictMetacluster in metaclustersList:
        if dictMetacluster != None:
            for metacluster,d2 in dictMetacluster.items():
                row = 'METACLUSTER: ' + str(metacluster) +' '+ str(len(metacluster.events)) +' '+ str(metacluster.ref) +' '+ str(metacluster.beg) +' '+ str(metacluster.end) +' '+ str(metacluster.intOrigin) + '\n'
                outFile.write(row)
                for event in metacluster.events:
                        if event.type == 'DISCORDANT':
                            row = str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' ' + str(event.identity) + ' ' + str(event.side) + ' ' + str(event.sample) + '\n'
                            outFile.write (row)
                        else:
                            row = str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' None ' + str(event.clippedSide) + ' ' + str(event.sample) + '\n'
                            outFile.write (row)
                for k,v in d2.items():
                    row = str(k) + ' = ' + str(v) + '\n'
                    outFile.write (row)  
         
    ## Close output file ##
    outFile.close()
