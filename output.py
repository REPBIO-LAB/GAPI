import structures


def write_INS(metaclustersBinDb, outDir):
    '''
    Write INS calls into a tsv file
    '''

    ## 1. Open output file 
    fileName = "INS_MEIGA.tsv"
    outFilePath = outDir + '/' + fileName
    outFile = open(outFilePath, 'w')

    ## 2. Write header 
    row = "#ref \t beg \t end \t binId \t insType \t family \t srcId \t status \t percResolved \t strand \t hits \t nbTotal \t nbTumour \t nbNormal \t nbINS \t nbDEL \t nbCLIPPING \t length \t insertSeq \n"
    outFile.write(row)

    ## 3. Write clusters 
    # For each binDb
    for binDb in metaclustersBinDb:
            
        binId = binDb.ref + ':' + str(binDb.beg) + '-' + str(binDb.end) 

        # For each INS in the binDb
        for INS in binDb.collect(['INS']):

            # Collect INS features
            nbTotal, nbTumour, nbNormal, nbINS, nbDEL, nbCLIPPING = INS.nbEvents()

            insType = INS.INS_features['insType'] if 'insType' in INS.INS_features else None
            family = INS.INS_features['family'] if 'family' in INS.INS_features else None
            srcId = INS.INS_features['srcId'] if 'srcId' in INS.INS_features else None
            status = INS.INS_features['status'] if 'status' in INS.INS_features else None
            percResolved = INS.INS_features['percResolved'] if 'percResolved' in INS.INS_features else None
            strand = INS.INS_features['strand'] if 'strand' in INS.INS_features else None
            hits = INS.INS_features['hits'] if 'hits' in INS.INS_features else None

            # Write INS call into output file
            row = "\t".join([INS.ref, str(INS.beg), str(INS.end), binId, str(insType), str(family), str(srcId), str(status), str(percResolved), str(strand), str(hits), str(nbTotal), str(nbTumour), str(nbNormal), str(nbINS), str(nbDEL), str(nbCLIPPING), str(INS.consensus.length), INS.consensus.pick_insert(), "\n"])
            outFile.write(row)

    ## Close output file ##
    outFile.close()
