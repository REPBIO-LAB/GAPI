import structures


def writeClusters(clusters, outDir):
    '''
    Write structural variation clusters into a tsv file
    '''

    ## Open output file ##
    fileName = "clusters.tsv"
    outFilePath = outDir + '/' + fileName
    outFile = open(outFilePath, 'w')

    ## Write header ##
    row = "#ref \t beg \t end \t binId \t clusterType \t insType \t family \t srcId \t status \t percResolved \t strand \t hits \t nbTotal \t nbTumour \t nbNormal \t nbOutliers \t consensusLen \t meanLen \t std \t cv \t failedFilters \t consensusIns \n"

    outFile.write(row)

    ## Write clusters ##
    # For each cluster type (NOTE: modify at one point since the row ordering is different in each execution)
    for clusterType in clusters.keys():

        # For each binDb
        for binDb in clusters[clusterType]:
            
            binId = binDb.ref + ':' + str(binDb.beg) + '-' + str(binDb.end) 

            # For each cluster in the binDb
            for cluster in binDb.collect(clusterType):

                # Collect cluster features
                nbTotal, nbTumour, nbNormal = cluster.nbEvents()
                mean, std, cv = cluster.meanLen()

                failedFilters = [f for f, v in cluster.filters.items() if not v]
                if len(failedFilters) == 0:
                    failedFilters="."

                # Write into output file
                row = "\t".join([cluster.ref, str(cluster.beg), str(cluster.end), binId, clusterType, str(cluster.insType), str(cluster.family), str(cluster.srcId), str(cluster.status), str(cluster.percResolved), str(cluster.strand), str(cluster.hits), str(nbTotal), str(nbTumour), str(nbNormal), str(cluster.nbOutliers), str(cluster.consensusLen), str(mean), str(std), str(cv), str(",".join(failedFilters)), str(cluster.consensusIns), "\n"])
                outFile.write(row)

         
        

    ## Close output file ##
    outFile.close()
