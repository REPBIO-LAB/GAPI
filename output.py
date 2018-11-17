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
    row = "#ref \t beg \t end \t clusterType \t nbTotal \t nbTumour \t nbNormal \t mean \t std \t cv \n"
    outFile.write(row)

    ## Write clusters ##
    # For each cluster type
    for clusterType in clusters.keys():

        # For each binDb
        for binDb in clusters[clusterType]:
            
            # For each cluster in the binDb
            for cluster in binDb.collect(clusterType):

                # Collect cluster features
                nbTotal, nbTumour, nbNormal = cluster.nbEvents()
                mean, std, cv = cluster.meanLen()

                # Write into output file
                row = "\t".join([cluster.ref, str(cluster.beg), str(cluster.end), clusterType, str(nbTotal), str(nbTumour), str(nbNormal), str(mean), str(std), str(cv), "\n"])
                outFile.write(row)

    ## Close output file ##
    outFile.close()
