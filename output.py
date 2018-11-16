import callers
import structures


def writeOutput(mode, bam, normalBam, confDict, outDir):
    '''
    FASTA file writer. Write data stored in the dictionary into a FASTA file
    '''

    clusters = []

    callerObj = callers.SVcaller_nano(mode, bam, normalBam, confDict, outDir)
    INS_clusters, DEL_clusters, left_CLIPPING_clusters, right_CLIPPING_clusters = callerObj.callSV()


    outfile = "architectTempOutput.txt"
    outfilepath = outDir + '/' + outfile
    output = open(outfilepath, 'w')

    output.write("ref \t beg \t end \t clusterType \t nbTotalEvents \t nbTumourEvents \t nbNormalEvents \t medianLength \n")

    clusters = [INS_clusters, DEL_clusters, left_CLIPPING_clusters, right_CLIPPING_clusters]
    eventTypes = ["INS-CLUSTER", "DEL-CLUSTER", "LEFT-CLIPPING-CLUSTER", "RIGHT-CLIPPING-CLUSTER"]

    for i in range (0,len(eventTypes)):
        for TYPE_cluster in clusters[i]:
            events = TYPE_cluster.collect(eventTypes[i])
            for event in events:
                print (event)
                nbTotalEvents, nbTumourEvents, nbNormalEvents = event.nbEvents()
                output.write(str(event.ref) + "\t" + str(event.beg) + "\t" + str(event.end) + "\t"  + eventTypes[i]  + "\t"  + str(nbTotalEvents) + "\t" + str(nbTumourEvents) + "\t" + str(nbNormalEvents) + "\n")

    output.close()
