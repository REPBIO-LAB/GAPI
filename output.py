import callers
import structures


def writeOutput(clusters, outDir):
    '''
    Write output file
    '''

    # Open output file
    outfile = "architectTempOutput.txt"
    outfilepath = outDir + '/' + outfile
    output = open(outfilepath, 'w')

    # Write header
    output.write("ref \t beg \t end \t clusterType \t nbTotalEvents \t nbTumourEvents \t nbNormalEvents \t mean \t standardDeviation \t variationConfficient \n")

    # Write core
    writeEventsFeatures(clusters, output)

    # Close output file
    output.close()

def writeEventsFeatures(clusters, output):
    '''
    Write output core
    '''

    # List containing all event types
    eventTypes = ["INS-CLUSTER", "DEL-CLUSTER", "LEFT-CLIPPING-CLUSTER", "RIGHT-CLIPPING-CLUSTER"]

    # For each list in cluster, pick each event and write its features
    for i in range (0,len(eventTypes)):

	    for TYPE_cluster in clusters[i]:
		    events = TYPE_cluster.collect(eventTypes[i])

		    for event in events:
			    nbTotalEvents, nbTumourEvents, nbNormalEvents = event.nbEvents()
			    mean, std, cv = event.meanLen()
			    output.write(str(event.ref) + "\t" + str(event.beg) + "\t" + str(event.end) + "\t" + str(eventTypes[i])  + "\t"  + str(nbTotalEvents) + "\t" + str(nbTumourEvents) + "\t" + str(nbNormalEvents) + "\t" + str(mean) + "\t" + str(std)  + "\t" + str(cv)+ "\n")
