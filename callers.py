'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import multiprocessing as mp

# Internal
import log
import unix
import databases
import bamtools
import structures
import clustering
import clusters
import filters
import output
## [SR CHANGE]
import virus
## [SR CHANGE]
import bkp

## FUNCTIONS ##

## CLASSES ##
class SV_caller():
    '''
    Structural variation (SV) caller 
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        self.mode = mode
        self.bam = bam
        self.normalBam = normalBam
        self.reference = reference
        self.refDir = refDir
        self.confDict = confDict
        self.outDir = outDir
        self.retrotransposonDb = None
        self.retrotransposonDbIndex = None


class SV_caller_long(SV_caller):
    '''
    Structural variation (SV) caller for long read sequencing data
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

    def call(self):
        '''
        Search for structural variants (SV) genome wide or in a set of target genomic regions
        '''
        ### 1. Create and index reference databases prior SV calling ##
        dbDir = self.outDir + '/databases'
        unix.mkdir(dbDir)

        self.retrotransposonDb, self.retrotransposonDbIndex = databases.buildRetrotransposonDb(self.refDir, self.confDict['transductionSearch'], dbDir)

        ### 2. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 3. Search for SV clusters in each bin ##
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        INS_clusters, DEL_clusters, left_CLIPPING_clusters, right_CLIPPING_clusters = zip(*pool.map(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        ### 4. Report clusters into output file
        # Create dictionary containing clusters 
        clusters = {}
        clusters['INS-CLUSTER'] = INS_clusters
        clusters['DEL-CLUSTER'] = DEL_clusters
        clusters['LEFT-CLIPPING-CLUSTER'] = left_CLIPPING_clusters 
        clusters['RIGHT-CLIPPING-CLUSTER']= right_CLIPPING_clusters
    
        # Write clusters
        output.writeClusters(clusters, self.outDir)

        ### 5. Do cleanup
        unix.rm([dbDir])

    def make_clusters_bin(self, window):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''

        ## 0. Set bin id and create bin directory ##
        ref, beg, end = window
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        binDir = self.outDir + '/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:
            eventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        step = 'COLLECT'
        SV_types = sorted(eventsDict.keys())

        counts = [str(len(eventsDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
    
        log.step(step, msg)

        ## 2. Organize all the SV events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior metaclustering'
        log.step(step, msg)

        ## Define bin database sizes 
        binSizes = [self.confDict['maxEventDist'], 1000, 10000, 100000, 1000000]

        ## Create bins
        eventsBinDb = structures.create_bin_database(ref, beg, end, eventsDict, binSizes)

        ## 3. Group events into SV metaclusters ##
        metaclustersBinDb = clusters.create_metaclusters(eventsBinDb, self.confDict)
        
        step = 'META-CLUSTERING'
        msg = 'Number of created metaclusters: ' + str(metaclustersBinDb.nbEvents()[0])
        log.step(step, msg)

        
class SV_caller_short(SV_caller):
    '''
    Structural variation (SV) caller for Illumina short read sequencing data
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

    def call(self):
        '''
        Search for structural variants (SV) genome wide or in a set of target genomic regions
        '''
        ### 1. Create and index reference databases prior SV calling ##
        dbDir = self.outDir + '/databases'
        unix.mkdir(dbDir)

        self.viralDb, self.viralDbIndex = databases.buildVirusDb(self.refDir, dbDir)

        ### 2. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 3. Search for SV clusters in each bin ##
        # Genomic bins will be distributed into X processes
        # TODO: mirar que pasa cuando tienes 2 dictionarios
        pool = mp.Pool(processes=self.confDict['processes'])
        metaclustersList = pool.map(self.make_clusters_bin, bins)
        pool.close()
        pool.join()

        '''
        for dictMetacluster in metaclustersList:
            for metacluster,d2 in dictMetacluster.items():
                print('METACLUSTER: ', str(metacluster) +' '+ str(len(metacluster.events)) +' '+ str(metacluster.ref) +' '+ str(metacluster.beg) +' '+ str(metacluster.end) +' '+ str(metacluster.intOrigin))
                for event in metacluster.events:
                        if event.type == 'DISCORDANT':
                            print (str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' ' + str(event.identity) + ' ' + str(event.side) + ' ' + str(event.sample))
                        else:
                            print (str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' None ' + str(event.clippedSide) + ' ' + str(event.sample))
                for k,v in d2.items():
                    print (str(k) + ' = ' + str(v))   
                    '''     
    
        # Write metaclusters
        output.writeMetaclusters(metaclustersList, self.outDir)

        ### 5. Do cleanup
        #unix.rm([dbDir])

    def make_clusters_bin(self, window):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''

        ## 0. Set bin id and create bin directory ##
        ref, beg, end = window
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        binDir = self.outDir + '/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            discordantEventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

        # b) Paired sample mode (tumour & matched normal)
        else:
            discordantEventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        step = 'COLLECT'
        SV_types = sorted(discordantEventsDict.keys())

        counts = [str(len(discordantEventsDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)

        '''
        ## PRINT 
        for events in discordantEventsDict.values():
            for event in events:
                print (str(event.type) +' '+ str(event.ref) +' '+ str(event.beg) +' '+ str(event.end) +' '+ str(event.mateSeq))
                DISCORDANT 19 33770596 33770697 None
                '''

        ## 2. Mark those discordant events corresponding to a RT insertion ##
        ## TODO

        ## 3. Mark those discordant events corresponding to a viral insertion ##

        # Create a list containing all discordant events:
        discordantEvents = discordantEventsDict['PLUS-DISCORDANT'] + discordantEventsDict['MINUS-DISCORDANT']

        # TODO: ADD RT!!!!

        # Initialize dictionary
        discordantEventsIdent = {}

        # a) Single sample mode
        if self.mode == "SINGLE":
            discordantEventsIdent = virus.is_virusSR(discordantEvents, self.bam, None, binDir, self.viralDbIndex)

        # b) Paired sample mode (tumour & matched normal)
        else:
            discordantEventsIdent = virus.is_virusSR(discordantEvents, self.bam, self.normalBam, binDir, self.viralDbIndex)

        ## 4. Organize identified events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior metaclustering'
        log.step(step, msg)

        ## Define bin database sizes 
        ## Big window sizes are needed for SR (see comments)
        binSizes = [self.confDict['maxEventDist'], 10000, 100000, 1000000]

        ## Create bins
        eventsBinDb = structures.create_bin_database(ref, beg, end, discordantEventsIdent, binSizes)

        '''
        ## PRINT 
        for a,b in eventsBinDb.data.items():
            print ('1 '+str(a) +' 2 '+str(b))
            for c,d in b.items():
                print ('3 '+ str(c) +' 4 '+ str(d))
                for e,f in d.items():
                    print ('5 '+ str(e) +' 6 '+ str(f))
                    
                    Above print yields the following:
                    1 100 2 {337718: {'PLUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f1b5060fc88>}}
                    3 337718 4 {'PLUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f1b5060fc88>}
                    5 PLUS-DISCORDANT-Hepatitis 6 <structures.events_bin object at 0x7f1b5060fc88>
                    1 10000 2 {3377: {'PLUS-DISCORDANT-Cowpox': <structures.events_bin object at 0x7f1b5060f6d8>, 'PLUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f1b5060f9e8>, 'PLUS-DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7f1b5060f908>}}
                    3 3377 4 {'PLUS-DISCORDANT-Cowpox': <structures.events_bin object at 0x7f1b5060f6d8>, 'PLUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f1b5060f9e8>, 'PLUS-DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7f1b5060f908>}
                    5 PLUS-DISCORDANT-Cowpox 6 <structures.events_bin object at 0x7f1b5060f6d8>
                    5 PLUS-DISCORDANT-Hepatitis 6 <structures.events_bin object at 0x7f1b5060f9e8>
                    5 PLUS-DISCORDANT-UNVERIFIED: 6 <structures.events_bin object at 0x7f1b5060f908>
                    '''

        ## 5. Group discordant events into clusters based on their identity ##
        discordantClustersDict = clusters.create_discordantClusters(eventsBinDb, self.confDict)

        discordantClustersBinDb = structures.create_bin_database(ref, beg, end, discordantClustersDict, binSizes)
        
        '''
        for a,b in discordantClustersBinDb.data.items():
            print ('1 '+str(a) +' 2 '+str(b))
            for c,d in b.items():
                print ('3 '+ str(c) +' 4 '+ str(d))
                for e,f in d.items():
                    print ('5 '+ str(e) +' 6 '+ str(f))
                    for cluster in f.events:
                        print ('7 '+ str(cluster))
                        for event in cluster.events:
                            print ('8 '+ event.type)
                            print ('9 '+ event.side)
                            '''
                            

        '''                 
        Above print yields the following:                    
        1 100 2 {}
        1 10000 2 {10545: {'PLUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f551ec92438>, 'PLUS-DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7f551ec92390>, 'MINUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f551ec92128>, 'MINUS-DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7f551ec92630>, 'MINUS-DISCORDANT-HBV': <structures.events_bin object at 0x7f551ec925f8>}}
        3 10545 4 {'PLUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f551ec92438>, 'PLUS-DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7f551ec92390>, 'MINUS-DISCORDANT-Hepatitis': <structures.events_bin object at 0x7f551ec92128>, 'MINUS-DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7f551ec92630>, 'MINUS-DISCORDANT-HBV': <structures.events_bin object at 0x7f551ec925f8>}
        5 PLUS-DISCORDANT-Hepatitis 6 <structures.events_bin object at 0x7f551ec92438>
        7 <clusters.DISCORDANT_cluster object at 0x7f551ec92240>
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        7 <clusters.DISCORDANT_cluster object at 0x7f551ec92668>
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        7 <clusters.DISCORDANT_cluster object at 0x7f551ec92470>
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        7 <clusters.DISCORDANT_cluster object at 0x7f551ec924e0>
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        5 PLUS-DISCORDANT-UNVERIFIED: 6 <structures.events_bin object at 0x7f551ec92390>
        7 <clusters.DISCORDANT_cluster object at 0x7f551ec926a0>
        8 DISCORDANT
        9 PLUS
        8 DISCORDANT
        9 PLUS
        5 MINUS-DISCORDANT-Hepatitis 6 <structures.events_bin object at 0x7f551ec92128>
        7 <clusters.DISCORDANT_cluster object at 0x7f551ec920b8>
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        8 DISCORDANT
        9 MINUS
        '''


        step = 'DISCORDANT-CLUSTERING'
        msg = 'Number of created discordant clusters: ' + str(discordantClustersBinDb.nbEvents()[0])
        log.step(step, msg)

        ## 6. Filter discordant metaclusters ##

        filters.filterClusters(discordantClustersBinDb, ['DISCORDANT'], self.confDict, self.bam)

        ## Remove those clusters that fail in one or more filters
        newDiscordantClustersDict = filters.applyFilters(discordantClustersBinDb)

        # vuelvo a hacer la bindb que contiene ya solo los clusters que pasaron los filtros
        discordantClustersBinDb = structures.create_bin_database(ref, beg, end, newDiscordantClustersDict, binSizes)

        ## 7. Making reciprocal clusters ##
        # TODO: AJUSTAR ESTOS PARAMETROS!!! (PASARLOS SI ESO COMO OPCION EN LOS ARGUMENTOS)
        reciprocalEventsDict = clustering.reciprocal(discordantClustersBinDb, 1, 1, 300)

        ## 8. Put reciprocal and independent clusters into bins ##

        reciprocalEventsBinDb = structures.create_bin_database(ref, beg, end, reciprocalEventsDict, binSizes)
        '''
        for a,b in reciprocalEventsBinDb.data.items():
            print ('1 '+str(a) +' 2 '+str(b))
            for c,d in b.items():
                print ('3 '+ str(c) +' 4 '+ str(d))
                for e,f in d.items():
                    print ('5 '+ str(e) +' 6 '+ str(f))
                    for event in f.events:
                        print ('7 '+ str(event))
                        print ('8 '+ event.side)
                        print ('9 '+ event.type)
                        print ('10 '+ event.identity)
        1 100 2 {1054575: {'DISCORDANT-Hepatitis': <structures.events_bin object at 0x7fc20a58d2e8>}}
        3 1054575 4 {'DISCORDANT-Hepatitis': <structures.events_bin object at 0x7fc20a58d2e8>}
        5 DISCORDANT-Hepatitis 6 <structures.events_bin object at 0x7fc20a58d2e8>
        7 <events.DISCORDANT object at 0x7fc20a582278>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        1 10000 2 {10545: {'DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7fc20a58d710>, 'DISCORDANT-Hepatitis': <structures.events_bin object at 0x7fc20a58d208>, 'DISCORDANT-HBV': <structures.events_bin object at 0x7fc20a58d198>}}
        3 10545 4 {'DISCORDANT-UNVERIFIED:': <structures.events_bin object at 0x7fc20a58d710>, 'DISCORDANT-Hepatitis': <structures.events_bin object at 0x7fc20a58d208>, 'DISCORDANT-HBV': <structures.events_bin object at 0x7fc20a58d198>}
        5 DISCORDANT-UNVERIFIED: 6 <structures.events_bin object at 0x7fc20a58d710>
        7 <events.DISCORDANT object at 0x7fc20a57ca20>
        8 PLUS
        9 DISCORDANT
        10 UNVERIFIED:
        7 <events.DISCORDANT object at 0x7fc20a57ccf8>
        8 PLUS
        9 DISCORDANT
        10 UNVERIFIED:
        7 <events.DISCORDANT object at 0x7fc20a5824a8>
        8 MINUS
        9 DISCORDANT
        10 UNVERIFIED:
        7 <events.DISCORDANT object at 0x7fc20a582a20>
        8 MINUS
        9 DISCORDANT
        10 UNVERIFIED:
        5 DISCORDANT-Hepatitis 6 <structures.events_bin object at 0x7fc20a58d208>
        7 <events.DISCORDANT object at 0x7fc20a57c940>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57c940>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57c978>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57c9e8>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cac8>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cb00>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cb38>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cba8>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cc50>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57ccc0>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cd30>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cda0>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cdd8>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57ce10>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57ceb8>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cef0>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cf28>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cf60>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a57cf98>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582048>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5820f0>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582128>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582160>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582240>
        8 PLUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582320>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5822e8>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582320>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582358>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582390>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5823c8>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582438>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582518>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5824e0>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582550>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5825c0>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582630>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5826a0>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5826d8>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582710>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582780>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582828>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582860>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5828d0>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a5829b0>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582a58>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582b00>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582ba8>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        7 <events.DISCORDANT object at 0x7fc20a582c88>
        8 MINUS
        9 DISCORDANT
        10 Hepatitis
        5 DISCORDANT-HBV 6 <structures.events_bin object at 0x7fc20a58d198>
        7 <events.DISCORDANT object at 0x7fc20a582588>
        8 MINUS
        9 DISCORDANT
        10 HBV
        7 <events.DISCORDANT object at 0x7fc20a582588>
        8 MINUS
        9 DISCORDANT
        10 HBV
        1 100000 2 {}
        1 1000000 2 {}
        '''
        # HASTA AQUI YO CREO QUE ESTA TODO BIEN!!!! A PARTIR DE AQUI HAY QUE REPASAR!!

        ## 8. Create metaclusters from reciprocal and independent clusters ##

        metaclustersBinDb = clusters.create_metaclusters(reciprocalEventsBinDb, self.confDict, self.bam, self.normalBam, self.mode)

        '''
        LO ULTIMO QUE HICE FUE RETORNAR EVENTS DE LA RECIPROCAL EN VEZ DE CLUSTERS, Y FUNCIONA, PERO HAY EN ALGUN MOMENTO QUE SE MEZCLAN LOS DE DISTINTO TIPO AL HACER LA RECIPROCAL, ASI QUE TENGO QUE REPASARLO!
        
        {'DISCORDANT-Hepatitis': [<events.DISCORDANT object at 0x7f73aeb5c208>, <events.DISCORDANT object at 0x7f73aeb5c780>, <events.DISCORDANT object at 0x7f73aeb55780>, <events.DISCORDANT object at 0x7f73aeb55a58>], 'DISCORDANT-UNVERIFIED:': [<events.DISCORDANT object at 0x7f73aeb5c208>, <events.DISCORDANT object at 0x7f73aeb5c780>, <events.DISCORDANT object at 0x7f73aeb55780>, <events.DISCORDANT object at 0x7f73aeb55a58>], 'DISCORDANT-HBV': [<events.DISCORDANT object at 0x7f73aeb5c2b0>, <events.DISCORDANT object at 0x7f73aeb5c2b0>, <events.DISCORDANT object at 0x7f73aeb55780>, <events.DISCORDANT object at 0x7f73aeb55a58>]}
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6390> 2 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6390> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6390> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> 5 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:7:2108:9935:21524/1 2 105457573 DISCORDANT UNVERIFIED: MINUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae63c8> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6400> 2 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6400> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6400> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> 5 2 105457243 105457870
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:7:2108:9935:21524/1 2 105457573 DISCORDANT UNVERIFIED: MINUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6438> HWI-ST672:129:D0DF0ACXX:8:2306:14996:36787/1 2 105457769 DISCORDANT UNVERIFIED: MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae6470> 2 2 105457243 105457720
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6470> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae6470> HWI-ST672:120:D0CF5ACXX:8:2104:17314:11147/1 2 105457619 DISCORDANT HBV MINUS
        METACLUSTER:  <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> 5 2 105457243 105457720
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:120:D0CF5ACXX:8:2308:15074:2709/2 2 105457243 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:129:D0DF0ACXX:8:2107:19576:151218/2 2 105457376 DISCORDANT UNVERIFIED: PLUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:120:D0CF5ACXX:8:2104:17314:11147/1 2 105457619 DISCORDANT HBV MINUS
        <clusters.DISCORDANT_cluster object at 0x7f73aeae64a8> HWI-ST672:120:D0CF5ACXX:8:2104:17314:11147/1 2 105457619 DISCORDANT HBV MINUS

        '''

        '''
        dictMetaclustersLEFT = bkp.analizeBkp(metaclustersBinDb, self.viralDb, self.reference, 'LEFT', binDir)
        dictMetaclustersRIGHT = bkp.analizeBkp(metaclustersBinDb, self.viralDb, self.reference, 'RIGHT', binDir)

        
        print ('LEFT')
        print (dictMetaclustersLEFT)

        print ('RIGHT')
        print (dictMetaclustersRIGHT)
        '''

        dictMetaclusters = bkp.analyzeMetaclusters(metaclustersBinDb, self.confDict, self.bam, self.normalBam, self.mode, self.viralDb, self.viralDbIndex, binDir)

        ### Do cleanup
        #unix.rm([binDir])

        return dictMetaclusters


