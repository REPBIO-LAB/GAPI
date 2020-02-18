'''
Module 'callers' - Contains classes and functions for calling variants from next generation sequencing data
'''

## DEPENDENCIES ##
# External
import sys
import multiprocessing as mp
import os
import pysam
import time
import libyay
# External
import subprocess

# Internal
import log
import unix
import databases
import formats
import bamtools
import structures
import events
import clusters
import output
import annotation
import bkp
import filters
import alignment
import gRanges
import virus
import clustering
import os

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
        self.repeatsBinDb = None


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

        ### 1. Create SV clusters 
        msg = '1. Create SV clusters'
        log.header(msg)
        allMetaclusters = self.make_clusters()
        
        ### 2. Annotate SV clusters intervals  
        msg = '2. Annotate SV clusters intervals'
        log.header(msg)

        # Create output directory
        annotDir = self.outDir + '/ANNOT/'
        unix.mkdir(annotDir)

        # Reference lengths, needed for repeats annotation
        refLengths = bamtools.get_ref_lengths(self.bam)
        
        # Define annotation steps
        steps = ['REPEAT']

        if self.confDict['annovarDir'] is not None:
            steps.append('GENE')

        # For each cluster type
        for SV_type in allMetaclusters:
            
            metaclusters = allMetaclusters[SV_type]
            annotation.annotate(metaclusters, steps, refLengths, self.refDir, self.confDict['annovarDir'], self.confDict['processes'], annotDir)

        # Remove annotation directory
        unix.rm([annotDir])

        ### 3. Determine what type of sequence has been inserted for INS metaclusters
        msg = '3. Determine what type of sequence has been inserted for INS metaclusters'
        log.header(msg)

        # Create output directory
        outDir = self.outDir + '/INS_TYPE/'
        unix.mkdir(outDir)

        if 'INS' in allMetaclusters:

            ## Infer insertion type
            clusters.INS_type_metaclusters(allMetaclusters['INS'], self.reference, refLengths, self.refDir, self.confDict['transductionSearch'], 1, outDir)

        # Remove output directory
        unix.rm([outDir])
            
        ### 4. Resolve structure for solo, partnered and orphan transductions
        msg = '4. Resolve structure for solo, partnered and orphan transductions'
        log.header(msg)
        
        if 'INS' in allMetaclusters:
            consensus = self.refDir + '/consensusDb.fa'
            transduced = self.refDir + '/transducedDb.fa.masked'

            # Create output directory
            outDir = self.outDir + '/STRUCTURE/'
            unix.mkdir(outDir)

            # Structure inference
            allMetaclusters['INS'] = clusters.structure_inference_parallel(allMetaclusters['INS'], consensus, transduced, self.confDict['transductionSearch'], self.confDict['processes'], outDir)

            # Remove output directory
            unix.rm([outDir])
        
        ### 5. Identify BND junctions
        msg = '5. Identify BND junctions'
        log.header(msg)
        allJunctions = []

        if 'BND' in allMetaclusters:
        
            ### Search for repeat, transduction or viral bridges
            # Create output directory
            outDir = self.outDir + '/BRIDGES/'
            unix.mkdir(outDir)   

            allMetaclusters['BND'] = clusters.search4bridges_metaclusters_parallel(allMetaclusters['BND'], 10000, 80, self.confDict['minSupportingReads'], 25, refLengths, self.refDir, self.confDict['processes'], outDir)

            ### Search for BND junctions
            allJunctions = clusters.search4junctions_metaclusters(allMetaclusters['BND'], refLengths, self.confDict['processes'], self.confDict['minSupportingReads'], 25)
            
            # Remove output directory
            unix.rm([outDir])

        ### 6. Apply second round of filtering 
        msg = '6. Apply second round of filtering'
        log.header(msg)
        filters2Apply = ['PERC-RESOLVED']
        metaclustersPass, metaclustersFailed = filters.filter_metaclusters(allMetaclusters, filters2Apply, self.confDict)
                
        ### 7. Report SV calls into output files
        msg = '7. Report SV calls into output files'
        log.header(msg)
        
        ##  7.1 Report INS
        if 'INS' in metaclustersPass:
            outFileName = 'INS_MEIGA.PASS.tsv'
            output.write_INS(metaclustersPass['INS'], outFileName, self.outDir)

        if 'INS' in metaclustersFailed:
            outFileName = 'INS_MEIGA.FAILED.2.tsv'
            output.write_INS(metaclustersFailed['INS'], outFileName, self.outDir)

        ## 7.2 Report BND junctions
        if allJunctions:
            outFileName = 'BND_junctions_MEIGA.PASS.tsv'
            output.write_junctions(allJunctions, outFileName, self.outDir)
        
    def make_clusters(self):
        '''
        Search for structural variant (SV) clusters 

        Output:
            1. metaclustersPass: dictionary containing one key per SV type and the list of metaclusters identified of each given SV type
        '''
        ### 1. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])
        
        ### 2. Search for SV clusters in each bin ##
        # Create output directory
        unix.mkdir(self.outDir + '/CLUSTER/')

        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        metaclustersPassList, metaclustersFailedList = zip(*pool.starmap(self.make_clusters_bin, bins))
        pool.close()
        pool.join()

        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        ### 3. Collapse metaclusters in a single dict and report metaclusters that failed filtering
        metaclustersPass = structures.merge_dictionaries(metaclustersPassList)
        metaclustersFailed = structures.merge_dictionaries(metaclustersFailedList)

        if 'INS' in metaclustersFailed:
            outFileName = 'INS_MEIGA.FAILED.1.tsv'
            output.write_INS(metaclustersFailed['INS'], outFileName, self.outDir)
        
        return metaclustersPass

    def make_clusters_bin(self, ref, beg, end):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''
        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)
        start = time.time()

        binDir = self.outDir + '/CLUSTER/' + binId
        unix.mkdir(binDir)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None, True)

        # b) Paired sample mode (tumour & matched normal)
        else:
            eventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        SV_types = sorted(eventsDict.keys())
        counts = [str(len(eventsDict[SV_type])) for SV_type in SV_types]
       
        step = 'COLLECT'
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)
        
        ## 2. Organize all the SV events into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize all the SV events into genomic bins prior metaclustering'
        log.step(step, msg)

        ## Define bin database sizes 
        minBinSize = min([self.confDict['maxInsDist'], self.confDict['maxBkpDist']])
        binSizes = [minBinSize, 1000, 10000, 100000, 1000000]

        ## Create bins
        eventsBinDb = structures.create_bin_database_interval(ref, beg, end, eventsDict, binSizes)

        ## 3. Group events into clusters ##
        step = 'CLUSTERING'
        msg = 'Group events into clusters' 
        log.step(step, msg)
        clustersBinDb = clusters.create_clusters(eventsBinDb, self.confDict)
        
        msg = 'Number of created clusters: ' + str(clustersBinDb.nbEvents()[0])
        log.step(step, msg)

        ## 4. Polish clusters ##
        step = 'POLISH'
        msg = 'Polish SV clusters' 
        log.step(step, msg)
        clusters.polish_clusters(clustersBinDb, self.confDict['minClusterSize'])

        ## 5. Group events into metaclusters ##
        step = 'META-CLUSTERING'
        msg = 'Group events into metaclusters' 
        log.step(step, msg)

        metaclusters = clusters.create_metaclusters(clustersBinDb)        
        msg = 'Number of created metaclusters: ' + str(len(metaclusters)) 
        log.step(step, msg)

        ## 6. Infer structural variant type ##
        step = 'SV-TYPE'
        msg = 'Infer structural variant type' 
        log.step(step, msg)

        outDir = binDir + '/SV_TYPE/' 
        metaclustersSVType = clusters.SV_type_metaclusters(metaclusters, self.confDict['minINDELlen'], self.confDict['technology'], outDir)
        
        # Do cleanup
        unix.rm([outDir])

        ## 7. Filter metaclusters ##
        step = 'FILTER'
        msg = 'Filter out metaclusters' 
        log.step(step, msg)
        filters2Apply = ['MIN-NBREADS', 'MAX-NBREADS', 'CV', 'SV-TYPE']
        metaclustersSVType, metaclustersSVTypeFailed = filters.filter_metaclusters(metaclustersSVType, filters2Apply, self.confDict)

        ## 8. Generate consensus event for SV metaclusters ##
        step = 'CONSENSUS'
        msg = 'Generate consensus event for SV metaclusters' 
        log.step(step, msg)

        targetSV = ['INS']
        outDir = binDir + '/CONSENSUS/' 
        clusters.create_consensus(metaclustersSVType, self.confDict, self.reference, targetSV, outDir)       

        ## 9. Lighten up metaclusters  ##
        ## Metaclusters passing all the filters
        clusters.lighten_up_metaclusters(metaclustersSVType)

        ## Filtered metaclusters
        clusters.lighten_up_metaclusters(metaclustersSVTypeFailed)
        
        # Do cleanup
        unix.rm([outDir, binDir])

        ## Print time taken to process bin
        end = time.time()
        time_taken = end - start
        msg = 'SV calling in bin: ' + binId + ' finished in ' + str(time_taken)
        log.info(msg)

        return metaclustersSVType, metaclustersSVTypeFailed


class SV_caller_short(SV_caller):
    '''
    Structural variation (SV) caller for Illumina short read sequencing data
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

    def call(self):
        '''
        Search for integrations genome wide or in a set of target genomic regions
        '''

        # TODO: DESILENCE
        '''
        annotDir = self.outDir + '/ANNOT/'
        refLengths = bamtools.get_ref_lengths(self.bam)
        self.annotations = annotation.load_annotations(['REPEATS', 'TRANSDUCTIONS'], refLengths, self.refDir, self.confDict['processes'], annotDir)
        '''

        ### 1. Create integration clusters 
        msg = '1. Create integration clusters'
        log.header(msg)      
        allMetaclusters = self.make_clusters()

        # TODO: DESILENCE
        '''
        ### 2. Annotate SV clusters intervals
        msg = '2. Annotate SV clusters intervals'
        log.header(msg)

        ## 5. Check if annotated retrotransposon on the reference genome at cluster intervals ##
        # COMMENT: This is temporary and will be incorporated into the filtering function at one point
        step = 'ANNOTATE-REPEATS'
        msg = 'Check if annotated retrotransposon on the reference genome at cluster intervals'
        log.step(step, msg)


        ## Create a list containing all discordant read pair events:
        allDiscordantClusters = []

        for eventType in allMetaclusters.keys():
            allDiscordantClusters.extend(allMetaclusters[eventType])


        ## Annotate
        buffer = 100
        annotation.repeats_annotation(allDiscordantClusters, self.annotations['REPEATS'], buffer)
        
        ## 6. Perform gene-based annotation with ANNOVAR of discordant read pair clusters ##
        # Do gene-based annotation step if enabled
        if self.confDict['annovarDir'] is not None:

            step = 'ANNOTATE'
            msg = 'Perform gene-based annotation with ANNOVAR of discordant read-pair clusters'
            log.step(step, msg)

            ## Annotate
            annotation.gene_annotation(allDiscordantClusters, self.confDict['annovarDir'], annotDir)

        # Remove annotation directory
        unix.rm([annotDir])
        '''

        ## 1.2 Create and index viral database
        #self.viralDb, self.viralDbIndex = databases.buildVirusDb(self.refDir, dbDir)
        
  
        # Report integrations calls into output files
        #output.write_DISCORDANT(discordantClusters, self.outDir)
        metaclustersList = list(allMetaclusters.values())
        outFileName = 'metaclusters.PASS.tsv'
        output.writeMetaclusters(metaclustersList, outFileName, self.outDir)


    def make_clusters(self):

        # TODO SR: Make analyse RT or viruses options

        # If viruses option is select, collect read name and sequence of discordant low quality reads from all bam refs
        # TODO: DESILENCE
        '''
        bins = bamtools.makeGenomicBins(self.bam, self.confDict['binSize'], None)

        # TODO: Pass more arguments
        pool = mp.Pool(processes=self.confDict['processes'])
        pool.starmap(self.collectSeq, bins)
        pool.close()
        pool.join()

        #TODO: merge fastas:
        filenames = []
        for bine in bins:
            window = self.outDir + '/FASTAS/' + str(bine[0]) +"_"+ str(bine[1])+"_"+str(bine[2])+".fasta"
            filenames.append(window)
        '''
        

        allFastas = self.outDir + "/allFastas.fasta"
        # TODO: DESILENCE
        '''
        with open(allFastas, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)


        # Remove fastas:
        fastasDir = self.outDir + '/FASTAS/'
        unix.rm([fastasDir])
        '''

        # bwa allFastas vs viralDb keep only mapped
        # TODO: usar una funcion ya hecha (o hacer una) y mirar si el -T vale para algo
        BAM = self.outDir + '/' + 'viralAligment' + '.bam'

        # TODO: DESILENCE
        '''
        err = open(self.outDir + '/align.err', 'w') 
        # TODO: set processes as argument
        command = 'bwa mem -Y -t 5 ' + self.confDict['viralDb'] + ' ' + allFastas + ' | samtools view -F 4 -b | samtools sort -O BAM   > ' + BAM
        status = subprocess.call(command, stderr=err, shell=True)

        if status != 0:
            step = 'ALIGN'
            msg = 'Alignment failed' 
            log.step(step, msg)

        command = 'samtools index ' + BAM
        status = subprocess.call(command, stderr=err, shell=True)

        # TODO: borro allfastas
        #unix.rm([allFastas])
        '''


        # Read bwa result and store in a dictionary
        bamFile = pysam.AlignmentFile(BAM, 'rb')

        iterator = bamFile.fetch()
        
        self.viralSeqs= {}
        # For each read alignment
        for alignmentObj in iterator:
            self.viralSeqs[alignmentObj.query_name] = alignmentObj.reference_name
        
        ### 1. Define genomic bins to search for SV ##
        bins = bamtools.binning(self.confDict['targetBins'], self.bam, self.confDict['binSize'], self.confDict['targetRefs'])

        ### 2. Search for SV clusters in each bin ##
        # Create output directory
        unix.mkdir(self.outDir + '/CLUSTER/')

        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        metaclustersPassList, metaclustersFailedList = zip(*pool.starmap(self.make_clusters_bin, bins))
        pool.close()
        pool.join()
        #print (libyay.massadd(listMates))
        #print (os.getpid())
        #libyay.cerrarGlob()

        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        #print (libyay.massadd(listMates))

        ### 3. Collapse metaclusters in a single dict and report metaclusters that failed filtering
        metaclustersPass = structures.merge_dictionaries(metaclustersPassList)
        metaclustersFailed = structures.merge_dictionaries(metaclustersFailedList)

        if metaclustersFailed:
            metaclustersFailedList = list(metaclustersFailed.values())
            outFileName = 'metaclusters.FAILED.tsv'
            output.writeMetaclusters(metaclustersFailedList, outFileName, self.outDir)

        return metaclustersPass

    def collectSeq(self, ref, binBeg, binEnd):
        '''
        Collect read names and sequences from reads below maxMAPQ
        '''
        # TODO SR: PASS this variables as argument
        #filterDuplicates = True
        maxMAPQ = 20
        checkUnmapped = True
        supplementary = True

        ## Initialize dictionary to store SV events
        eventsSeqDict = {}

        ## Open BAM file for reading
        bamFile = pysam.AlignmentFile(self.bam, "rb")

        ## Extract alignments
        iterator = bamFile.fetch(ref, binBeg, binEnd)
        
        # For each read alignment
        for alignmentObj in iterator:

            ### 1. Filter out alignments based on different criteria:
            MAPQ = int(alignmentObj.mapping_quality) # Mapping quality

            ## No query sequence available
            if alignmentObj.query_sequence == None:
                continue

            ## Aligments with MAPQ < threshold
            if (MAPQ > maxMAPQ):
                continue

            # TODO SR: make this work
            ## Duplicates filtering enabled and duplicate alignment
            #if (confDict['filterDuplicates'] == True) and (alignmentObj.is_duplicate == True):
                #continue

            # Filter supplementary alignments if TRUE. (Neccesary to avoid pick supplementary clipping reads while adding to discordant clusters in short reads mode)
            if supplementary == False and alignmentObj.is_supplementary == True:
                continue
            
            ## 4. Collect DISCORDANT

            if not alignmentObj.is_proper_pair:

                # Pick only those sequences that are unmmapped or with mapping quality < maxMAPQ
                if checkUnmapped == True:
                    if (alignmentObj.is_unmapped == True) or (MAPQ < maxMAPQ):
                        eventsSeqDict[alignmentObj.query_name]=alignmentObj.query_sequence
                else:
                    if MAPQ < maxMAPQ:
                        eventsSeqDict[alignmentObj.query_name]=alignmentObj.query_sequence
            
        ## Close 
        bamFile.close()

        # Write FASTA:
        fastasDir = self.outDir + '/FASTAS/'
        unix.mkdir(fastasDir)

        seqsFastaObj= formats.FASTA()
        seqsFastaObj.seqDict = eventsSeqDict

        del eventsSeqDict


        outputFasta = self.outDir + '/FASTAS/' + str(ref) +"_"+ str(binBeg) +"_"+ str(binEnd) +".fasta"
        seqsFastaObj.write(outputFasta)

        # return sv candidates
        return



    def make_clusters_bin(self, ref, beg, end):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''

        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)
        start = time.time()

        binDir = self.outDir + '/CLUSTER/' + binId
        unix.mkdir(binDir)

        ## 1. Search for integration candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            discordantDict = bamtools.collectDISCORDANT(ref, beg, end, self.bam, self.confDict, None, True, self.viralSeqs)

        #   TODO: ADAPT FOR PAIRED!
        # b) Paired sample mode (tumour & matched normal)
        #else:
            #discordantDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        SV_types = sorted(discordantDict.keys())
        counts = [str(len(discordantDict[SV_type])) for SV_type in SV_types]
        
        step = 'COLLECT'
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)
        
        ## 2. Discordant read pair identity ##
        ## Determine identity
        # TODO: HACER EL DICT COMO ESTABA ANTES
        if self.mode == "SINGLE":
            # TODO: DESILENCE AND DO THIS OLNY FOR RT!!!
            #discordantsIdentity = events.determine_discordant_identity(discordantDict['DISCORDANT'], self.annotations['REPEATS'], self.annotations['TRANSDUCTIONS'],self.bam, None, binDir, self.confDict['viralDb'])
            discordantsIdentity = events.determine_discordant_identity(discordantDict['DISCORDANT'], None, None,self.bam, None, binDir)
        #else:
            # TODO: DESILENCE
            #discordantsIdentity = events.determine_discordant_identity(discordantDict['DISCORDANT'], self.annotations['REPEATS'], self.annotations['TRANSDUCTIONS'],self.bam, None, binDir, self.confDict['viralDb'])
            #discordantsIdentity = events.determine_discordant_identity(discordantDict['DISCORDANT'], None, None,self.bam, None, binDir, self.confDict['viralDb'])

        for discirdant in discordantDict['DISCORDANT']:
            print ('DISCORDAAAAAAAAANT' +' '+ str(discirdant.beg) +' '+ str(discirdant.end)  +' '+ str(discirdant.orientation) +' '+ str(discirdant.pair) +' '+ str(discirdant.readName)  +' '+ str(discirdant.identity))
        
        # TODO: FILTER LOS QUE SON NONE!!!
        step = 'IDENTITY'
        SV_types = sorted(discordantsIdentity.keys())
        counts = [str(len(discordantsIdentity[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events per identity in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)

        ## 3. Organize discordant read pairs into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize discordant read pairs into genomic bins prior clustering'
        log.step(step, msg)

        ## NOTE: bigger window sizes are needed for SR (see comments, ask Eva where are the comments?)
        ## Define bin database sizes 
        minBinSize = min([self.confDict['maxInsDist'], self.confDict['maxBkpDist']])
        binSizes = [minBinSize, 1000, 10000, 100000, 1000000]

        ## Create bins
        discordantsBinDb = structures.create_bin_database_interval(ref, beg, end, discordantsIdentity, binSizes)

        ## 4. Group discordant read pairs into clusters based on their mate identity ##
        buffer = 100
        discordantClustersDict = clusters.create_discordantClusters(discordantsBinDb, self.confDict['minClusterSize'], buffer)

        '''
        # TODO: Remove this print
        for dis in discordantClustersDict.values():
            for discirdant in dis:
                print ('DISCORDAAAAAAAAANT')
                print (discirdant.beg)
                print (discirdant.end)
                print (discirdant.clusterType)
        '''
    
        step = 'DISCORDANT-CLUSTERING'
        SV_types = sorted(discordantClustersDict.keys())
        counts = [str(len(discordantClustersDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of created discordant clusters in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)

        counter1 = 0
        
        for lista in discordantClustersDict.values():
            for clustereva in lista:
                for ebent in clustereva.events:
                    print (ebent.readName)
                    counter1 += 1 
        print ('discordantClustersDict0 ' + str(counter1))

        #if counts == []:
            #unix.rm([binDir])
            #return None

        '''
        ## 5. Check if annotated retrotransposon on the reference genome at cluster intervals ##
        # COMMENT: This is temporary and will be incorporated into the filtering function at one point
        step = 'ANNOTATE-REPEATS'
        msg = 'Check if annotated retrotransposon on the reference genome at cluster intervals'
        log.step(step, msg)

        ## Create a list containing all discordant read pair events:
        allDiscordantClusters = []

        for eventType in discordantClustersDict.keys():
            allDiscordantClusters.extend(discordantClustersDict[eventType])

        ## Annotate
        buffer = 100
        annotation.repeats_annotation(allDiscordantClusters, self.annotations['REPEATS'], buffer)
        
        ## 6. Perform gene-based annotation with ANNOVAR of discordant read pair clusters ##
        # Do gene-based annotation step if enabled
        if self.confDict['annovarDir'] is not None:

            step = 'ANNOTATE'
            msg = 'Perform gene-based annotation with ANNOVAR of discordant read-pair clusters'
            log.step(step, msg)

            ## Annotate
            annotDir = binDir + '/ANNOT/' 
            annotation.gene_annotation(allDiscordantClusters, self.confDict['annovarDir'], annotDir)
        '''

        ### Do cleanup
        #unix.rm([binDir])

        #return discordantClustersDict


        ## 5. Filter discordant clusters ##
        step = 'FILTER'
        msg = 'Filter out metaclusters' 
        log.step(step, msg)
        filters2Apply = ['MIN-NBREADS', 'MAX-NBREADS', 'AREAMAPQ', 'AREASMS']
        discordantClustersDict, discordantClustersDictFailed = filters.filter_clusters(discordantClustersDict, filters2Apply, self.confDict, self.bam)

        counter2 = 0
        for lista in discordantClustersDict.values():
            for clustereva in lista:
                for ebent in clustereva.events:
                    print (ebent.readName)
                    counter2 += 1
        print ('discordantClustersDictAfterFilter ' + str(counter2))

        counter3 = 0
        for lista1 in discordantClustersDictFailed.values():
            for clustereva1 in lista1:
                for ebent1 in clustereva1.events:
                    print (ebent1.readName)
                    counter3 += 1
        print ('discordantClustersDictFailed ' + str(counter3))

        '''
        # TODO: Remove this print
        for dis in discordantClustersDict.values():
            for discirdant in dis:
                print ('DISCORDAAAAAAAAANT AFTER FILTERRRIIIING')
                print (discirdant.beg)
                print (discirdant.end)
                print (discirdant.clusterType)
        '''
        
        '''
        Eva will further polish next steps!

        ## 6. Filter discordant metaclusters ##

        filters.filterClusters(discordantClustersBinDb, ['DISCORDANT'], self.confDict, self.bam)

        ## Remove those clusters that fail in one or more filters
        newDiscordantClustersDict = filters.applyFilters(discordantClustersBinDb)

        step = 'DISCORDANT-FILTERING'
        SV_types = sorted(discordantClustersDict.keys())

        counts = [str(len(discordantClustersDict[SV_type])) for SV_type in SV_types]
        msg = '[DISCORDANT-FILTERING] Number of created discordant clusters after filtering in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.subHeader(msg)

        if counts == []:
            unix.rm([binDir])
            return None
        '''


        ## 6. Organize discordant clusters in bin database structure ##
        discordantClustersBinDb = structures.create_bin_database_interval(ref, beg, end, discordantClustersDict, binSizes)
        discordantClustersFailedBinDb = structures.create_bin_database_interval(ref, beg, end, discordantClustersDictFailed, binSizes)

        ## 7. Make reciprocal clusters ##
        # TODO: AJUSTAR ESTOS PARAMETROS!!! (PASARLOS SI ESO COMO OPCION EN LOS ARGUMENTOS)
        reciprocalClustersDict = clustering.reciprocal(discordantClustersBinDb, 1, 1, 300)
        reciprocalClustersFailedDict = clustering.reciprocal(discordantClustersFailedBinDb, 1, 1, 300)

        counter4 = 0
        for lista2 in reciprocalClustersDict.values():
            for clustereva2 in lista2:
                for ebent2 in clustereva2.events:
                    print (ebent2.readName)
                    counter4 += 1
        print ('reciprocalClustersDict ' + str(counter4))

        '''
        # TODO: Remove this print
        for dis in reciprocalClustersDict.values():
            for discirdant in dis:
                print ('DISCORDAAAAAAAAANT AFTER RECIPROCAAAL')
                print (discirdant.beg)
                print (discirdant.end)
                print (discirdant.clusterType)
        '''

        ## 8. Organize reciprocal and independent discordant clusters in bin database structure ##

        metaclusters=[]
        metaclustersFailed=[]
        reciprocalClustersBinDb = structures.create_bin_database_interval(ref, beg, end, reciprocalClustersDict, binSizes)
        reciprocalClustersFailedBinDb = structures.create_bin_database_interval(ref, beg, end, reciprocalClustersFailedDict, binSizes)
        buffer=300

        identities = set([iden.split('-')[2] for iden in reciprocalClustersBinDb.eventTypes])
        identitiesFailed = set([iden.split('-')[2] for iden in reciprocalClustersFailedBinDb.eventTypes])


        # Create metaclusters taking into account viral identities eventType.
        for iden in identities:
            currentEventTypes = [eventType for eventType in reciprocalClustersBinDb.eventTypes if (iden in eventType)]
            metaclusters.extend(clusters.create_discordant_metaclusters(reciprocalClustersBinDb, currentEventTypes))


        for idenF in identitiesFailed:
            currentEventTypes = [eventType for eventType in reciprocalClustersFailedBinDb.eventTypes if (idenF in eventType)]
            metaclustersFailed.extend(clusters.create_discordant_metaclusters(reciprocalClustersFailedBinDb, currentEventTypes))


        # HASTA AQUI YO CREO QUE ESTA TODO BIEN!!!! A PARTIR DE AQUI HAY QUE REPASAR!!

        ## 9. Create metaclusters from reciprocal and independent discordant clusters ##

        # Mirar aqui pq habra que ajustar varios parametros

        #metaclusters = clusters.create_metaclusters(reciprocalClustersBinDb)
        #metaclustersFailed = clusters.create_metaclusters(reciprocalClustersFailedBinDb)


        step = 'META-CLUSTERING'
        #msg = '[META-CLUSTERING] Number of created metaclusters: ' + str(metaclustersBinDb.nbEvents()[0])
        #log.subHeader(msg)

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



        
        dictMetaclustersLEFT = bkp.analizeBkp(metaclustersBinDb,'/mnt/lustre/scratch/home/usc/mg/eal/results/MEIGA_ShortReads/consensusHBV_goodHeaders.fa', self.reference, 'LEFT', binDir)
        dictMetaclustersRIGHT = bkp.analizeBkp(metaclustersBinDb, '/mnt/lustre/scratch/home/usc/mg/eal/results/MEIGA_ShortReads/consensusHBV_goodHeaders.fa', self.reference, 'RIGHT', binDir)


        print ('LEFT')
        print (dictMetaclustersLEFT)

        print ('RIGHT')
        print (dictMetaclustersRIGHT)
        '''

        ## 10. Analyse metaclusters features and add supporting clipping reads ##
        dictMetaclusters = bkp.analyzeMetaclusters(metaclusters, self.confDict, self.bam, self.normalBam, self.mode, binDir)
        dictMetaclustersFailed = bkp.analyzeMetaclusters(metaclustersFailed, self.confDict, self.bam, self.normalBam, self.mode, binDir)

        metaclustersSVType = {}
        metaclustersSVTypeFailed = {}

        # TODO: Mirar si aqui solo pueden ser discordant
        if dictMetaclusters:
            metaclustersSVType['DISCORDANT'] = list(dictMetaclusters.keys())
        else:
            metaclustersSVType['DISCORDANT'] = []

        if dictMetaclustersFailed:
            metaclustersSVTypeFailed['DISCORDANT'] = list(dictMetaclustersFailed.keys())
        else:
            metaclustersSVTypeFailed['DISCORDANT'] = []

        '''
        # TODO: Remove this print
        for metacluster in dictMetaclusters.keys():
                print ('METACLSUTERSSSSSS AFTER BKP!!!')
                print (metacluster.beg)
                print (metacluster.end)
                print (metacluster.events)
                print (metacluster.subclusters)
        '''

        ### Do cleanup
        unix.rm([binDir])

        ## 11. Lighten up metaclusters  ##
        ## Metaclusters passing all the filters
        clusters.lighten_up_metaclusters(metaclustersSVType)

        ## Filtered metaclusters
        clusters.lighten_up_metaclusters(metaclustersSVTypeFailed)

        return metaclustersSVType, metaclustersSVTypeFailed

class SV_caller_sureselect(SV_caller):
    '''
    Structural variation (SV) caller for Illumina sureselect data targetering source element´s downstream regions
    '''
    def __init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir):

        SV_caller.__init__(self, mode, bam, normalBam, reference, refDir, confDict, outDir)

    def call(self):
        '''
        Search for structural variants (SV) genome wide or in a set of target genomic regions
        '''

        ### 1. Create bed file containing source element´s transduced intervals 
        tdDir = self.outDir + '/TRANSDUCED/'
        unix.mkdir(tdDir)
        sourceBed = self.refDir + '/srcElements.bed'
        transducedPath = databases.create_transduced_bed(sourceBed, 10000, tdDir)
                
        ### 2. Define genomic bins to search for SV (will correspond to transduced areas)
        bins = bamtools.binning(transducedPath, None, None, None)

        ## Organize bins into a dictionary
        self.rangesDict = gRanges.rangeList2dict(bins)

        ### 3. Associate to each bin the src identifier
        BED = formats.BED()
        BED.read(transducedPath, 'List', None)   

        for index, coords in enumerate(bins):
            coords.append(BED.lines[index].optional['cytobandId'])
        
        unix.rm([tdDir])

        ### 4. Search for SV clusters in each bin 
        # Genomic bins will be distributed into X processes
        pool = mp.Pool(processes=self.confDict['processes'])
        clusterPerSrc = pool.starmap(self.make_clusters_bin, bins)
        pool.close()
        pool.join()

        # Convert into dictionary
        clusterPerSrcDict = {srcId:clusters for srcId,clusters in clusterPerSrc}

        ### 5. Write calls to file
        ## 5.1 Transduction counts per source element
        output.write_tdCounts_surelect(clusterPerSrcDict, self.outDir)

        ## 5.2 Transduction calls
        output.write_tdCalls_surelect(clusterPerSrcDict, self.outDir)

    def make_clusters_bin(self, ref, beg, end, srcId):
        '''
        Search for structural variant (SV) clusters in a genomic bin/window
        '''

        ## 0. Set bin id and create bin directory ##
        binId = '_'.join([str(ref), str(beg), str(end)])
        msg = 'SV calling in bin: ' + binId
        log.subHeader(msg)

        ## 1. Search for SV candidate events in the bam file/s ##
        # a) Single sample mode
        if self.mode == "SINGLE":
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None, True)

        # b) Paired sample mode (tumour & matched normal)
        else:
            eventsDict = bamtools.collectSV_paired(ref, beg, end, self.bam, self.normalBam, self.confDict)

        step = 'COLLECT'
        SV_types = sorted(eventsDict.keys())
        counts = [str(len(eventsDict[SV_type])) for SV_type in SV_types]
        msg = 'Number of SV events in bin (' + ','.join(['binId'] + SV_types) + '): ' + '\t'.join([binId] + counts)
        log.step(step, msg)
        
        ## 2. Organize discordant read pairs into genomic bins prior clustering ##
        step = 'BINNING'
        msg = 'Organize discordant read pairs into genomic bins prior clustering'
        log.step(step, msg)

        ## Define bin database sizes 
        ## Note: bigger window sizes are needed for SR (see comments, ask Eva where are the comments?)
        binSizes = [500, 1000, 10000, 100000, 1000000]

        ## Create bins
        discordantsBinDb = structures.create_bin_database_interval(ref, beg, end, eventsDict, binSizes)

        ## 3. Group discordant read pairs into clusters based on reciprocal overlap ##
        step = 'CLUSTERING'
        msg = 'Group discordant read pairs into clusters based on reciprocal overlap'
        log.step(step, msg)

        buffer = 150
        discordantClustersDict = clusters.create_discordantClusters(discordantsBinDb, self.confDict['minClusterSize'], buffer)

        ## 4. Do an extra clustering step based on mate position ##
        step = 'GROUP-BY-MATE'
        msg = 'Group discordant read pairs based on mate position'
        log.step(step, msg)

        ## Compute reference lengths
        refLengths = bamtools.get_ref_lengths(self.bam)

        ## Make groups
        discordantClustersDict = clusters.extra_clustering_by_matePos(discordantClustersDict, refLengths, self.confDict['minClusterSize'])

        ## 5. Filter out those clusters over NOT target reference ##
        step = 'FILTER-REF'
        msg = 'Filter out those clusters over NOT target reference'
        log.step(step, msg)
        filteredDiscordants = filters.filter_discordant_mate_ref(discordantClustersDict['DISCORDANT'], self.confDict['targetRefs'])

        ## 6. Filter out those clusters whose mates aligns over any source element downstream region ##
        step = 'FILTER-DOWSTREAM'
        msg = 'Filter out those clusters whose mates aligns over any source element downstream region'
        log.step(step, msg)
        filteredDiscordants = filters.filter_discordant_mate_position(filteredDiscordants, self.rangesDict, 10000)        
        
        ## 7. Filter out those clusters based on average MAPQ for mate alignments ##
        step = 'FILTER-MATE-MAPQ'
        msg = 'Filter out those clusters based on average MAPQ for mate alignments'
        log.step(step, msg)
        filteredDiscordants = filters.filter_discordant_mate_MAPQ(filteredDiscordants, 20, self.bam)

        return [srcId, filteredDiscordants]  

