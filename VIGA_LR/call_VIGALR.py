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
import subprocess
import statistics
import itertools
from collections import Counter
# TEMP BORRAR
import Bio.SeqUtils
from Bio.SeqUtils import lcc
from VIGA_SR import bamtools_VIGASR
from VIGA_LR import virusLR

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
import assembly

#from modules.callers import SV_caller
from callers import SV_caller

## FUNCTIONS ##
# Multiprocessing lock as global variable. Useful for safely writting from different processes to the same output file.
def init(l):
    global lock
    lock = l

## CLASSES ##
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

        # TODO: Put annotations as an option as it's done in short reads
        '''
        ### 2. Annotate SV clusters intervals  
        msg = '2. Annotate SV clusters intervals'
        log.header(msg)

        # Load annotations
        self.load_annotations()

        # Empty list harcoded, because it is not needed por viruses and it crashes otherwise.
        self.annotations['TRANSDUCTIONS'] = []

        # Create output directory
        annotDir = self.outDir + '/ANNOT/'
        unix.mkdir(annotDir)
        
        # Define annotation steps
        steps = ['REPEAT']

        if self.confDict['annovarDir'] is not None:
            steps.append('GENE')

        # For each cluster type
        for SV_type in allMetaclusters:
            
            metaclusters = allMetaclusters[SV_type]
            annotation.annotate(metaclusters, steps, self.annotations, self.confDict['annovarDir'], annotDir)

        # Remove annotation directory
        unix.rm([annotDir])
        '''
        # NOTE: For the moment, annotations is a empty dictionary. Change this!!!
        self.annotations = {}
        self.annotations['REPEATS'] = []
        self.annotations['TRANSDUCTIONS'] = []
        self.annotations['EXONS'] = []

        # 1.a Write sequences of all events of metaclusters together in a fasta file:
        # If SV_type == INS, write only the inserted sequence of each event
        # If SV_type == BND, write only the clipped sequence of each event

        selectMETAoutDir = self.outDir + '/SELECT_METACLUSTERS/'
        unix.mkdir(selectMETAoutDir)

        tupleList = []

        for SV_type, metaclusters in allMetaclusters.items():
            for metacluster in metaclusters:
            
                ## Add to the list
                fields = (metacluster, selectMETAoutDir)
                tupleList.append(fields)
        
        l = mp.Lock()
        
        # Create FASTA file with INS or CLIPPING sequence of each event of each metacluster.

        pool = mp.Pool(processes=self.confDict['processes'], initializer=init, initargs=(l,))
        pool.starmap(virusLR.supportingEventsFasta, tupleList)
        pool.close()
        pool.join()

        # 1.b Alig against viral db and filter hits.
        # Have a dictionary with only those events that have identity after the filters

        eventsFasta = selectMETAoutDir + "/events.fasta"
        eventsIdentity, allHits_viral = virusLR.checkEventsIdentity(eventsFasta, self.confDict['viralDb'], self.confDict['processes'], selectMETAoutDir)

        '''
        # No se si esto me hara falta en algun momento pero creo que de momento no...
        for SV_type, metaclusters in allMetaclusters.items():
            for metacluster in metaclusters:
                # TODO: PONER EN PARALELO!!!!
                metacluster.insertHits = []
                for query in allHits_viral.keys():
                    if any(query == event.readName for event in metacluster.events):
                        metacluster.insertHits.append(allHits_viral[query])
                        # Asi ya no tengo que hacer mas alineamientos
                        # No se si esto me valdra para algo...
        '''

        unix.rm([selectMETAoutDir])

        # 1.d Put identity to only those metaclusters that are mapped and PASS filters, all the others identity will remain None

        tupleList = []

        for SV_type, metaclusters in allMetaclusters.items():
            for metacluster in metaclusters:
            
                ## Add to the list
                fields = (SV_type, metacluster, eventsIdentity, allHits_viral, self.outDir)
                tupleList.append(fields)
        
        l = mp.Lock()
        pool = mp.Pool(processes=self.confDict['processes'])
        # Ponerle la identity a los eventos
        allMetaclustersIdentitySeparated, fastaDict = zip(*pool.starmap(virusLR.check_metaclusters_identity_allHits, tupleList))
        pool.close()
        pool.join()

        #print ('fastaDict ' + str(fastaDict))

        del allMetaclusters

        allMetaclustersIdentity = structures.merge_dictionaries(allMetaclustersIdentitySeparated)

        
        # Analyse those with identity proportion < 0.70 >= 0.50
        # Merge dictionary
        allfastaDict = {}
        for dicti in fastaDict:
            if dicti != {}:
                for key,value in dicti.items():
                    allfastaDict[key] = value
                

        # Perform bwa and minimap2 alignments for those sequences in fastaDict.
        # Write fasta
        polishDir = self.outDir + '/lowIdentitiesProportion'
        unix.mkdir(polishDir)
        fastaObj = formats.FASTA()
        fastaObj.seqDict = allfastaDict
        lowIdentityFasta = polishDir + '/lowidentity.fa'
        fastaObj.write(lowIdentityFasta)

        ## 3.3 Align consensus inserted sequences into the viral database
        msg = '3.3 Align consensus inserted sequences into the viral database'
        log.info(msg)  
        
        # bwa alignment
        SAM_viral = alignment.alignment_bwa(lowIdentityFasta, self.confDict['viralDb'], 'lowidentity', self.confDict['processes'], polishDir)
        print ('SAM_viral ' + str(SAM_viral))

        BAM_sorted = bamtools.SAM2BAM(SAM_viral, polishDir)
        # NOTE: Vuelve a pasar los filtros pero esta vez de la consenso. Vuelvo a filtrar (pero ahora la consenso)
        minParcialMatchVirus = 0
        minTotalMatchVirus = 125 # >= 125
        maxMatchCheckMAPQVirus = 45 # Cambio 28/08 15:54
        minMAPQVirus = 0
        maxBasePercVirus = 60 # <= 60 (no mucho pero algo filtra)
        minLccVirus = 1.57 # > 1.7 Cambio 28/08 16:00

        # Filter bwa alignment
        # Keys: readName_metacluster.beg = identity (solo aquellos que pasaron los filtros)
        if pysam.idxstats(BAM_sorted) != '*\t0\t0\t0\n':
            eventsIdentity2 = bamtools_VIGASR.filterBAM2FastaDict_LR(BAM_sorted, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='LR')
        else:
            eventsIdentity2 = {}

        ## Convert SAM to PAF
        PAF_viral = alignment.sam2paf(SAM_viral, 'lowidentity', polishDir)
        
        # minimap2 alignment
        SAM_viral_minimap = virusLR.alignment_minimap2_SAM(lowIdentityFasta, self.confDict['viralDb'], 'lowidentity_minimap', self.confDict['processes'], polishDir)
        print ('SAM_viral_minimap ' + str(SAM_viral_minimap))

        BAM_sorted_minimap = bamtools.SAM2BAM(SAM_viral_minimap, polishDir)

        # NOTE: Vuelve a pasar los filtros pero esta vez de la consenso. Vuelvo a filtrar (pero ahora la consenso)
        minParcialMatchVirus = 0
        minTotalMatchVirus = 125 # >= 125
        maxMatchCheckMAPQVirus = 83 # cambio 28/08 12:53, cambio 28/08 15:29
        minMAPQVirus = 0
        maxBasePercVirus = 60 # <= 60 (no mucho pero algo filtra)
        minLccVirus = 1.7 # > 1.7

        PAF_viral_minimap = alignment.sam2paf(SAM_viral_minimap, 'lowidentity_minimap', polishDir) # cambio 12:10

        allHits_viral_minimap = alignment.organize_hits_paf(PAF_viral_minimap)


        # Filter minimap2 alignment
        # Keys: readName_metacluster.beg = identity (solo aquellos que pasaron los filtros)
        if pysam.idxstats(BAM_sorted_minimap) != '*\t0\t0\t0\n':
            eventsIdentity3 = bamtools_VIGASR.filterBAM2FastaDict_LR(BAM_sorted_minimap, minTotalMatchVirus, minParcialMatchVirus, maxMatchCheckMAPQVirus, minMAPQVirus, maxBasePercVirus, minLccVirus, mode='LR', allHits_viral=allHits_viral_minimap)
        else:
            eventsIdentity3 = {}



        ## Organize hits according to their corresponding metacluster
        allHits_viral = alignment.organize_hits_paf(PAF_viral)


        # Put identity tu metacluster if the consensus match with db.
        # TODO: Save the consensus since it doesnt have to be done again
        print ('allMetaclustersIdentity0 ' + str(allMetaclustersIdentity))
        print ('eventsIdentity2 ' + str(eventsIdentity2))
        print ('eventsIdentity3 ' + str(eventsIdentity3))

        # Based in their consensus and in previous alignments with bwa and minimap, give identity to these metaclusters that don't have it:
        for metaclusterList in allMetaclustersIdentity.values():
            for metacluster in metaclusterList:
                print ('meta1 ' + str(metacluster))
                # TODO: Or identity == None.
                if 'IDENTITY' not in metacluster.SV_features.keys():
                    print ('IDENTITY1')
                    SV_TypeTested = None
                    identity = []
                    specificIdentity = []
                    bkpProximity = None
                    #if any([event.readName in allHits_viral.keys() for event in metacluster.events]) or any([event.readName in allHits_viral_minimap.keys() for event in metacluster.events]): # aqui estaran solo los que tienen la low identity
                    # If the alignment of the event.read passed the filters:
                    if any([(event.readName + '_'+ str(metacluster.beg)) in eventsIdentity2.keys() for event in metacluster.events]) or any([(event.readName + '_'+ str(metacluster.beg)) in eventsIdentity3.keys() for event in metacluster.events]):
                        
                        print ('IDENTITY2')
                        '''
                        identities = list(itertools.chain(*[event.identity for event in metacluster.events if event.identity != None]))
                        if identities:
                            countIdentities = Counter(identities)
                            metacluster.SV_features['IDENTITY'] = countIdentities.most_common(3)
                            specificIdentities = list(itertools.chain(*[event.specificIdentity for event in metacluster.events if event.identity != None]))
                            if specificIdentities:
                                countSpecificIdentities = Counter(specificIdentities)
                                metacluster.SV_features['SPECIDENTITY'] = countSpecificIdentities.most_common(3)
                        '''
                        # TODO: Do this in another way to find the event that is the key
                        # Find the event that was used as consensus name
                        # TODO: maybe previous if inst neccesary, as eventsIdentity is checked here as well.
                        for event in metacluster.events:
                            nombre =  event.readName + '_'+ str(metacluster.beg)
                            # Check bwa alignment
                            if (nombre in allHits_viral.keys() and nombre in eventsIdentity2.keys()):
                                for alig in allHits_viral[nombre].alignments:
                                    splitedtName = alig.tName.split('|')
                                    identity.append(splitedtName[0])
                                    if len(splitedtName) > 1:
                                        specificIdentity.append(splitedtName[1])
                                
                                # Check if it's near bkp
                                if 'BND' in allMetaclustersIdentity.keys():
                                    print ('BND1')
                                    
                                    if metacluster in allMetaclustersIdentity['BND']:
                                        #begsEnds = []
                                        SV_TypeTested = True
                                        print ('BND2')

                                        buffer = 651
                                        bkpProximity = virusLR.checkBkpProximity_wSide(event, nombre, allHits_viral, buffer)

                                        if bkpProximity != None:
                                            metacluster.SV_features['IDENTITY']  = list(set(identity))
                                            metacluster.SV_features['SPECIDENTITY'] = list(set(specificIdentity))
                                            print ('metacluster.SV_featuresmarte1 ' + str(metacluster.SV_features['IDENTITY']) + ' ' + str(metacluster.beg))

                                # Check if it's near bkp
                                if 'INS_noOverlap' in allMetaclustersIdentity.keys():
                                    print ('allMetaclustersIdentity1 ' + str(allMetaclustersIdentity['INS_noOverlap']))
                                    if metacluster in allMetaclustersIdentity['INS_noOverlap']:
                                        begsEnds = []
                                        SV_TypeTested = True

                                        buffer = 250
                                        bkpProximity = virusLR.checkBkpProximity_wSide(event, nombre, allHits_viral, buffer)
                                        event.bkpProximity = bkpProximity
                                        
                                        if bkpProximity != None:
                                            #if (event.clippedSide == 'right' and bkpProximity < 250) or (event.clippedSide == 'left' and bkpProximity >= (ali_query_Len - 250)):
                                            metacluster.SV_features['IDENTITY']  = list(set(identity))
                                            metacluster.SV_features['SPECIDENTITY'] = list(set(specificIdentity))
                                            print ('metacluster.SV_featuresmarte2 ' + str(metacluster.SV_features['IDENTITY']) + ' ' + str(metacluster.beg) +' '+ str(event.readName))

                                # If it's another SV_type dont need to check if bkp is near
                                if not SV_TypeTested:
                                    metacluster.SV_features['IDENTITY']  = list(set(identity))
                                    metacluster.SV_features['SPECIDENTITY'] = list(set(specificIdentity))
                                    print ('metacluster.SV_featuresmarte3 ' + str(metacluster.SV_features['IDENTITY']) + ' ' + str(metacluster.beg))

                            # if bwa aligment wasnt successful, check minimap2 alignment
                            elif nombre in allHits_viral_minimap.keys() and nombre in eventsIdentity3.keys():
                                for alig in allHits_viral_minimap[nombre].alignments:
                                    splitedtName = alig.tName.split('|')
                                    identity.append(splitedtName[0])
                                    if len(splitedtName) > 1:
                                        specificIdentity.append(splitedtName[1])

                                # Check if it's near bkp
                                if 'BND' in allMetaclustersIdentity.keys():
                                    print ('BND3')
                                    if metacluster in allMetaclustersIdentity['BND']:
                                        begsEnds = []
                                        SV_TypeTested = True
                                        print ('BND4')
                                        buffer = 251
                                        bkpProximity = virusLR.checkBkpProximity_woSide(event, nombre, allHits_viral_minimap, buffer)

                                        if bkpProximity != None:
                                            metacluster.SV_features['IDENTITY']  = list(set(identity))
                                            metacluster.SV_features['SPECIDENTITY'] = list(set(specificIdentity))
                                            print ('metacluster.SV_featuresmarteMINIMAP1 ' + str(metacluster.SV_features['IDENTITY']) + ' ' + str(metacluster.beg))

                                # Repito con INs_noOverlap, hacer una sola:
                                if 'INS_noOverlap' in allMetaclustersIdentity.keys():
                                    print ('allMetaclustersIdentity2 ' + str(allMetaclustersIdentity['INS_noOverlap']))
                                    if metacluster in allMetaclustersIdentity['INS_noOverlap']:
                                        begsEnds = []
                                        SV_TypeTested = True

                                        buffer = 251
                                        bkpProximity = virusLR.checkBkpProximity_wSide(event, nombre, allHits_viral_minimap, buffer)
                                        event.bkpProximity = bkpProximity

                                        if bkpProximity != None:
                                            ## NUEVO!!!:
                                            #if (event.clippedSide == 'right' and bkpProximity < 250) or (event.clippedSide == 'left' and bkpProximity >= (ali_query_Len - 250)):
                                            metacluster.SV_features['IDENTITY']  = list(set(identity))
                                            metacluster.SV_features['SPECIDENTITY'] = list(set(specificIdentity))
                                            print ('metacluster.SV_featuresmarteMINIMAP2 ' + str(metacluster.SV_features['IDENTITY']) + ' ' + str(metacluster.beg) +' '+ str(event.readName) +' '+ str(event.clippedSide))

                                # If it's another SV_type dont need to check if bkp is near
                                if not SV_TypeTested:
                                    metacluster.SV_features['IDENTITY']  = list(set(identity))
                                    metacluster.SV_features['SPECIDENTITY'] = list(set(specificIdentity))
                                    print ('metacluster.SV_featuresmarte4 ' + str(metacluster.SV_features['IDENTITY']) + ' ' + str(metacluster.beg))
                            
        
        unix.rm([polishDir])
        del allMetaclustersIdentitySeparated
        # Ahora tengo el mismo dictionario que antes, pero marcada la identidad en los eventos que la tienen (en el caso de los clipping, ademas aquellos que hacen clipping cerca del bkp)

        # Split metaclusters dictionary in those that have identity and those with no identity

        allMetaclustersOrganised = {}
        for SV_type, metaclusters in allMetaclustersIdentity.items():
            for metacluster in metaclusters:
                SV_type_identity = SV_type + '_identity'
                SV_type_unknown = SV_type + '_unknown'

                if 'IDENTITY' in metacluster.SV_features.keys():
                    if metacluster.SV_features['IDENTITY']:
                        if SV_type_identity in allMetaclustersOrganised:
                            allMetaclustersOrganised[SV_type_identity].append(metacluster)
                        else:
                            allMetaclustersOrganised[SV_type_identity] = []
                            allMetaclustersOrganised[SV_type_identity].append(metacluster)
                else:
                    if SV_type_unknown in allMetaclustersOrganised:
                        allMetaclustersOrganised[SV_type_unknown].append(metacluster)
                    else:
                        allMetaclustersOrganised[SV_type_unknown] = []
                        allMetaclustersOrganised[SV_type_unknown].append(metacluster)

        del allMetaclustersIdentity

        # allMetaclustersOrganised: Es un dictionario en el que los metaclsuters estan ordenados en funcion de si sus eventos tienen identidad o no.


        ### 3. Determine what type of sequence has been inserted for INS metaclusters
        msg = '3. Determine what type of sequence has been inserted for INS metaclusters'
        log.header(msg)

        # Create output directory
        outDir = self.outDir + '/INS_TYPE/'
        unix.mkdir(outDir)

        if 'INS_identity' in allMetaclustersOrganised:

            ## Infer insertion type
            clusters.INS_type_metaclusters(allMetaclustersOrganised['INS_identity'], self.reference, self.annotations, 1, self.confDict['viralDb'], outDir)

        # Remove output directory
        unix.rm([outDir])
            
        ### 4. Resolve structure for solo, partnered and orphan transductions
        msg = '4. Resolve structure for solo, partnered and orphan transductions'
        log.header(msg)
        
        if 'INS_identity' in allMetaclustersOrganised:
            consensus = self.refDir + '/consensusDb.fa'
            transduced = self.refDir + '/transducedDb.fa.masked'

            # Create output directory
            outDir = self.outDir + '/STRUCTURE/'
            unix.mkdir(outDir)

            # Structure inference
            # TODO: Etso esta solo hecho para RT creo!!
            allMetaclustersOrganised['INS_identity'] = clusters.structure_inference_parallel(allMetaclustersOrganised['INS_identity'], consensus, transduced, self.confDict['transductionSearch'], self.confDict['processes'], outDir)

            # Remove output directory
            unix.rm([outDir])     

        # Resolve identity for INS_noOverlap
        '''
        if 'INS_noOverlap_identity' in allMetaclustersOrganised:
            for metacluster in allMetaclustersOrganised['INS_noOverlap_identity']:
                metacluster.consRightSeq
                metacluster.consLefttSeq
        '''
        ### 5. Identify BND junctions
        msg = '5. Identify BND junctions'
        log.header(msg)
        allJunctions = []

        if 'BND_identity' in allMetaclustersOrganised:

            # 1. Analyse only those BND matching with virus:
            outDir = self.outDir + '/BND_JUNCTIONS/'
            unix.mkdir(outDir)
            
            # From now on, only work with those that have identity
            # TODO: Put as an option to work with all BNDs, as it is done in VIGA_SR
                
            ### Search for repeat, transduction or viral bridges
            # Create output directory
            allMetaclustersOrganised['BND_identity'] = clusters.search4bridges_metaclusters_parallel(allMetaclustersOrganised['BND_identity'], 10000, 80, self.confDict['minReads'], 25, self.annotations, self.refDir, self.confDict['viralDb'], self.confDict['processes'], outDir)

            ### Search for BND junctions
            allJunctions = clusters.search4junctions_metaclusters(allMetaclustersOrganised['BND_identity'], self.refLengths, self.confDict['processes'], self.confDict['minReads'], 25, self.reference, self.refDir, self.confDict['viralDb'], outDir)

            if allJunctions:
                # Analyse BNDjunction structure

                # Initialize variables:
                tupleList = []

                for junction in allJunctions:
                    
                    ## Add to the list
                    fields = (junction, self.reference, self.confDict['viralDb'], outDir)
                    tupleList.append(fields)
                
                l = mp.Lock()
                pool = mp.Pool(processes=self.confDict['processes'])
                allJunctions = pool.starmap(clusters.analyse_BNDjunction_structure, tupleList)
                pool.close()
                pool.join()
                

            # NOTE 2020: New June 2020. For keeping those BNDs without pair
            # TODO: From now on, only work with those that have identity!!!!! Put unknown as option.
            # TODO: Cambiar y poner solo-BND_identity

            for metaclusterBND in allMetaclustersOrganised['BND_identity']:
                if metaclusterBND.beg not in [junction.metaclusterA.beg for junction in allJunctions]:
                    if metaclusterBND.beg not in [junction.metaclusterB.beg for junction in allJunctions]:
                        if 'solo-BND_identity' not in allMetaclustersOrganised:
                            allMetaclustersOrganised['solo-BND_identity'] = []
                        allMetaclustersOrganised['solo-BND_identity'].append(metaclusterBND)
            
            # Make consensus of solo-BND:
            
            if 'solo-BND_identity' in allMetaclustersOrganised.keys():

                
                # Determine soloBND type

                # Initialize variables:
                tupleList = []
                allHits_genome = {}
                groupedEntries = {}

                for metacluster in allMetaclustersOrganised['solo-BND_identity']:

                    ## Add to the list
                    fields = (metacluster, 'BND', self.confDict['processes'], outDir)
                    tupleList.append(fields)
                
                l = mp.Lock()
                pool = mp.Pool(processes=self.confDict['processes'])
                allMetaclustersOrganised['solo-BND_identity'] = pool.starmap(virusLR.get_consensus_nonClassisSVTypes, tupleList)
                pool.close()
                pool.join()

            # DESILENCE
            # Remove output directory
            unix.rm([outDir])
            
        ### 6. Apply second round of filtering 
        msg = '7. Apply second round of filtering'
        log.header(msg)
        filters2Apply = ['PERC-RESOLVED','IDENTITY']
        print ('allMetaclustersOrganised ' + str(allMetaclustersOrganised))
        metaclustersPass, metaclustersFailed = filters.filter_metaclusters(allMetaclustersOrganised, filters2Apply, self.confDict)

        ### 7. Report SV calls into output files
        msg = '8. Report SV calls into output files'
        log.header(msg)

        ## 7.1 Report INS
        if 'INS_identity' in metaclustersPass:
            outFileName = 'INS_MEIGA.PASS'
            virusLR.INS2VCF(metaclustersPass['INS_identity'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        if 'INS_unknown' in metaclustersFailed:
            outFileName = 'INS_MEIGA.FAILED.2'
            virusLR.INS2VCF(metaclustersFailed['INS_unknown'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        ## 7.2 Report BND junctions
        if allJunctions:
            outFileName = 'BND_MEIGA.tsv'
            output.write_junctions(allJunctions, outFileName, self.outDir)

        ## 7.3 Report INS_noOverlap junctions
        if 'INS_noOverlap_identity' in metaclustersPass.keys():
            outFileName = 'INS_noOverlap_MEIGA'
            virusLR.INS2VCF_junction(metaclustersPass['INS_noOverlap_identity'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        if 'INS_noOverlap_unknown' in metaclustersFailed.keys():
            outFileName = 'INS_noOverlap_MEIGA.FAILED'
            virusLR.INS2VCF_junction(metaclustersFailed['INS_noOverlap_unknown'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        ## 7.4 Report solo-BND junctions
        if 'solo-BND_identity' in metaclustersPass.keys():
            outFileName = 'soloBND_MEIGA'
            virusLR.INS2VCF_junction(metaclustersPass['solo-BND_identity'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

        if 'solo-BND_unknown' in metaclustersFailed.keys():
            outFileName = 'soloBND_MEIGA.FAILED'
            virusLR.INS2VCF_junction(metaclustersFailed['solo-BND_unknown'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)
        
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
        # NOTE 2020: New 2020
        #metaclustersPass = pool.starmap(self.make_clusters_bin, bins)
        pool.close()
        pool.join()

        # Remove output directory
        unix.rm([self.outDir + '/CLUSTER/'])

        ### 3. Collapse metaclusters in a single dict and report metaclusters that failed filtering
        # NOTE 2020: In 2020:
        #metaclustersPass = structures.merge_dictionaries(metaclustersPass)

        metaclustersPass = structures.merge_dictionaries(metaclustersPassList)
        metaclustersFailed = structures.merge_dictionaries(metaclustersFailedList)
        
        if 'INS' in metaclustersFailed:
            outFileName = 'INS_MEIGA.FAILED.1.tsv'
            output.write_INS(metaclustersFailed['INS'], outFileName, self.outDir)

        if 'INS_noOverlap' in metaclustersFailed:
            outFileName = 'INS_noOverlap_MEIGA.FAILED.1.tsv'
            output.write_INS(metaclustersFailed['INS_noOverlap'], outFileName, self.outDir)

        # Output FAILED BND:
        if 'BND' in metaclustersFailed:
            outFileName = 'BND_MEIGA_FAILED_1'
            virusLR.INS2VCF_junction(metaclustersFailed['BND'], self.minimap2_index(), self.refLengths, self.confDict['source'], self.confDict['build'], self.confDict['species'], outFileName, self.outDir)

                
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
            eventsDict = bamtools.collectSV(ref, beg, end, self.bam, self.confDict, None)

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
        # Get consensus event and FASTA for INS and DEL.
        # Get bkpPos and supplClusters for BND.
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
        # NOTE 2020: New 2020
        metaclustersSVType, metaclustersSVTypeFailed = filters.filter_metaclusters(metaclustersSVType, filters2Apply, self.confDict)
        #metaclustersSVType = filters.filter_metaclusters(metaclustersSVType, filters2Apply, self.confDict)[0]

        ## 8. Generate consensus event for SV metaclusters ##
        step = 'CONSENSUS'
        msg = 'Generate consensus event for SV metaclusters' 
        log.step(step, msg)

        targetSV = ['INS']
        outDir = binDir + '/CONSENSUS/' 
        clusters.create_consensus(metaclustersSVType, self.confDict, self.reference, targetSV, outDir)
        # TODO: Create consensus also for metaclustersSVTypeFailed

        # Consensus for INS_noOverlap
        for SV, metaclusters in metaclustersSVType.items():
            if SV == 'INS_noOverlap':
                for metacluster in metaclusters:
                    virusLR.get_consensus_nonClassisSVTypes(metacluster, 'INS_noOverlap', self.confDict['processes'], outDir)

        for SV, metaclusters in metaclustersSVTypeFailed.items():
            if SV == 'INS_noOverlap':
                for metacluster in metaclusters:
                    virusLR.get_consensus_nonClassisSVTypes(metacluster, 'INS_noOverlap', self.confDict['processes'], outDir)

        ## 9. Lighten up metaclusters  ##
        ## Metaclusters passing all the filters
        clusters.lighten_up_metaclusters(metaclustersSVType)
        
        ## Filtered metaclusters
        clusters.lighten_up_metaclusters(metaclustersSVTypeFailed)

        # Do cleanup
        #([outDir, binDir])

        ## Print time taken to process bin
        end = time.time()
        time_taken = end - start
        msg = 'SV calling in bin: ' + binId + ' finished in ' + str(time_taken)
        log.info(msg)

        return metaclustersSVType, metaclustersSVTypeFailed
        #return metaclustersSVType