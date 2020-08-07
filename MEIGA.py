#!/usr/bin/env python
#coding: utf-8

####################
## Import modules ##
####################

# import basic internal module
import check_dependencies as cd

if __name__ == '__main__':
	
	## Check dependencies, true if there are some missing dependencies ##
	missing_dependencies=cd.missing_python_dependencies() or cd.missing_program_dependencies()

	####################
	## Import modules ##
	####################

	# External
	import argparse
	import sys
	import os
	import multiprocessing as mp

	#internal modules import delayed until dependency check. ()
	
	mp.set_start_method('spawn')
    
	######################
	## Get user's input ##
	######################

	## 1. Define parser ##
	### Mandatory arguments
	parser = argparse.ArgumentParser(description='Call mobile element insertions (MEI) and viral integrations from second and third generation sequencing data. Two running modes: 1) SINGLE: individual sample; 2) PAIRED: tumour and matched normal sample')
	parser.add_argument('bam', help='Input bam file. Will correspond to the tumour sample in the PAIRED mode')
	parser.add_argument('technology', help='Sequencing technology used to generate the data (NANOPORE, PACBIO, ILLUMINA or SURESELECT)')
	parser.add_argument('reference', help='Reference genome in fasta format. An index of the reference generated with samtools faidx must be located in the same directory')
	parser.add_argument('refDir', help='Directory containing reference databases (consensus sequences, source elements...)')

	### Optional arguments
	## General
	parser.add_argument('--species', default='Homo sapiens', dest='species', help='Target species. Default: Homo sapiens')
	parser.add_argument('--build', default='GRCh37', dest='build', help='Reference genome build. Default: GRCh37')
	parser.add_argument('--normalBam', default=None, dest='normalBam', help='Matched normal bam file. If provided MEIGA will run in PAIRED mode')
	parser.add_argument('--germlineMEI', default=None, dest='germlineMEI', help='Bed file containing set of known germline MEI')
	parser.add_argument('--transduction-search', action="store_true", default=False, dest='transductionSearch', help='Enable transduction search. If not enabled only solo events will be identified')
	parser.add_argument('--source-families', default=None, dest='srcFamilies', type=str, help='Comma separated list of possible families for source elements mediating transductions. Default: None. Mandatory if transduction search enabled')
	parser.add_argument('--gene-annot-dir', default=None, dest='annovarDir', help='Directory containing annovar reference files for gene-based annotation of MEI breakpoints. If not provided gene annotation step will be skipped')
	parser.add_argument('--polishing-rounds', default=1, dest='rounds', type=int, help='Number of polishing rounds to be attempted. Default: 1')
	parser.add_argument('-p', '--processes', default=1, dest='processes', type=int, help='Number of processes. Default: 1')
	parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='Output directory. Default: current working directory')

	## BAM processing
	parser.add_argument('--targetBins', default=None, dest='targetBins', type=str, help='Bed file containing target genomic bins for SV calling. Overrides --binSize and --refs. Default: None')
	parser.add_argument('-bS', '--binSize', default=1000000, dest='binSize', type=int, help='Input bams will be analised in genomic bins of this size. Default: 1000000')
	parser.add_argument('--refs', default="ALL", dest='refs', type=str, help='Comma separated list of target references to call SV (i.e. 1,2,3,X). Default: All references included in the bam file')
	parser.add_argument('--target-SV', default="INS", dest='targetSV', type=str, help='Comma separated list of SV types to be detected (INS AND/OR BND). Default: INS')
	parser.add_argument('--no-duplicates', action="store_true", default=False, dest='filterDuplicates', help='Filter out reads marked as duplicates if filter enabled')
	parser.add_argument('--minMAPQ', default=20, dest='minMAPQ', type=int, help='Minimum mapping quality required for each read. Default: 20')
	parser.add_argument('--readFilters', default="SMS", dest='readFilters', type=str, help='Comma separated list of read filters to apply (SMS)')
	parser.add_argument('--readOverhang', default=5000, dest='overhang', type=int, help='Number of flanking base pairs around the SV event to be collected from the supporting read sequence. Default: 5000')
	parser.add_argument('--minINDELlen', default=50, dest='minINDELlen', type=int, help='Minimum indel length. Default: 50')
	parser.add_argument('--minCLIPPINGlen', default=None, dest='minCLIPPINGlen', type=int, help='Minimum clipped sequence length for each read. Default [ILLUMINA, SURESELECT]: 20; [NANOPORE, PACBIO]: 500')
	
	## Clustering
	parser.add_argument('--INSdist', default=250, dest='maxInsDist', type=int, help='Maximum distance bewteen two adjacent INS to be clustered together (Between 0-999). Default: 250')
	parser.add_argument('--BKPdist', default=50, dest='maxBkpDist', type=int, help='Maximum distance bewteen two adjacent breakpoints for CLIPPING clustering (Between 0-999). Default: 250')
	parser.add_argument('--minPercOverlap', default=70, dest='minPercRcplOverlap', type=int, help='Minimum percentage of reciprocal overlap for DEL clustering. Default: 70')
	parser.add_argument('--equalOrientBuffer', default=200, dest='equalOrientBuffer', type=int, help='Distance between reads that are equally oriented. Default: 200')
	parser.add_argument('--oppositeOrientBuffer', default=600, dest='oppositeOrientBuffer', type=int, help='Distance between reads that are opposite oriented. Default: 600')
	parser.add_argument('--libraryReadLength', default=151, dest='libraryReadLength', type=int, help='Illumina library read length. Default: 151')

	## Databases
	parser.add_argument('--viralDb', default=None, dest='viralDb', type=str, help='Viral database in fasta format or minimap index.')

	## Filtering
	parser.add_argument('--minClusterSize', default=4, dest='minClusterSize', type=int, help='Minimum number of reads composing a cluster. Default: 4')
	parser.add_argument('--maxClusterSize', default=500, dest='maxClusterSize', type=int, help='Maximum number of reads composing a metacluster. Default: 500')
	parser.add_argument('--minReads', default=4, dest='minReads', type=int, help='Minimum number of reads supporting a SV. Default: 4')
	parser.add_argument('--minNormalReads', default=2, dest='minNormalReads', type=int, help='Minimum number of reads supporting a SV in normal sample. Default: 2')

	# Long
	parser.add_argument('--maxClusterCV', default=40, dest='maxClusterCV', type=int, help='Maximum coefficient of variation of a metacluster. Default: 40')
	parser.add_argument('--targetStatus', default='resolved,partially_resolved', dest='targetStatus', type=str, help='Filter out those insertions with an status not included in the list. Default: resolved,partially_resolved')

	# Short
	parser.add_argument('--minReadsRegionMQ', default=10, dest='minReadsRegionMQ', type=int, help='Surrounding reads above this MQ are considered low MQ reads. Default: 10')
	parser.add_argument('--maxRegionlowMQ', default=0.3, dest='maxRegionlowMQ', type=int, help='Maximum percentage of lowMAPQ/nbReads in cluster´s region. Default: 0.3')
	parser.add_argument('--maxRegionSMS', default=0.15, dest='maxRegionSMS', type=int, help='Maximum percentage of SMS clipping reads in cluster´s region. Default: 0.15')
	parser.add_argument('--INT2Search', default="ME,VIRUS", dest='INT2Search', type=str, help='Comma separated list of insertion types to collect (Mobile Elements (ME),VIRUS). Default: ME,VIRUS')
	parser.add_argument('--komplexityThreshold', default=0.4, dest='komplexityThreshold', type=float, help='Threshold for filtering mates sequence with komplexity tool. Default: 0.4')
	parser.add_argument('--discordantMatesMaxMAPQ', default=20, dest='discordantMatesMaxMAPQ', type=int, help='Maximum mapping quality used for collecting dicordant read mates. Default: 20')
	parser.add_argument('--no-discordantMatesCheckUnmapped', action="store_false", default=True, dest='discordantMatesCheckUnmapped', help='If not selected, when a dicordant read mate is unmapped, collect it no matter its MAPQ. If selected, mapping state is not checked.')
	parser.add_argument('--no-discordantMatesSupplementary', action="store_false", default=True, dest='discordantMatesSupplementary', help='When selected, avoid collecting dicordant read mates that are supplementary alignments.')
	parser.add_argument('--discordantMatesMaxBasePerc', default=85, dest='discordantMatesMaxBasePerc', type=int, help='Maximum base percentage of discordant read mates sequences. Default: 85')
	parser.add_argument('--discordantMatesMinLcc', default=1.49, dest='discordantMatesMinLcc', type=float, help='Minimum local complexity of discordant read mates sequences. Default: 1.49')
	parser.add_argument('--MEDB', default=None, dest='MEDB', type=str, help='Path to ME database used for analysing metacluster bkp.')
	parser.add_argument('--filtersBfClip', default="MAX-NBREADS,AREAMAPQ,AREASMS,IDENTITY", dest='filtersBfClip', type=str, help='Comma separated list of filters to apply before adding clippings to metaclusters. Default: MAX-NBREADS,AREAMAPQ,AREASMS,IDENTITY')
	parser.add_argument('--filtersAfClip', default="MIN-NBREADS,MAX-NBREADS", dest='filtersAfClip', type=str, help='Comma separated list of filters to apply after adding clippings to metaclusters. Default: MIN-NBREADS, MAX-NBREADS')
	parser.add_argument('--analyseFiltered', action="store_true", default=False, dest='analyseFiltered', help='If selected, add clippings and analyse breakpoint of those metaclusters that do not PASS selected filters.')
	
 	# RetroTest
	parser.add_argument('--blatClip', action="store_true", default=False, dest='blatClip', help='When selected, blat realignment will be performed with clippings. Higher sensitivity, but time-consuming. Default: False')
	parser.add_argument('--retroTestWGS', action="store_true", default=False, dest='retroTestWGS', help='Apply Retrotest method on WGS data. Default: False')
	
	# Filtering viral bam
	parser.add_argument('--minTotalMatchVirus', default=40, dest='minTotalMatchVirus', type=int, help='Minimum total matches of a read against viral DB. Default: 40.')
	parser.add_argument('--minParcialMatchVirus', default=15, dest='minParcialMatchVirus', type=int, help='Minimum length of a match against viral DB. Default: 15.')
	parser.add_argument('--maxMatchCheckMAPQVirus', default=60, dest='maxMatchCheckMAPQVirus', type=int, help='MAPQ of viral DB matches above this threshold is not interrogated. Default: 60.')
	parser.add_argument('--minMAPQVirus', default=0, dest='minMAPQVirus', type=int, help='MAPQ of viral DB matches has to be above this threshold. Default: 0.')
	parser.add_argument('--maxBasePercVirus', default=85, dest='maxBasePercVirus', type=int, help='Maximum base percentage of a match against viralDB. Default: 85.')
	parser.add_argument('--minLccVirus', default=1.49, dest='minLccVirus', type=float, help='Minimum local complexity (lcc) of a match against viralDB. Default: 1.49.')

	# Output
	parser.add_argument('--VCFInfoFields', default="VTYPE,NBTOTAL,NBTUMOR,NBNORMAL,LEN,NBDISC,NBCLIP,IDENT,ORIENT,BKP2,DISCMAPQ,CLIPMAPQ,SPECIDENT,DISCDUP,CLIPDUP,REP,REPSUB,DIST,REGION,GENE,INTBKP,INTBKP2,CLIPTYPE,REFSeq,BKPCSEQ,BKP2CSEQ,BKPRSEQ,BKP2RSEQ,CLIPTYPE2,DISC,CLIP", dest='VCFInfoFields', type=str, help='Comma separated list of INFO fields to display in output VCF (VTYPE,NBTOTAL,NBTUMOR,NBNORMAL,LEN,NBDISC,NBCLIP,IDENT,ORIENT,BKP2,DISCMAPQ,CLIPMAPQ,SPECIDENT,DISCDUP,CLIPDUP,REP,REPSUB,DIST,REGION,GENE,INTBKP,INTBKP2,CLIPTYPE,REFSeq,BKPCSEQ,BKP2CSEQ,BKPRSEQ,BKP2RSEQ,CLIPTYPE2,DISC,CLIP) *NOTE that REFSeq is incompatible with --no-VCFREF. Default: All are included')
	parser.add_argument('--no-annotRepeats', action="store_false", default=True, dest='annotRepeats', help='If selected, not show annotated repeats on the reference genome at insertion interval. Only works with VIRUS mode. If ME are analysed, annotation repeats step is always performed.')
	parser.add_argument('--no-VCFREF', action="store_false", default=True, dest='VCFREF', help='If selected, not show REF field in VCF output file (it consumes ~5Gb).')
	parser.add_argument('--no-consensusBkpSeq', action="store_false", default=True, dest='consBkpSeq', help='If selected, a representative read is selected for breakpoint sequence. Otherwise, make consensus of breakpoint sequence. Time and space saver.')
	parser.add_argument('--keepIdentDb', action="store_true", default=False, dest='keepIdentDb', help='If selected, do not delete the fasta file that is created with detected identities in input sample.')
	parser.add_argument('--no-FAILED-VCF', action="store_false", default=True, dest='printFiltered', help='If selected, not print VCF with those results than do not PASS the filters.')

	## 2. Parse user´s input and initialize variables ##
	args = parser.parse_args()

	### Mandatory arguments
	bam = args.bam
	technology = args.technology
	reference = args.reference
	refDir = args.refDir

	### Optional arguments
	## General
	species = args.species
	build = args.build
	normalBam = args.normalBam
	germlineMEI = args.germlineMEI
	transductionSearch = args.transductionSearch
	srcFamilies = args.srcFamilies
	annovarDir = args.annovarDir
	rounds = args.rounds
	processes = args.processes
	outDir = args.outDir

	## BAM processing
	targetBins = args.targetBins
	binSize = args.binSize
	refs = args.refs
	targetSV = args.targetSV
	filterDuplicates = args.filterDuplicates
	minMAPQ = args.minMAPQ
	readFilters = args.readFilters
	overhang = args.overhang
	minINDELlen = args.minINDELlen
	minCLIPPINGlen = args.minCLIPPINGlen

	## Clustering
	maxInsDist = args.maxInsDist
	maxBkpDist = args.maxBkpDist
	minPercRcplOverlap = args.minPercRcplOverlap
	equalOrientBuffer = args.equalOrientBuffer
	oppositeOrientBuffer = args.oppositeOrientBuffer
	libraryReadLength = args.libraryReadLength
	
	## Databases
	viralDb = args.viralDb

	## Filtering thresholds
	# Long
	minClusterSize = args.minClusterSize
	maxClusterSize = args.maxClusterSize
	maxClusterCV = args.maxClusterCV
	minReads = args.minReads
	minNormalReads = args.minNormalReads
	targetStatus = args.targetStatus

	# Short
	minReadsRegionMQ = args.minReadsRegionMQ
	maxRegionlowMQ = args.maxRegionlowMQ
	maxRegionSMS = args.maxRegionSMS
	INT2Search = args.INT2Search
	komplexityThreshold = args.komplexityThreshold
	discordantMatesMaxMAPQ = args.discordantMatesMaxMAPQ
	discordantMatesCheckUnmapped = args.discordantMatesCheckUnmapped
	discordantMatesSupplementary = args.discordantMatesSupplementary
	discordantMatesMaxBasePerc = args.discordantMatesMaxBasePerc
	discordantMatesMinLcc = args.discordantMatesMinLcc
	MEDB = args.MEDB
	filtersBfClip = args.filtersBfClip
	filtersAfClip = args.filtersAfClip
	analyseFiltered = args.analyseFiltered

	minTotalMatchVirus = args.minTotalMatchVirus
	minParcialMatchVirus = args.minParcialMatchVirus
	maxMatchCheckMAPQVirus = args.maxMatchCheckMAPQVirus
	minMAPQVirus = args.minMAPQVirus
	maxBasePercVirus = args.maxBasePercVirus
	minLccVirus = args.minLccVirus

	# RetroTest
	blatClip = args.blatClip
	retroTestWGS = args.retroTestWGS
	
	# Ouput
	VCFInfoFields = args.VCFInfoFields
	annotRepeats = args.annotRepeats
	VCFREF = args.VCFREF
	consBkpSeq = args.consBkpSeq
	keepIdentDb = args.keepIdentDb
	printFiltered = args.printFiltered
	
    ## Check file dependencies (list with the paths) ##
	missing_dependencies = missing_dependencies or  cd.missing_needed_files((bam,reference,refDir,normalBam,annovarDir,outDir))

	if missing_dependencies:
		exit()

	###################################################################
	# lazy load of the modules to avoid problems with the dependences #
	###################################################################
	
	# load Internal
	import callers
	import log
	import bamtools

    ## Check file dependencies (list with the paths) ##
	missing_dependencies = missing_dependencies or  cd.missing_needed_files((bam,reference,refDir,normalBam,annovarDir,outDir))

	if missing_dependencies:
		exit()

	###################################################################
	# lazy load of the modules to avoid problems with the dependences #
	###################################################################
	
	# load Internal
	import callers
	import log
	import bamtools

	# If no reference is specified, get all that are present in the bam file.
	if refs == 'ALL':
		refs = bamtools.get_refs(bam)

	# Convert comma-separated string inputs into lists:
	targetSV = targetSV.split(',')
	targetRefs = refs.split(',')
	targetINT2Search = INT2Search.split(',')
	targetVCFInfoFields = VCFInfoFields.split(',')
	filtersBfClipList = filtersBfClip.split(',')
	filtersAfClipList = filtersAfClip.split(',')

	## Determine running mode:
	mode = 'SINGLE' if normalBam == None else 'PAIRED'

	## If unknown technology provided raise an error and exit 
	if technology not in ['NANOPORE', 'PACBIO', 'ILLUMINA', 'SURESELECT']:
		log.info('[ERROR] Abort execution as ' + technology + ' technology not supported')
		sys.exit(1)

	## Abort if SV unknown SV type provided

	for SV_type in targetSV:		
		if SV_type not in ['INS', 'BND']:
			log.info('[ERROR] Abort execution as ' + SV_type + ' SV type not supported')
			sys.exit(1)
			

	## Abort if transduction search enabled and target families for source elements not provided
	if (transductionSearch) and (srcFamilies is None):
		log.info('[ERROR] Abort execution as transduction search enabled (--transduction-search) and target families for source elements not provided (--source-families)')
		sys.exit(1)	

	##############################################
	## Display configuration to standard output ##
	##############################################
	scriptName = os.path.basename(sys.argv[0])
	scriptName = os.path.splitext(scriptName)[0]
	version='0.24.0'

	print()
	print('***** ', scriptName, version, 'configuration *****')
	print('*** Mandatory arguments ***')
	print('bam: ', bam)
	print('technology: ', technology)
	print('reference: ', reference)
	print('refDir: ', refDir, "\n")

	print('*** Optional arguments ***')
	print('** General **')
	print('mode: ', mode)
	print('species: ', species)
	print('build: ', build)
	print('normalBam: ', normalBam)
	print('germlineMEI: ', germlineMEI)
	print('transduction-search: ', transductionSearch)
	print('source-families: ', srcFamilies)
	print('gene-annot-dir: ', annovarDir)
	print('polishing-rounds: ', rounds)
	print('processes: ', processes)
	print('outDir: ', outDir, "\n")

	print('** BAM processing **')
	print('targetBins: ', targetBins)
	print('binSize: ', binSize)
	print('targetRefs: ', refs)
	print('targetSV: ', targetSV)
	print('filterDuplicates: ', filterDuplicates)
	print('minMAPQ: ', minMAPQ)
	print('readFilters: ', readFilters)
	print('overhang: ', overhang)
	print('minINDELlength: ', minINDELlen)
	print('minCLIPPINGlength: ', minCLIPPINGlen)
	print('targetINT2Search ', INT2Search, "\n")

	print('** Clustering **')
	print('maxInsDist: ', maxInsDist)
	print('maxBkpDist: ', maxBkpDist)
	print('minPercOverlap: ', minPercRcplOverlap)
	print('equalOrientBuffer: ', equalOrientBuffer)
	print('oppositeOrientBuffer: ', oppositeOrientBuffer)
	print('libraryReadLength: ', libraryReadLength, "\n")

	print('** Databases **')
	print('viralDb: ', viralDb, "\n")

	print('** Filtering **')
	print('minClusterSize: ', minClusterSize)
	print('maxClusterSize: ', maxClusterSize)
	print('maxClusterCV: ', maxClusterCV)
	print('minReads: ', minReads)
	print('minNormalReads: ', minNormalReads)
	print('targetStatus: ', targetStatus)
	print('minReadsRegionMQ: ', minReadsRegionMQ)
	print('maxRegionlowMQ: ', maxRegionlowMQ)
	print('maxRegionSMS: ', maxRegionSMS)
	print('komplexityThreshold: ', komplexityThreshold)
	print('discordantMatesMaxMAPQ: ', discordantMatesMaxMAPQ)
	print('discordantMatesCheckUnmapped: ', discordantMatesCheckUnmapped)
	print('discordantMatesSupplementary: ', discordantMatesSupplementary)
	print('discordantMatesMaxBasePerc: ', discordantMatesMaxBasePerc)
	print('discordantMatesMinLcc: ', discordantMatesMinLcc)
	print('MEDB: ', MEDB)
	print('minTotalMatchVirus', minTotalMatchVirus)
	print('minParcialMatchVirus', minParcialMatchVirus)
	print('maxMatchCheckMAPQVirus', maxMatchCheckMAPQVirus)
	print('minMAPQVirus', minMAPQVirus)
	print('maxBasePercVirus', maxBasePercVirus)
	print('minLccVirus', minLccVirus)
	print('filtersBfClip', filtersBfClip)
	print('filtersAfClip', filtersAfClip)
	print('analyseFiltered', analyseFiltered, "\n")

	print('** Output format**')
	print ('VCFInfoFields: ', VCFInfoFields)
	print ('annotRepeats: ', annotRepeats)
	print ('VCFREF: ', VCFREF)
	print ('consBkpSeq: ', consBkpSeq)
	print ('keepIdentDb: ', keepIdentDb)
	print ('printFiltered: ', printFiltered, "\n")

	print('***** Executing ', scriptName, '.... *****', "\n")

	##########
	## CORE ## 
	##########

	## 1. Create configuration dictionary
	#######################################
	confDict = {}

	## Mandatory
	confDict['technology'] = technology

	## General
	confDict['source'] = 'MEIGA-' + version
	confDict['species'] = species
	confDict['build'] = build
	confDict['germlineMEI'] = germlineMEI
	confDict['transductionSearch'] = transductionSearch
	confDict['srcFamilies'] = srcFamilies.split(',') if srcFamilies is not None else []
	confDict['annovarDir'] = annovarDir
	confDict['rounds'] = rounds
	confDict['processes'] = processes

	## BAM processing
	confDict['targetBins'] = targetBins
	confDict['binSize'] = binSize
	confDict['targetRefs'] = targetRefs
	confDict['filterDuplicates'] = filterDuplicates	
	confDict['readFilters'] = readFilters
	confDict['overhang'] = overhang
	confDict['minMAPQ'] = minMAPQ
	confDict['minINDELlen'] = minINDELlen
	confDict['minCLIPPINGlen'] = minCLIPPINGlen

	## Target SV events to search for
	confDict['targetSV'] = targetSV

	## Technology specific parameters
	# a) WGS illumina data (WGS)
	if (confDict['technology'] == 'ILLUMINA'):
		confDict['targetEvents'] = ['DISCORDANT']
		if confDict['minCLIPPINGlen'] == None:
			confDict['minCLIPPINGlen'] = 20

	# b) Sureselect illumina data
	elif (confDict['technology'] == 'SURESELECT'):
		confDict['blatClip'] = blatClip
		confDict['retroTestWGS'] = retroTestWGS
		confDict['targetEvents'] = ['DISCORDANT', 'CLIPPING']
		confDict['minNbDISCORDANT'] = minClusterSize
		confDict['minNbCLIPPING'] = minClusterSize
		if confDict['minCLIPPINGlen'] == None:
			confDict['minCLIPPINGlen'] = 20

	# c) Long read sequencing data -> INS or INS + BND
	elif ('INS' in confDict['targetSV']):
		confDict['targetEvents'] = ['INS', 'CLIPPING']
		if confDict['minCLIPPINGlen'] == None:
			confDict['minCLIPPINGlen'] = 500

	# d) BND alone
	elif ('BND' in confDict['targetSV']):
		confDict['targetEvents'] = ['CLIPPING']
		if confDict['minCLIPPINGlen'] == None:
			confDict['minCLIPPINGlen'] = 500
	
	## Clustering
	confDict['maxInsDist'] = maxInsDist
	confDict['maxBkpDist'] = maxBkpDist
	confDict['minPercRcplOverlap'] = minPercRcplOverlap
	confDict['equalOrientBuffer'] = equalOrientBuffer
	confDict['oppositeOrientBuffer'] = oppositeOrientBuffer
	confDict['libraryReadLength'] = libraryReadLength
		
	## Databases
	confDict['viralDb'] = viralDb

	## Filtering thresholds
	# Long
	confDict['minClusterSize'] = minClusterSize
	confDict['maxClusterSize'] = maxClusterSize
	confDict['maxClusterCV'] = maxClusterCV
	confDict['minReads'] = minReads
	confDict['minNormalReads'] = minNormalReads
	confDict['targetStatus'] = targetStatus.split(',')
	# TODO: maybe add special for viruses
	confDict['minPercResolved'] = 40

	## Filtering thresholds short reads
	confDict['minReadsRegionMQ'] = minReadsRegionMQ
	confDict['maxRegionlowMQ'] = maxRegionlowMQ
	confDict['maxRegionSMS'] = maxRegionSMS
	confDict['targetINT2Search'] = targetINT2Search
	confDict['komplexityThreshold'] = komplexityThreshold
	confDict['discordantMatesMaxMAPQ'] = discordantMatesMaxMAPQ
	confDict['discordantMatesCheckUnmapped'] = discordantMatesCheckUnmapped
	confDict['discordantMatesSupplementary'] = discordantMatesSupplementary
	confDict['discordantMatesMaxBasePerc'] = discordantMatesMaxBasePerc
	confDict['discordantMatesMinLcc'] = discordantMatesMinLcc
	confDict['MEDB'] = MEDB
	confDict['filtersBfClip'] = filtersBfClipList
	confDict['filtersAfClip'] = filtersAfClipList
	confDict['analyseFiltered'] = analyseFiltered
	confDict['minTotalMatchVirus'] = minTotalMatchVirus
	confDict['minParcialMatchVirus'] = minParcialMatchVirus
	confDict['maxMatchCheckMAPQVirus'] = maxMatchCheckMAPQVirus
	confDict['minMAPQVirus'] = minMAPQVirus
	confDict['maxBasePercVirus'] = maxBasePercVirus
	confDict['minLccVirus'] = minLccVirus

	# Output
	confDict['VCFInfoFields'] = targetVCFInfoFields
	confDict['annotRepeats'] = annotRepeats
	confDict['VCFREF'] = VCFREF
	confDict['consBkpSeq'] = consBkpSeq
	confDict['keepIdentDb'] = keepIdentDb
	confDict['printFiltered'] = printFiltered
	
	## 2. Execute structural variation caller
	###########################################

	### Create caller
	# A) Pacbio or Nanopore long reads 
	if confDict['technology'] in ['NANOPORE', 'PACBIO']:
		caller = callers.SV_caller_long(mode, bam, normalBam, reference, refDir, confDict, outDir)

	# B) Illumina short reads
	elif confDict['technology'] == 'ILLUMINA':
		caller = callers.SV_caller_short(mode, bam, normalBam, reference, refDir, confDict, outDir)

	# C) Source elements sureselect data
	elif confDict['technology'] == 'SURESELECT':
		caller = callers.SV_caller_sureselect(mode, bam, normalBam, reference, refDir, confDict, outDir)

	### Do calling
	caller.call()

	print('***** Finished! *****')
	print()
