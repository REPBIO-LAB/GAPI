#!/usr/bin/env python
#coding: utf-8

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
	parser.add_argument('--normalBam', default=None, dest='normalBam', help='Matched normal bam file. If provided MEIGA will run in PAIRED mode')
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
	parser.add_argument('--minMAPQ', default=10, dest='minMAPQ', type=int, help='Minimum mapping quality required for each read. Default: 10')
	parser.add_argument('--readFilters', default="SMS", dest='readFilters', type=str, help='Comma separated list of read filters to apply (SMS)')
	parser.add_argument('--readOverhang', default=5000, dest='overhang', type=int, help='Number of flanking base pairs around the SV event to be collected from the supporting read sequence. Default: 5000')
	parser.add_argument('--minINDELlen', default=50, dest='minINDELlen', type=int, help='Minimum indel length. Default: 50')
	parser.add_argument('--minCLIPPINGlen', default=500, dest='minCLIPPINGlen', type=int, help='Minimum clipped sequence length for each read. Default: 500')

	## Clustering
	parser.add_argument('--INSdist', default=250, dest='maxInsDist', type=int, help='Maximum distance bewteen two adjacent INS to be clustered together (Between 0-999). Default: 250')
	parser.add_argument('--BKPdist', default=50, dest='maxBkpDist', type=int, help='Maximum distance bewteen two adjacent breakpoints for CLIPPING clustering (Between 0-999). Default: 250')
	parser.add_argument('--minPercOverlap', default=70, dest='minPercRcplOverlap', type=int, help='Minimum percentage of reciprocal overlap for DEL clustering. Default: 70')

	## Filtering
	# Long
	parser.add_argument('--minClusterSize', default=2, dest='minClusterSize', type=int, help='Minimum number of reads composing a cluster. Default: 2')
	parser.add_argument('--maxClusterSize', default=500, dest='maxClusterSize', type=int, help='Maximum number of reads composing a metacluster. Default: 500')
	parser.add_argument('--maxClusterCV', default=30, dest='maxClusterCV', type=int, help='Maximum coefficient of variation of a metacluster. Default: 30')
	parser.add_argument('--minSupportingReads', default=3, dest='minSupportingReads', type=int, help='Minimum number of reads supporting a SV. Default: 3')
	parser.add_argument('--minNormalSupportingReads', default=2, dest='minNormalSupportingReads', type=int, help='Minimum number of reads supporting a SV in normal sample. Default: 2')
	parser.add_argument('--targetStatus', default='resolved,partially_resolved', dest='targetStatus', type=str, help='Filter out those insertions with an status not included in the list. Default: resolved,partially_resolved')

	# Short
	parser.add_argument('--minReadsRegionMQ', default=10, dest='minReadsRegionMQ', type=int, help='Surrounding reads above this MQ are considered low MQ reads. Default: 10')
	parser.add_argument('--maxRegionlowMQ', default=0.3, dest='maxRegionlowMQ', type=int, help='Maximum percentage of lowMAPQ/nbReads in cluster´s region. Default: 0.3')
	parser.add_argument('--maxRegionSMS', default=0.15, dest='maxRegionSMS', type=int, help='Maximum percentage of SMS clipping reads in cluster´s region. Default: 0.15')

	## 2. Parse user´s input and initialize variables ##
	args = parser.parse_args()

	### Mandatory arguments
	bam = args.bam
	technology = args.technology
	reference = args.reference
	refDir = args.refDir

	### Optional arguments
	## General
	normalBam = args.normalBam
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

	## Filtering thresholds
	# Long
	minClusterSize = args.minClusterSize
	maxClusterSize = args.maxClusterSize
	maxClusterCV = args.maxClusterCV
	minSupportingReads = args.minSupportingReads
	minNormalSupportingReads = args.minNormalSupportingReads
	targetStatus = args.targetStatus

	# Short
	minReadsRegionMQ = args.minReadsRegionMQ
	maxRegionlowMQ = args.maxRegionlowMQ
	maxRegionSMS = args.maxRegionSMS

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
	version='0.13.4'

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
	print('normalBam: ', normalBam)
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
	print('minCLIPPINGlength: ', minCLIPPINGlen, "\n")

	print('** Clustering **')
	print('maxInsDist: ', maxInsDist)
	print('maxBkpDist: ', maxBkpDist)
	print('minPercOverlap: ', minPercRcplOverlap, "\n")

	print('** Filtering **')
	print('minClusterSize: ', minClusterSize)
	print('maxClusterSize: ', maxClusterSize)
	print('maxClusterCV: ', maxClusterCV)
	print('minSupportingReads: ', minSupportingReads)
	print('minNormalSupportingReads: ', minNormalSupportingReads)
	print('targetStatus: ', targetStatus)
	print('minReadsRegionMQ: ', minReadsRegionMQ)
	print('maxRegionlowMQ: ', maxRegionlowMQ)
	print('maxRegionSMS: ', maxRegionSMS, "\n")
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

	# a) Illumina data (WGS or sureselect)
	if (confDict['technology'] == 'ILLUMINA') or (confDict['technology'] == 'SURESELECT'):
		confDict['targetEvents'] = 'DISCORDANT'
	
	# b) Long read sequencing data -> INS or INS + BND
	elif ('INS' in confDict['targetSV']):
		confDict['targetEvents'] = 'INS,CLIPPING'

	# c) BND alone
	elif ('BND' in confDict['targetSV']):
		confDict['targetEvents'] = 'CLIPPING'

	## Clustering
	confDict['maxInsDist'] = maxInsDist
	confDict['maxBkpDist'] = maxBkpDist
	confDict['minPercRcplOverlap'] = minPercRcplOverlap

	## Filtering thresholds
	# Long
	confDict['minClusterSize'] = minClusterSize
	confDict['maxClusterSize'] = maxClusterSize
	confDict['maxClusterCV'] = maxClusterCV
	confDict['minSupportingReads'] = minSupportingReads
	confDict['minNormalSupportingReads'] = minNormalSupportingReads
	confDict['targetStatus'] = targetStatus.split(',')
	confDict['minPercResolved'] = 40

	## Filtering thresholds short reads
	confDict['minReadsRegionMQ'] = minReadsRegionMQ
	confDict['maxRegionlowMQ'] = maxRegionlowMQ
	confDict['maxRegionSMS'] = maxRegionSMS
	
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
