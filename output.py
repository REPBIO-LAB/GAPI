## DEPENDENCIES ##
# External
import pybedtools
import mappy as mp
import itertools
import statistics

# Internal
import structures
import formats

def VCFMetaclustersFields(metaclusters):
    '''
    Get from metaclusters those fields that should be printed in VCF output

    Input:
        1. metaclusters: list containing list of metaclusters

    Output:
        1. metaclustersFields: List of lists: containing VCF fields for each metacluster. Structure -> [[CHROM, POS, ID, ALT, QUAL, FILTER, INFO], [CHROM, POS, ID, ALT, QUAL, FILTER, INFO], ...]
    '''
    metaclustersFields = []

    ## 3.2 Iterate over metaclusters
    for metacluster in metaclusters:

        posCluster = None
        BKP2 = None

        ## Collect insertion basic features
        CHROM = metacluster.ref
        # TODO SR: fix this
        #POS, CIPOS = metacluster.mean_pos()
        # Give always most left position first
        # If there are 2 bkps, show left first
        if metacluster.refRightBkp and metacluster.refLeftBkp:
            POS = min([metacluster.refRightBkp, metacluster.refLeftBkp])
            BKP2 = max([metacluster.refRightBkp, metacluster.refLeftBkp])
            # Check what side is POS right or left
            if POS == metacluster.refRightBkp:
                posCluster = 'RIGHT'
            elif POS == metacluster.refLeftBkp:
                posCluster = 'LEFT'
        
        # If there is only one clipping bkp
        else:
            POS = metacluster.refLeftBkp if metacluster.refLeftBkp != None else metacluster.refRightBkp
            # If it's reciprocal even though there is only one clipping bkp
            if metacluster.orientation == 'RECIPROCAL' and POS == metacluster.refLeftBkp and POS != metacluster.refLeftBkp:
                POS = min([metacluster.refLeftBkp, metacluster.end])
                BKP2 = max([metacluster.refLeftBkp, metacluster.end])
                if POS == metacluster.refLeftBkp:
                    posCluster = 'LEFT'
                elif POS == metacluster.end:
                    posCluster = 'RIGHT'
            elif metacluster.orientation == 'RECIPROCAL' and POS == metacluster.refRightBkp and POS != metacluster.refLeftBkp:
                POS = min([metacluster.refRightBkp, metacluster.beg])
                BKP2 = max([metacluster.refRightBkp, metacluster.beg])
                if POS == metacluster.refRightBkp:
                    posCluster = 'RIGHT'
                elif POS == metacluster.beg:
                    posCluster = 'LEFT'             

        # If there are no clipping bkp
        if POS == None:
            if metacluster.orientation == 'PLUS':
                POS = metacluster.end
                posCluster = 'RIGHT'
            elif metacluster.orientation == 'MINUS':
                POS = metacluster.beg
                posCluster = 'LEFT'
            elif metacluster.orientation == 'RECIPROCAL':
                POS = min([metacluster.beg, metacluster.end])
                BKP2 = max([metacluster.beg, metacluster.end])

        # NOTE SR: ID, ALT and QUAL are same for all metaclusters. So it could be possible to write them at the end, instead of in each metcluster
        ID = '.'
        ALT = '<INS>'
        QUAL = '.'

        # Set FILTER field
        # Check if there are FAILED filters at cluster level:
        filters = []
        for cluster in metacluster.rawSubclusters:
            if cluster.failedFilters:
                for failedFilter in cluster.failedFilters:
                    failAppend = failedFilter + '_CLUSTER'
                    if failAppend not in filters:
                        filters.append(failAppend)

        # Check if there are FAILED filters at metacluster level:
        if metacluster.failedFilters:
            for failedFilterMeta in metacluster.failedFilters:
                    failMetaAppend = failedFilterMeta + '_META'
                    if failMetaAppend not in filters:
                        filters.append(failMetaAppend)
        
        FILTER = 'PASS' if not filters else ','.join(filters)
        
        ## Collect extra insertion features to include at info field
        INFO = {}
        repeats = metacluster.repeatAnnot if hasattr(metacluster, 'repeatAnnot') else []

        discordants = [ discordant.readName for discordant in metacluster.events if discordant.type == 'DISCORDANT' ]
        clippings = [ clipping.readName for clipping in metacluster.events if clipping.type == 'CLIPPING' ]

        discordantsMAPQ = [ discordant.mapQual for discordant in metacluster.events if discordant.type == 'DISCORDANT' ]
        clippingsMAPQ = [ clipping.mapQual for clipping in metacluster.events if clipping.type == 'CLIPPING' ]

        discordantsMAPQMean = int(statistics.mean(discordantsMAPQ)) if discordantsMAPQ else None
        clippingsMAPQMean = int(statistics.mean(clippingsMAPQ)) if clippingsMAPQ else None

        # Set specific identities
        checkSpecIdent = {}
        for event in metacluster.events:
            if event.type != 'CLIPPING' and event.specificIdentity != None:

                if event.specificIdentity in checkSpecIdent.keys():
                    checkSpecIdent[event.specificIdentity] += 1
                else:
                    checkSpecIdent[event.specificIdentity] = 1
        
        if not checkSpecIdent:
            SPECIDENT = None
        else: 
            SPECIDENT = ",".join("{}:{}".format(k, checkSpecIdent[k]) for k in sorted(checkSpecIdent))
        
        # TODO SR: Avoid non-existing fields in another way and add interesting fields
        INFO['VTYPE'] = metacluster.mutOrigin
        # TODO SR: Put metacluster.identity as metacluster.SV_features['INS_TYPE'], so this INFO field can be deleted.
        INFO['IDENT'] = metacluster.identity
        #INFO['ITYPE'] = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        #INFO['MECHANISM'] = metacluster.SV_features['MECHANISM'] if 'MECHANISM' in metacluster.SV_features else None        
        #INFO['FAM'] = ','.join(metacluster.SV_features['FAMILY']) if ('FAMILY' in metacluster.SV_features and metacluster.SV_features['FAMILY']) else None
        #INFO['SUBFAM'] = ','.join(metacluster.SV_features['SUBFAMILY']) if ('SUBFAMILY' in metacluster.SV_features and metacluster.SV_features['SUBFAMILY']) else None
        #INFO['CIPOS'] = str(CIPOS[0]) + ',' + str(CIPOS[1]) 
        #INFO['CYTOID'] = ','.join(metacluster.SV_features['CYTOBAND']) if ('CYTOBAND' in metacluster.SV_features and metacluster.SV_features['CYTOBAND']) else None
        #INFO['NBEXONS'] = metacluster.SV_features['NB_EXONS'] if 'NB_EXONS' in metacluster.SV_features else None
        #INFO['SRCGENE'] = ','.join(metacluster.SV_features['SOURCE_GENE']) if 'SOURCE_GENE' in metacluster.SV_features else None
        #INFO['STRAND'] = metacluster.SV_features['STRAND'] if 'STRAND' in metacluster.SV_features else None
        #INFO['REGION'], INFO['GENE'] = metacluster.geneAnnot if hasattr(metacluster, 'geneAnnot') else (None, None)
        #INFO['REP'] = ','.join([repeat['family'] for repeat in repeats]) if repeats else None 
        #INFO['REPSUB'] = ','.join([repeat['subfamily'] for repeat in repeats]) if repeats else None   
        #INFO['DIST'] = ','.join([str(repeat['distance']) for repeat in repeats]) if repeats else None
        INFO['NBTOTAL'], INFO['NBTUMOR'], INFO['NBNORMAL'] = metacluster.nbEvents()[0], metacluster.nbEvents()[1], metacluster.nbEvents()[2]
        #INFO['NBSPAN'], INFO['NBCLIP'] = str(metacluster.nbINS), str(metacluster.nbCLIPPING)
        INFO['LEN'] = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        #INFO['CV'] = metacluster.cv
        #INFO['RTLEN'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        #INFO['TRUN5LEN'] = metacluster.SV_features['TRUNCATION_5_LEN'] if 'TRUNCATION_5_LEN' in metacluster.SV_features else None
        #INFO['TRUN3LEN'] = metacluster.SV_features['TRUNCATION_3_LEN'] if 'TRUNCATION_3_LEN' in metacluster.SV_features else None
        #INFO['FULL'] = metacluster.SV_features['IS_FULL'] if 'IS_FULL' in metacluster.SV_features else None
        #INFO['TDLEN'] = metacluster.SV_features['TRANSDUCTION_LEN'] if 'TRANSDUCTION_LEN' in metacluster.SV_features else None
        #INFO['INVLEN'] = metacluster.SV_features['INVERSION_LEN'] if 'INVERSION_LEN' in metacluster.SV_features else None
        #INFO['PERCR'] = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        #INFO['QHITS'] = None if metacluster.insertHits is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertHits.alignments ])
        #INFO['THITS'] = None if metacluster.insertHits is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertHits.alignments ])
        #INFO['RTCOORD'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        #INFO['POLYA'] = metacluster.SV_features['POLYA'] if 'POLYA' in metacluster.SV_features else None
        #INFO['INSEQ'] = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None
        INFO['DISC'] = ','.join(discordants) if discordants else None
        INFO['CLIP'] = ','.join(clippings) if clippings else None
        INFO['NBDISC'] = len(discordants)
        INFO['NBCLIP'] = len(clippings)
        INFO['ORIENT'] = metacluster.orientation
        INFO['BKP2'] = BKP2
        INFO['DISCMAPQ'] = discordantsMAPQMean
        INFO['CLIPMAPQ'] = clippingsMAPQMean
        INFO['SPECIDENT'] = SPECIDENT
        INFO['DISCDUP'] = 0 if metacluster.percDuplicates()[0] == 0 else "{:.2f}".format(metacluster.percDuplicates()[0])
        INFO['CLIPDUP'] = 0 if metacluster.percDuplicates()[1] == 0 else "{:.2f}".format(metacluster.percDuplicates()[1])
        # NOTE SR: Necessary for MEs annotation step.
        # TODO SR: When MEs step is futher, maybe this INFO['END'] = metacluster.end can be replaced by BKP2
        INFO['END'] = metacluster.end
        # NOTE SR: Necessary for keeping only MEs in annotation step.
        # NOTE SR: Perform repeats_annotation for both, MEIs and VIRUSES.
        #INFO['INTERNAL_ELEMENT'] = metacluster.events[0].element

        # Avoid showing BKP2 fields when there is not BKP2:
        # If integration is reciprocal (there are 2 bkps), put VCF tags of corresponding clusters
        if BKP2:
            if posCluster == 'RIGHT':
                INFO['BKPRSEQ'] = metacluster.repreRightSeq if metacluster.repreRightSeq != None else None
                INFO['BKP2RSEQ'] = metacluster.repreLeftSeq if metacluster.repreLeftSeq != None else None
                INFO['BKPCSEQ'] = metacluster.consRightSeq if metacluster.consRightSeq != None else None
                INFO['BKP2CSEQ'] = metacluster.consLeftSeq if metacluster.consLeftSeq != None else None
                INFO['INTBKP'] = ",".join("{}:{}".format(k, metacluster.intRightBkp[k]) for k in sorted(metacluster.intRightBkp)) if metacluster.intRightBkp else None
                INFO['INTBKP2'] = ",".join("{}:{}".format(k, metacluster.intLeftBkp[k]) for k in sorted(metacluster.intLeftBkp)) if metacluster.intLeftBkp else None
                INFO['CLIPTYPE'] = metacluster.rightClipType if metacluster.rightClipType != None else None
                INFO['CLIPTYPE2'] = metacluster.leftClipType if metacluster.leftClipType != None else None
            elif posCluster == 'LEFT':
                INFO['BKPRSEQ'] = metacluster.repreLeftSeq if metacluster.repreLeftSeq != None else None
                INFO['BKP2RSEQ'] = metacluster.repreRightSeq if metacluster.repreRightSeq != None else None
                INFO['BKPCSEQ'] = metacluster.consLeftSeq if metacluster.consLeftSeq != None else None
                INFO['BKP2CSEQ'] = metacluster.consRightSeq if metacluster.consRightSeq != None else None
                INFO['INTBKP'] = ",".join("{}:{}".format(k, v) for k, v in metacluster.intLeftBkp.items()) if metacluster.intLeftBkp else None
                INFO['INTBKP2'] = ",".join("{}:{}".format(k, v) for k, v in metacluster.intRightBkp.items()) if metacluster.intRightBkp else None
                INFO['CLIPTYPE'] = metacluster.leftClipType if metacluster.leftClipType != None else None
                INFO['CLIPTYPE2'] = metacluster.rightClipType if metacluster.rightClipType != None else None

        # If integration is not reciprocal (there are only one bkp)
        else:
            # Try with right bkp
            INFO['BKPRSEQ'] = metacluster.repreRightSeq if metacluster.repreRightSeq != None else None
            INFO['BKPCSEQ'] = metacluster.consRightSeq if metacluster.consRightSeq != None else None
            INFO['INTBKP'] = ",".join("{}:{}".format(k, v) for k, v in metacluster.intRightBkp.items()) if metacluster.intRightBkp else None
            INFO['CLIPTYPE'] = metacluster.rightClipType if metacluster.rightClipType != None else None
            # If there are not info at all of right bkp, show the info of left bkp.
            if INFO['BKPRSEQ'] == None and INFO['BKPCSEQ'] == None and INFO['INTBKP'] == None:
                INFO['BKPRSEQ'] = metacluster.repreLeftSeq if metacluster.repreLeftSeq != None else None
                INFO['BKPCSEQ'] = metacluster.consLeftSeq if metacluster.consLeftSeq != None else None
                INFO['INTBKP'] = ",".join("{}:{}".format(k, v) for k, v in metacluster.intLeftBkp.items()) if metacluster.intLeftBkp else None
                INFO['CLIPTYPE'] = metacluster.leftClipType if metacluster.leftClipType != None else None

        ## Create VCF variant object
        fields = [CHROM, POS, ID, ALT, QUAL, FILTER, INFO]

        metaclustersFields.append(fields)

    return metaclustersFields

def INS2VCF_SR(metaclustersFields, index, refLengths, source, build, species, VCFInfoFields, VCFREF, outName, outDir):
    '''
    Write INS calls into a VCF file

    Input:
        1. metaclustersFields: List of lists: containing VCF fields for each metacluster. Structure -> [[CHROM, POS, ID, ALT, QUAL, FILTER, INFO], [CHROM, POS, ID, ALT, QUAL, FILTER, INFO], ...]
        2. index: minimap2 index for the reference genome 
        3. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        4. source: software version used to generate the insertion calls
        5. build: reference genome build
        6. species: specie
        7. VCFInfoFields: List of INFO fields to display in output VCF.
        8. outName: Output file name
        9. outDir: Output directory

    Output: vcf file containing identified metaclusters
    '''

    ## 1. Initialize VCF 
    VCF = formats.VCF()

    ## 2. Create header
    ## Define info 
    info = {'VTYPE': ['.', 'String', 'Type of variant'], \
            'ITYPE': ['.', 'String', 'Type of structural variant'], \
            'MECHANISM': ['.', 'String', 'Insertion mechanism'], \
            'FAM': ['.', 'String', 'Repeat family'], \
            'SUBFAM': ['.', 'String', 'Repeat subfamily'], \
            'CIPOS': ['2', 'Integer', 'Confidence interval around POS for imprecise variants'], \
            'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
            'NBEXONS': ['1', 'Integer', 'Number of exons for a processed pseudogene insertion'], \
            'SRCGENE': ['.', 'String', 'Source gene for a processed psendogene insertion'], \
            'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
            'REGION': ['.', 'String', 'Genomic region where insertion occurs'], \
            'GENE': ['.', 'String', 'HUGO gene symbol'], \
            'REP': ['.', 'String', 'Families for annotated repeats at the insertion region'], \
            'REPSUB': ['.', 'String', 'Subfamilies for annotated repeats at the insertion region'], \
            'DIST': ['.', 'Integer', 'Distance between insertion breakpoint and annotated repeat'], \
            'NBTOTAL': ['1', 'Integer', 'Total number of insertion supporting reads'], \
            'NBTUMOR': ['1', 'Integer', 'Number of insertion supporting reads in the tumour'], \
            'NBNORMAL': ['1', 'Integer', 'Number of insertion supporting reads in the normal'], \
            'NBSPAN': ['1', 'Integer', 'Number of spanning supporting reads'], \
            'NBCLIP': ['1', 'Integer', 'Number of clipping supporting reads'], \
            'LEN': ['1', 'Integer', 'Insertion length'], \
            'CV': ['1', 'Float', 'Length coefficient of variation'], \
            'RTLEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
            'TRUN5LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
            'TRUN3LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
            'FULL': ['0', 'Flag', 'Full length mobile element'], \
            'TDLEN': ['1', 'Integer', 'Transduction length'], \
            'INVLEN': ['1', 'Integer', '5-inversion length'], \
            'PERCR': ['1', 'Float', 'Percentage of inserted sequence that has been resolved'], \
            'QHITS': ['.', 'String', 'Coordinates for inserted sequence hits on the reference'], \
            'THITS': ['.', 'String', 'Inserted sequence hits on the reference'], \
            'RTCOORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
            'POLYA': ['0', 'Flag', 'PolyA tail identified'], \
            'INSEQ': ['.', 'String', 'Inserted sequence'], \
            'DISC': ['.', 'String', 'Discordant reads supporting the insertion'], \
            'CLIP': ['.', 'String', 'Clipping reads supporting the insertion'], \
            'NBDISC': ['1', 'Integer', 'Number of discordant reads supporting the insertion'], \
            'NBCLIP': ['1', 'Integer', 'Number of clipping reads supporting the insertion'], \
            'DISCMAPQ': ['1', 'Integer', 'Average of MAPQ of discordant reads supporting the insertion'], \
            'CLIPMAPQ': ['1', 'Integer', 'Average of MAPQ of clipping reads supporting the insertion'], \
            'IDENT': ['.', 'String', 'Identity of the insertion'], \
            'ORIENT': ['.', 'String', 'Orientation of the insertion: RECIPROCAL, PLUS or MINUS'], \
            'BKP2': ['1', 'Integer', 'Reference position of second breakpoint of the insertion'], \
            'SPECIDENT': ['.', 'String', 'Specific identities and number of discordant reads supporting it. SpecificIdentity:#reads'], \
            'DISCDUP': ['1', 'Float', 'Percentage of discordant reads that are labeled as duplicates in input bam file.'], \
            'CLIPDUP': ['1', 'Float', 'Percentage of clipping reads that are labeled as duplicates in input bam file.'], \
            'REP': ['.', 'String', 'Families for annotated repeats at the insertion region'], \
            'REPSUB': ['.', 'String', 'Subfamilies for annotated repeats at the insertion region'], \
            'DIST': ['.', 'Integer', 'Distance between insertion breakpoint and annotated repeat'], \
            'REGION': ['.', 'String', 'Genomic region where insertion occurs'], \
            'GENE': ['.', 'String', 'HUGO gene symbol'], \
            'REFSeq': ['.', 'String', 'Reference sequence at insertion point (+- 10bp from insertion bkp).'], \
            'BKPCSEQ': ['.', 'String', 'Consensus sequence at BKP.'], \
            'BKP2CSEQ': ['.', 'String', 'Consensus sequence at BKP2.'], \
            'BKPRSEQ': ['.', 'String', 'Consensus sequence at BKP.'], \
            'BKP2RSEQ': ['.', 'String', 'Consensus sequence at BKP2.'], \
            'INTBKP': ['.', 'String', 'Coordinates of inserted sequence at BKP.'], \
            'INTBKP2': ['.', 'String', 'Coordinates of inserted sequence at BKP2.'], \
            'CLIPTYPE': ['.', 'String', 'Way that clipping where collected at BKP. DISC: Based on discordant clippings. BLAT: Based on BLAT search. REG: Based on cluster position.'], \
            'CLIPTYPE2': ['.', 'String', 'Way that clipping where collected at BKP2. DISC: Based on discordant clippings. BLAT: Based on BLAT search. REG: Based on cluster position.'], \
            }
            
    ## Create header
    VCF.create_header(source, build, species, refLengths, info, {}, [])
    
    ## 3. Add insertion calls to the VCF

    ## 3.1 Load reference index
    if VCFREF:
        reference = mp.Aligner(fn_idx_in=index) # comment as time consuming

    ## 3.2 Iterate over metaclusters fields
    for fields in metaclustersFields:

        # Get reference sequence (+- 5 bases)
        if VCFREF:
            REF = reference.seq(fields[0], int(fields[1]), int(fields[1]) + 1)
            fields[-1]['REFSeq'] = reference.seq(fields[0], int(fields[1]) - 10, int(fields[1]) + 10)
        else:
            REF = '.'
            fields[-1]['REFSeq'] = None

        # Insert REF to fields list [CHROM, POS, ID, ALT, QUAL, FILTER, INFO] -> [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO]
        fields.insert(3, REF)
        # Add empty dictionary because it is not multiVCF
        fields.append({})

        ## Add variant to the VCF
        INS = formats.VCF_variant(fields)
        
        VCF.add(INS)

    ## 4. Sort VCF
    VCF.sort()

    ## 5. Write VCF in disk
    IDS = VCFInfoFields

    VCF.write(IDS, {}, outName, outDir)

def INS2VCF(metaclusters, index, refLengths, source, build, species, outName, outDir):
    '''
    Write INS calls into a VCF file

    Input:
        1. metaclusters: list containing list of INS metaclusters
        2. index: minimap2 index for the reference genome 
        3. refLengths: Dictionary containing reference ids as keys and as values the length for each reference
        4. source: software version used to generate the insertion calls
        5. build: reference genome build
        6. species: specie
        7. outName: Output file name
        8. outDir: Output directory

    Output: vcf file containing identified metaclusters
    '''
    ## 1. Initialize VCF 
    VCF = formats.VCF()

    ## 2. Create header
    ## Define info 
    info = {'VTYPE': ['.', 'String', 'Type of variant'], \
            'ITYPE': ['.', 'String', 'Type of structural variant'], \
            'MECHANISM': ['.', 'String', 'Insertion mechanism'], \
            'FAM': ['.', 'String', 'Repeat family'], \
            'SUBFAM': ['.', 'String', 'Repeat subfamily'], \
            'GERMDB': ['.', 'String', 'List of germline variation databases where the variant is reported'], \
            'CIPOS': ['2', 'Integer', 'Confidence interval around POS for imprecise variants'], \
            'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
            'NBEXONS': ['1', 'Integer', 'Number of exons for a processed pseudogene insertion'], \
            'SRCGENE': ['.', 'String', 'Source gene for a processed psendogene insertion'], \
            'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
            'REGION': ['.', 'String', 'Genomic region where insertion occurs'], \
            'GENE': ['.', 'String', 'HUGO gene symbol'], \
            'REP': ['.', 'String', 'Families for annotated repeats at the insertion region'], \
            'REPSUB': ['.', 'String', 'Subfamilies for annotated repeats at the insertion region'], \
            'DIST': ['.', 'Integer', 'Distance between insertion breakpoint and annotated repeat'], \
            'NBTOTAL': ['1', 'Integer', 'Total number of insertion supporting reads'], \
            'NBTUMOR': ['1', 'Integer', 'Number of insertion supporting reads in the tumour'], \
            'NBNORMAL': ['1', 'Integer', 'Number of insertion supporting reads in the normal'], \
            'NBSPAN': ['1', 'Integer', 'Number of spanning supporting reads'], \
            'NBCLIP': ['1', 'Integer', 'Number of clipping supporting reads'], \
            'LEN': ['1', 'Integer', 'Insertion length'], \
            'CV': ['1', 'Float', 'Length coefficient of variation'], \
            'RTLEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
            'TRUN5LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
            'TRUN3LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
            'FULL': ['0', 'Flag', 'Full length mobile element'], \
            'TDLEN': ['1', 'Integer', 'Transduction length'], \
            'INVLEN': ['1', 'Integer', '5-inversion length'], \
            'PERCR': ['1', 'Float', 'Percentage of inserted sequence that has been resolved'], \
            'QHITS': ['.', 'String', 'Coordinates for inserted sequence hits on the reference'], \
            'THITS': ['.', 'String', 'Inserted sequence hits on the reference'], \
            'RTCOORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
            'POLYA': ['0', 'Flag', 'PolyA tail identified'], \
            'INSEQ': ['.', 'String', 'Inserted sequence'], \
            }
            
    ## Create header
    VCF.create_header(source, build, species, refLengths, info, {}, None)

    ## 3. Add insertion calls to the VCF
    ## 3.1 Load reference index
    reference = mp.Aligner(fn_idx_in=index) # comment as time consuming

    ## 3.2 Iterate over INS metaclusters
    for metacluster in metaclusters:

        ## Collect insertion basic features
        CHROM = metacluster.ref
        POS, CIPOS = metacluster.mean_pos()
        ID = '.'
        REF = reference.seq(CHROM, POS, POS + 1)
        ALT = '<INS>'
        QUAL = '.'
        FILTER = 'PASS' if not metacluster.failedFilters else ','.join(metacluster.failedFilters)
        
        ## Collect extra insertion features to include at info field
        INFO = {}
        repeats = metacluster.repeatAnnot if hasattr(metacluster, 'repeatAnnot') else []        

        INFO['VTYPE'] = metacluster.mutOrigin
        INFO['ITYPE'] = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        INFO['MECHANISM'] = metacluster.SV_features['MECHANISM'] if 'MECHANISM' in metacluster.SV_features else None        
        INFO['FAM'] = ','.join(metacluster.SV_features['FAMILY']) if ('FAMILY' in metacluster.SV_features and metacluster.SV_features['FAMILY']) else None
        INFO['SUBFAM'] = ','.join(metacluster.SV_features['SUBFAMILY']) if ('SUBFAMILY' in metacluster.SV_features and metacluster.SV_features['SUBFAMILY']) else None
        INFO['GERMDB'] = metacluster.germlineDb if hasattr(metacluster, 'germlineDb') else None       
        INFO['CIPOS'] = str(CIPOS[0]) + ',' + str(CIPOS[1]) 
        INFO['CYTOID'] = ','.join(metacluster.SV_features['CYTOBAND']) if ('CYTOBAND' in metacluster.SV_features and metacluster.SV_features['CYTOBAND']) else None
        INFO['NBEXONS'] = metacluster.SV_features['NB_EXONS'] if 'NB_EXONS' in metacluster.SV_features else None
        INFO['SRCGENE'] = ','.join(metacluster.SV_features['SOURCE_GENE']) if 'SOURCE_GENE' in metacluster.SV_features else None
        INFO['STRAND'] = metacluster.SV_features['STRAND'] if 'STRAND' in metacluster.SV_features else None
        INFO['REGION'], INFO['GENE'] = metacluster.geneAnnot if hasattr(metacluster, 'geneAnnot') else (None, None)
        INFO['REP'] = ','.join([repeat['family'] for repeat in repeats]) if repeats else None 
        INFO['REPSUB'] = ','.join([repeat['subfamily'] for repeat in repeats]) if repeats else None   
        INFO['DIST'] = ','.join([str(repeat['distance']) for repeat in repeats]) if repeats else None
        INFO['NBTOTAL'], INFO['NBTUMOR'], INFO['NBNORMAL'] = str(metacluster.nbTotal), str(metacluster.nbTumour), str(metacluster.nbNormal) 
        INFO['NBSPAN'], INFO['NBCLIP'] = str(metacluster.nbINS), str(metacluster.nbCLIPPING)
        INFO['LEN'] = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        INFO['CV'] = metacluster.cv
        INFO['RTLEN'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['TRUN5LEN'] = metacluster.SV_features['TRUNCATION_5_LEN'] if 'TRUNCATION_5_LEN' in metacluster.SV_features else None
        INFO['TRUN3LEN'] = metacluster.SV_features['TRUNCATION_3_LEN'] if 'TRUNCATION_3_LEN' in metacluster.SV_features else None
        INFO['FULL'] = metacluster.SV_features['IS_FULL'] if 'IS_FULL' in metacluster.SV_features else None
        INFO['TDLEN'] = metacluster.SV_features['TRANSDUCTION_LEN'] if 'TRANSDUCTION_LEN' in metacluster.SV_features else None
        INFO['INVLEN'] = metacluster.SV_features['INVERSION_LEN'] if 'INVERSION_LEN' in metacluster.SV_features else None
        INFO['PERCR'] = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        INFO['QHITS'] = None if metacluster.insertHits is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['THITS'] = None if metacluster.insertHits is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertHits.alignments ])
        INFO['RTCOORD'] = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        INFO['POLYA'] = metacluster.SV_features['POLYA'] if 'POLYA' in metacluster.SV_features else None
        INFO['INSEQ'] = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None

        ## Create VCF variant object
        fields = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, {}]

        ## Add variant to the VCF
        INS = formats.VCF_variant(fields)
        VCF.add(INS)
        
    ## 4. Sort VCF
    VCF.sort()

    ## 5. Write VCF in disk
    infoIds = ['VTYPE', 'ITYPE', 'MECHANISM', 'FAM', 'SUBFAM', 'GERMDB', 'CIPOS', 'CYTOID', \
           'NBEXONS', 'SRCGENE', 'STRAND', 'REGION', 'GENE', 'REP', 'REPSUB', 'DIST', \
           'NBTOTAL', 'NBTUMOR', 'NBNORMAL', 'NBSPAN', 'NBCLIP', 'LEN', 'CV', 'RTLEN', \
           'TRUN5LEN', 'TRUN3LEN', 'FULL', 'TDLEN', 'INVLEN', 'PERCR', \
           'QHITS', 'THITS', 'RTCOORD', 'POLYA', 'INSEQ']

    VCF.write(infoIds, [], outName, outDir)

def write_INS(INS_metaclusters, outFileName, outDir):
    '''
    Write INS calls into a tsv file

    Input:
        1. INS_metaclusters: list containing list of INS metaclusters
        2. outFileName: Output file name
        3. outDir: Output directory

    Output: tsv file containing identified metaclusters
    '''
    ## 1. Open output file 
    outFilePath = outDir + '/' + outFileName
    outFile = open(outFilePath, 'w')

    ## 2. Write header 
    row = "#ref \t beg \t end \t filters \t mutOrigin \t insType \t mechanism \t family \t subfamily \t cytobandId \t nbExons \t srcGene \t strand \t region \t gene \t families \t subfamilies \t distances \t nbTotal \t nbTumour \t nbNormal \t nbINS \t nbDEL \t nbCLIPPING \t length \t cv \t retroLen \t truncation5len \t truncation3len \t full \t transductionLen \t invLen \t percResolved \t qHits \t tHits \t retroCoord \t polyA \t insertSeq \n"
    outFile.write(row)

    ## 3. Write INS metaclusters 
    # For each metacluster
    for metacluster in INS_metaclusters:
                    
        ## General features 
        filters = 'PASS' if not metacluster.failedFilters else ','.join(metacluster.failedFilters)
        insType = metacluster.SV_features['INS_TYPE'] if 'INS_TYPE' in metacluster.SV_features else None
        mechanism = metacluster.SV_features['MECHANISM'] if 'MECHANISM' in metacluster.SV_features else None        
        strand = metacluster.SV_features['STRAND'] if 'STRAND' in metacluster.SV_features else None
        length = metacluster.consensusEvent.length if metacluster.consensusEvent is not None else None
        percResolved = metacluster.SV_features['PERC_RESOLVED'] if 'PERC_RESOLVED' in metacluster.SV_features else None
        qHits = None if metacluster.insertHits is None else ','.join([ 'insertedSeq' + ':' + str(alignment.qBeg) + '-' + str(alignment.qEnd) for alignment in metacluster.insertHits.alignments ])
        tHits = None if metacluster.insertHits is None else ','.join([ alignment.tName + ':' + str(alignment.tBeg) + '-' + str(alignment.tEnd) for alignment in metacluster.insertHits.alignments ])
        insert = metacluster.consensusEvent.pick_insert() if metacluster.consensusEvent is not None else None
        polyA = metacluster.SV_features['POLYA'] if 'POLYA' in metacluster.SV_features else None

        ## Repeat specific features
        family = ','.join(metacluster.SV_features['FAMILY']) if ('FAMILY' in metacluster.SV_features and metacluster.SV_features['FAMILY']) else None
        subfamily = ','.join(metacluster.SV_features['SUBFAMILY']) if ('SUBFAMILY' in metacluster.SV_features and metacluster.SV_features['SUBFAMILY']) else None
        retroCoord = metacluster.SV_features['RETRO_COORD'] if 'RETRO_COORD' in metacluster.SV_features else None
        
        ## Transduction specific features
        cytobandId = ','.join(metacluster.SV_features['CYTOBAND']) if ('CYTOBAND' in metacluster.SV_features and metacluster.SV_features['CYTOBAND']) else None

        ## Length features
        retroLen = metacluster.SV_features['RETRO_LEN'] if 'RETRO_LEN' in metacluster.SV_features else None
        full = metacluster.SV_features['IS_FULL'] if 'IS_FULL' in metacluster.SV_features else None
        transductionLen = metacluster.SV_features['TRANSDUCTION_LEN'] if 'TRANSDUCTION_LEN' in metacluster.SV_features else None
        invLen = metacluster.SV_features['INVERSION_LEN'] if 'INVERSION_LEN' in metacluster.SV_features else None
        truncation5len = metacluster.SV_features['TRUNCATION_5_LEN'] if 'TRUNCATION_5_LEN' in metacluster.SV_features else None
        truncation3len = metacluster.SV_features['TRUNCATION_3_LEN'] if 'TRUNCATION_3_LEN' in metacluster.SV_features else None

        ## Pseudogene specific features
        nbExons = metacluster.SV_features['NB_EXONS'] if 'NB_EXONS' in metacluster.SV_features else None
        srcGene = ','.join(metacluster.SV_features['SOURCE_GENE']) if 'SOURCE_GENE' in metacluster.SV_features else None

        ## Insertion region annotation
        repeats = metacluster.repeatAnnot if hasattr(metacluster, 'repeatAnnot') else []        
        families = ','.join([repeat['family'] for repeat in repeats]) if repeats else None 
        subfamilies = ','.join([repeat['subfamily'] for repeat in repeats]) if repeats else None   
        distances = ','.join([str(repeat['distance']) for repeat in repeats]) if repeats else None         
        region, gene = metacluster.geneAnnot if hasattr(metacluster, 'geneAnnot') else (None, None)

        # Write INS call into output file
        row = "\t".join([metacluster.ref, str(metacluster.beg), str(metacluster.end), str(filters), str(metacluster.mutOrigin), str(insType), str(mechanism), str(family), str(subfamily), str(cytobandId), str(nbExons), str(srcGene), str(strand), str(region), str(gene), str(families), str(subfamilies), str(distances), str(metacluster.nbTotal), str(metacluster.nbTumour), str(metacluster.nbNormal), str(metacluster.nbINS), str(metacluster.nbDEL), str(metacluster.nbCLIPPING), str(length), str(metacluster.cv), str(retroLen), str(truncation5len), str(truncation3len), str(full), str(transductionLen), str(invLen), str(percResolved), str(qHits), str(tHits), str(retroCoord), str(polyA), str(insert), "\n"])
        outFile.write(row)

    ## Close output file ##
    outFile.close()

def write_DISCORDANT(discordantClusters, outDir):
    '''
    Write DISCORDANT read pair calls into a tsv file
    '''

    ## 1. Open output file 
    fileName = "DISCORDANT_MEIGA.tsv"
    outFilePath = outDir + '/' + fileName
    outFile = open(outFilePath, 'w')

    ## 2. Write header 
    row = "#ref \t beg \t end \t clusterType \t family \t nbTotal \t nbTumour \t nbNormal \t repeats \t region \t gene \n"
    outFile.write(row)

    ## 3. Write clusters 
    # Iterate over the bins
    for discordantDict in discordantClusters:
        
        # Iterate over the cluster types
        for key, clusterList in discordantDict.items():

            key=key.replace("Simple_repeat", "Simple-repeat") #TODO find a more elegant solution. (javi: Meiga13)
            key=key.replace("Low_complexity","Low-complexity")
        
            orientation, clusterType, family = key.split('_')
            clusterType = orientation + '_' + clusterType

            # For each cluster from a given cluster type
            for DISCORDANT in clusterList:

                ## Collect info
                nbTotal, nbTumour, nbNormal = DISCORDANT.nbEvents()
                region, gene = DISCORDANT.geneAnnot if hasattr(DISCORDANT, 'geneAnnot') else ("None", "None")
            
                ## TEMPORARY: Only report discordant clusters that:
                # - Mate do not aligns on a retrotransposon
                if (family != 'None'):

                    # Write DISCORDANT cluster into output file
                    row = "\t".join([DISCORDANT.ref, str(DISCORDANT.beg), str(DISCORDANT.end), clusterType, family, str(nbTotal), str(nbTumour), str(nbNormal), str(DISCORDANT.repeatAnnot), region, gene, "\n"])
                    outFile.write(row)

def writeMetaclusters(metaclustersList, outDir):
    '''
    Write structural variation clusters into a tsv file
    '''

    ## Open output file ##
    fileName = "metaclusters.tsv"
    outFilePath = outDir + '/' + fileName
    outFile = open(outFilePath, 'w')

    for dictMetacluster in metaclustersList:
        if dictMetacluster != None:
            for metacluster,d2 in dictMetacluster.items():
                row = 'METACLUSTER: ' + str(metacluster) +' '+ str(len(metacluster.events)) +' '+ str(metacluster.ref) +' '+ str(metacluster.beg) +' '+ str(metacluster.end) +' '+ str(metacluster.intOrigin) + '\n'
                outFile.write(row)
                for event in metacluster.events:
                        if event.type == 'DISCORDANT':
                            row = str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' ' + str(event.identity) + ' ' + str(event.side) + ' ' + str(event.sample) + '\n'
                            outFile.write (row)
                        else:
                            row = str(metacluster) + ' ' + str(event.readName) + ' ' + str(event.ref) + ' ' + str(event.beg) + ' ' + str(event.type) + ' None ' + str(event.clippedSide) + ' ' + str(event.sample) + '\n'
                            outFile.write (row)
                for k,v in d2.items():
                    row = str(k) + ' = ' + str(v) + '\n'
                    outFile.write (row)  
         
    ## Close output file ##
    outFile.close()


def write_junctions(junctions, outFileName, outDir):
    '''
    Write BND junction calls into a tsv file

    Input:
        1. junctions: list containing list of BND junction objects 
        2. outFileName: Output file name
        3. outDir: Output directory
        
    Output: tsv file containing identified BND junctions
    '''
    ## 1. Open output file 
    outFilePath = outDir + '/' + outFileName
    outFile = open(outFilePath, 'w')

    ## 2. Write header 
    row = "#refA \t bkpA \t refB \t bkpB \t junctionType \t nbReadsTotal \t nbReadsTumour \t nbReadsNormal \t bridgeType \t family \t srcId \t supportTypes \t bridgeLen \t nbReadsBridge \t bridgeSeq \n"
    outFile.write(row)

    ## 3. Write BND junctions  
    # For each junction
    for junction in junctions:
        
        ## Collect into to write into file:
        refA = junction.metaclusterA.ref
        bkpA = junction.metaclusterA.bkpPos
        refB = junction.metaclusterB.ref
        bkpB = junction.metaclusterB.bkpPos
        junctionType = junction.junctionType()
        nbReadsTotal, nbReadsTumour, nbReadsNormal = junction.supportingReads()

        # Bridge info
        bridgeType = junction.bridge.bridgeType if junction.bridge is not None else None
        family = junction.bridge.family if junction.bridge is not None else None
        srcId = junction.bridge.srcId if junction.bridge is not None else None
        supportTypes = ','.join(junction.bridge.support_types()) if junction.bridge is not None else None
        bridgeLen = int(junction.bridge.bridgeLen) if junction.bridge is not None else None
        nbReadsBridge = junction.bridge.nbReads() if junction.bridge is not None else None       
        bridgeSeq = junction.bridge.bridgeSeq if junction.bridge is not None else None

        # Write BND junction call into output file
        row = "\t".join([refA, str(bkpA), refB, str(bkpB), junction.junctionType(), str(nbReadsTotal), str(nbReadsTumour), str(nbReadsNormal), str(bridgeType), str(family), str(srcId), str(supportTypes), str(bridgeLen), str(nbReadsBridge), str(bridgeSeq), "\n"])
        outFile.write(row)

    ## Close output file ##
    outFile.close()


def write_tdCounts_sureselect(clusterPerSrc, outDir):
    '''
    Compute transduction counts per source element and write into file

    Input:
        1. clusterPerSrc: 
        2. outDir: Output directory
        
    Output: tsv file containing transduction counts per source element
    '''
    ## 1. Compute number of transductions per source element
    counts = []

    for srcId, clusters in clusterPerSrc.items():
        nbTd = len(clusters)
        counts.append((srcId, nbTd))
    
    ## 2. Sort source elements in decreasing order of counts
    sortedCounts = sorted(counts, key=lambda x: x[1], reverse=True)

    ## 3. Write header 
    outFilePath = outDir + '/transduction_counts.tsv'
    outFile = open(outFilePath, 'w')

    row = "#srcId \t nbTransductions \n"
    outFile.write(row)

    ## 4. Write counts 
    for srcId, nbTd in sortedCounts:
        row = "\t".join([srcId, str(nbTd), "\n"])
        outFile.write(row)       

def write_tdCalls_sureselect(clustersPerSrc, outDir):
    '''
    Compute transduction counts per source element and write into file

    Input:
        1. clustersPerSrc: Dictionary containing one key per source element and the list of identified transduction clusters 
        2. outDir: Output directory
        
    Output: tsv file containing transduction counts per source element ordered by chromosome and then by start position
    '''
    ## 1. Write header 
    outFilePath = outDir + '/transduction_calls.tsv'
    outFile = open(outFilePath, 'w')
    row = "#ref \t beg \t end \t srcId \t nbReads \t nbDiscordant \t nbClipping \t readIds \n"
    outFile.write(row)

    ## 2. Generate list containing transduction calls
    call = None 
    calls = []

    # For each source element
    for srcId, clusters in clustersPerSrc.items(): 

        # For each cluster
        for cluster in clusters:
            readIds = ','.join(cluster.supportingReads()[3])
            call = [cluster.ref, str(cluster.beg), str(cluster.end), srcId, str(cluster.supportingReads()[0]), str(cluster.nbDISCORDANT()), str(cluster.nbSUPPLEMENTARY()), readIds]
            calls.append(call)

    ## 3. Sort transduction calls first by chromosome and then by start position
    calls.sort(key=lambda x: (x[0], int(x[1])))

    ## 4. Write transduction calls into output file
    for transduction in calls:
        row = "\t".join(transduction) + "\n"
        outFile.write(row)
    
    outFile.close()
    
    ## 5. Collapse calls when pointing to the same MEI. It happens when source elements are too close.
    
    if call is not None:
    
    	outFile = pybedtools.BedTool(outFilePath)
    
    	# Columns to collapse (without ref, beg and end columns)
    	colList = list(range(4, len(call)+1))
    	colFormat = ['collapse'] * (len(call) - 3)

    	mergedOutput = outFile.merge(c=colList, o=colFormat, header=True)
    	mergedOutput.saveas(outFilePath)
