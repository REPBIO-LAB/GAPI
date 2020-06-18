# MEIGA et all (Mobile Element Insertion Genome Analizer) and virus 

### Preliminary Meiga documentation

Method for the detection of somatic retrotransposition insertions from long-read data adn more ...

Evaluation of the algorithm shows that MEIGA is able to detect MEIs with a sensitivity and specicity over 90%, outperforming state-of-the-art MEI callers. Then, we used
MEIGA to identify somatic retrotranspositions in the cancer cell-lines. The number of L1 mobilizations increased by 60% as compared with short-read data and we observed striking high Alu activity rates in one of the samples while it is assumed that Alu repeats are rarely active in somatic cells. In addition, based on the analysis of diagnostic nucleotides we were able to infer for 18% of the solo L1 insertions their corresponding L1 source element. This analysis revealed novel hot L1s that mediate a reduced number of transductions. Finally, we have incorporated an additional module in MEIGA for the identication of viral insertions. Analysis of ONT sequencing data for one hepatocarcinoma cell-line identies novel viral insertions not previously identied with short-read approaches.


* Identification of MEI and viral insertions from long-read data 
* Integrates split-read and gapped alignment information 
* Single sample (germline) and tumour/matched normal (somatic) modes 
* Ready to use with any specie with reference genome available (human, mouse, dog and arabidopsis) 
* Fast and low memory usage (2h and 15Gb for 55x Pacbio data) 
* Compatible with mainstream long-read aligners (minimap2 and ngmlr) 

Three modes.

Long-reads 
* MEI insertion long read caller. Still need refinements
* Prototype for SV calling (including bridges)
* Database containing viral sequences for long and short read viral insertion calling
* Module for viral insertion detection from long reads

Short-reads

Virus detection (shor reads)


Surselect (short reads)
* Apply filters: Duplicates, FPs, germlines (report in the output)  

Dependencies.

* python 3
* mimimap2
* bwa
* samtools
* racon
* ncbi-blast
* bedtools
* perl
* annovar (perl annotate_variation.pl)

Python packages

* numpy
* scipy
* pysam
* mappy
* cigar
* pybedtools


Examples 

In linux console 

# 1. SINGLE mode for high coverage 2087 tumour cell line (Only chr22)
#######################################################################
program_dir=
working_dir=
MEIGA=$program_dir/MEIGA.py
tumour=$working_dir/H2009/H2009.all_merged.bam
normal=$working_dir/BL2009/BL2009.all_merged.sorted.bam
technology="NANOPORE"
reference=$working_dir/reference/0_homo_sapiens/hg_19_hs37d5/hs37d5.fa
databases=$working_dir/databases/Homo_sapiens/hg19/
annovarDir=$working_dir/databases/Homo_sapiens/hg19/annovarDb/
targetBins=$working_dir/minibam_DEL/L1_mediated_deletions_h2009.ok.bed
outDir=$working_dir/MEIGA/test/0.13.0/PAIRED_TARGET_NANOPORE_DELS/

mkdir $outDir
cd $outDir

$MEIGA $bam $technology $reference $databases --transduction-search --source-families 'L1' --refs '
22' -p 10 --gene-annot-dir $annovarDir -o $outDir 1> $outDir/MEIGA.out 2> $outDir/MEIGA.err" >> $outDir/launch.sh



# 2. PAIRED mode for L1-mediated deletions from 2009 tumour cell line (PCAWG paper)
####################################################################################
program_dir=
working_dir=
MEIGA=$program_dir/MEIGA.py
tumour=$working_dir/H2009/H2009.all_merged.bam
normal=$working_dir/BL2009/BL2009.all_merged.sorted.bam
technology="NANOPORE"
reference=$working_dir/reference/0_homo_sapiens/hg_19_hs37d5/hs37d5.fa
databases=$working_dir/databases/Homo_sapiens/hg19/
annovarDir=$working_dir/databases/Homo_sapiens/hg19/annovarDb/
targetBins=$working_dir/minibam_DEL/L1_mediated_deletions_h2009.ok.bed
outDir=$working_dir/MEIGA/test/0.13.0/PAIRED_TARGET_NANOPORE_DELS/

mkdir $outDir
cd $outDir

$MEIGA $tumour $technology $reference $databases --normalBam $normal --minReads 2 --minNo
rmalSupportingReads 1 --transduction-search --source-families 'L1' --refs '22' -p 10 --gene-annot-dir $annovarDir --targetBins $targetBins -o $outDir
 1> $outDir/MEIGA.out 2> $outDir/MEIGA.err"

# 3. SINGLE mode for genome in a bottle sample PACBIO data
###########################################################
# Genome in a bottle sample sequenced at 28X with Minion
MEIGA=$working_dir/MEIGA/MEIGA.py
bam=$working_dir/NA12878/pacbio/minimap2/pacbio_55x_minimap2_NA12878.sorted.bam
technology="PACBIO"
reference=$working_dir/reference/0_homo_sapiens/hg_19_hs37d5/hs37d5.fa
databases=$working_dir/MEIGA/databases/Homo_sapiens/hg19/
annovarDir=$working_dir/MEIGA/databases/Homo_sapiens/hg19/annovarDb/
refs='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'
outDir=$working_dir/MEIGA/test/0.9.0/SINGLE_WG_PACBIO_GIA

mkdir $outDir
cd $outDir

$MEIGA $bam $technology $reference $databases --transduction-search --source-families 'L1' --refs $
refs -p 22 --gene-annot-dir $annovarDir -o $outDir 1> $outDir/MEIGA.out 2> $outDir/MEIGA.err" >> $outDir/launch.sh





In a SLURM CLUSTER

# 1. PAIRED mode for L1-mediated deletions from 2009 tumour cell line (PCAWG paper)
####################################################################################
program_dir=
working_dir=
MEIGA=$program_dir/MEIGA.py
tumour=$working_dir/H2009/H2009.all_merged.bam
normal=$working_dir/BL2009/BL2009.all_merged.sorted.bam
technology="NANOPORE"
reference=$working_dir/reference/0_homo_sapiens/hg_19_hs37d5/hs37d5.fa
databases=$working_dir/databases/Homo_sapiens/hg19/
annovarDir=$working_dir/databases/Homo_sapiens/hg19/annovarDb/
targetBins=$working_dir/minibam_DEL/L1_mediated_deletions_h2009.ok.bed
outDir=$working_dir/MEIGA/test/0.13.0/PAIRED_TARGET_NANOPORE_DELS/

mkdir $outDir
cd $outDir

echo '#!/bin/sh' > $outDir/launch.sh
#echo "#SBATCH -A uscmg # ask for priority" >> $outDir/launch.sh
echo "#SBATCH -N 1 # (single node)" >> $outDir/launch.sh
echo "#SBATCH -n 10 # (10 processes in total)" >> $outDir/launch.sh
echo "#SBATCH -c 1 # (1 core per process)" >> $outDir/launch.sh
echo "#SBATCH -p cola-corta # (Ask for a partition. Several partitions can be specified separated by commas.)" >> $outDir/launch.sh
echo "#SBATCH -t 01:00:00 # (Maximum of 1h)" >> $outDir/launch.sh
echo "#SBATCH --error $outDir/job.err" >> $outDir/launch.sh
echo "#SBATCH --output $outDir/job.out" >> $outDir/launch.sh
echo 'source $working_dir/MEIGA/bin/activate' >> $outDir/launch.sh
echo "mprof run --include-children --multiprocess $MEIGA $tumour $technology $reference $databases --normalBam $normal --minReads 2 --minNo
rmalSupportingReads 1 --transduction-search --source-families 'L1' --refs '22' -p 10 --gene-annot-dir $annovarDir --targetBins $targetBins -o $outDir
 1> $outDir/MEIGA.out 2> $outDir/MEIGA.err" >> $outDir/launch.sh
chmod a+x $outDir/launch.sh

sbatch $outDir/launch.sh