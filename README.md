## Genome Analysis API (GAPI)

GAPI is a set of python modules containing classes and functions for a wide range of purposes in genome analysis. These include, the processing of read alignments in BAM format, the handling of various standalone formats (including BED, VCF, FASTA, FASTQ, PAF, etc...) and more complex task, such as clustering of genome features based on their genomic positions and intersection of genomic intervals.
 
GAPI has been used for as the core library for the implementation of SVAN and MEIGA, two computational pipelines for the annotation and detection of structural variants and Mobile Element Insertions:

SVAN - Schloissnig et al., “Structural variation in 1,019 diverse humans based on long-read sequencing”, Nature, July 23, 2025, https://doi.org/10.1038/s41586-025-09290-7.

MEIGA - Zumalave et al., “Synchronous L1 retrotransposition events promote chromosomal crossover early in human tumorigenesis”, bioRxiv, August 27, 2024, https://doi.org/10.1101/2024.08.27.596794.

## Download 
Two different ways:

* Go to the releases tab and download the latest release. 

* Clone the git repository in case you want the latest version of the code:

```
# Move into the folder in which you want to clone the repositoy.
$ cd ~/apps
# Clone it.
$ git clone https://github.com/REPBIO-LAB/GAPI.git 
```

GAPI does not require any further installation step. It is written in Python and can be run as a standalone application on diverse Linux systems. 

## Requirements
1. Hardware:

    * 64-bits CPU

2. Software:

    * 64-bit Linux System
    * Python v3.5.4 or higher
    * paftools v.r755 (https://github.com/lh3/minimap2/tree/master/misc)
    * bwa-mem v0.7.17-r1188 (https://github.com/lh3/bwa)
    * minimap2 v2.10-r764-dirty (https://github.com/lh3/minimap2)

3. Python libraries 
    * pysam 
    * cigar
    * itertools
    * Biopython
    * subprocess
    * pandas
    * scipy
    * numpy
      
## Developers
GAPI was initially developed by Bernardo Rodríguez Martín and Eva García Álvarez, with contributions from Sonia Zumalave, Javier Temes and Martin Santamarina, in the Mobile Genomes Group (https://mobilegenomes.es) at the Center for Research in Molecular Medicine and Chronic Diseases (CiMUS) (2018–2022). Since 2024, it has been maintained and further developed at the Repetitive DNA Biology Lab (https://www.crg.eu/en/content/research/independent-fellow/bernardo-rodriguez-martin) at the Centre for Genomic Regulation (CRG).

## License
GAPI is distributed under GPL-3.0 License.

## Contact
Please open a case on the Github page for problems.

