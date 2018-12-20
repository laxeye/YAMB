## YAMB  (Yet Another Metagenome Binner) - semi-automatic pipeline for metagenomic contigs binnig.
It's based on t-Distributed Stochastic Neighbor Embedding (t-SNE) (van der Maaten, 2014) technique for dimensionality reduction, that uses tetranucleotide frequency distribution in metagenomic contigs and mean contig coverage, and HDBSCAN clustering method (Campello et.al., 2013).

## Requirements
* Bash
* Perl >= 5.10
* R >= 3.2
* R packages: ggplot2, Rtsne, getopt, dbscan
* samtools
* Read mapping software: **bowtie2** (by default) or minimap or bwa.
* esl-sfetch from easel for sequence extracting

To install R packages just run R environment and execute command:
*install.packages(c("ggplot2","Rtsne","getopt","dbscan"))*
To install samtools use your OS package manager (apt-get, yum etc.) or bioconda.

YAMB was tested on Ubuntu 16.04 and Ubuntu 18.04.

## How to use
At first, add the folder with executables to your PATH variable.

Run yamb.sh providing multifasta file with contigs, sequencing reads and output folder (optional, new folder *contigs filename*-yamb will be created by default):
./yamb.sh contigs.fasta R1.fastq\[.gz\] R2.fastq\[.gz\] \[output folder\]

In output folder you will find several files 
yamb-pp-\*-hdbscan.csv
where contig number, t-SNE perplexity and cluster count presented in the filenames.
These files are tab-separated sheets with columns:
Contig ID, 1st tSNE component, 2nd tSNE component, cluster #, normalised coverage, contig length
Corersponding files *yamb-pp-\*-hdbscan.png* are two-dimensional visualizations of contig fragments. Bins are coloured by different colours.

Contig acquisition is performed by easel utility "esl-sfetch". You can find metagenomic bins-yamb-N, where N is perplexity parameter.

It's strongly recommended to estimate completeness and contamination of bins.


## Validation

The pipeline was validated on synthetic methagenome, which ressembled an AMD community and included 7 bacterial and archeal species. Reads were produced by SimSeq and assembled by SPAdes. Five bins showed completeness between 94.49% and 100%, contamination was bellow 1.82%. One archael bin with mean coverage near 3x showed completeness 71.95% and 0% contamination, one bacterial bin 21.69% and 0% respectively. Low completeness of two bins can be explained by random read generation and low initial coverage (proportionally to their shares in community) which caused messing of genome parts. Binning results were compared to taxonomy assignment using best hits of nucleotide BLAST alignment against reference genomes: precision equaled 98.56%, recall equaled 98.13%.

YAMB binnnig results are comparable to the CONCOCT binning results on the real metagenomic data in terms of completeness and contamination of metagenomics bins (Korzhenkov A., unpublished data).


## Realisation

A. Korzhenkov, 2017-2018


## How to cite

The paper is still pending. 


## References

1. van der Maaten L.J.P. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15(Oct):3221-3245, 2014.
2. Campello, Ricardo JGB, Davoud Moulavi, and Jörg Sander. "Density-based clustering based on hierarchical density estimates." Pacific-Asia conference on knowledge discovery and data mining. Springer, Berlin, Heidelberg, 2013.
3. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043–1055.
