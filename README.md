# YAMB  (Yet Another Metagenome Binner) - semi-automatic pipeline for metagenomic contigs binning.
It's based on t-Distributed Stochastic Neighbor Embedding (t-SNE) (van der Maaten, 2014) technique for dimensionality reduction, that uses tetranucleotide frequency distribution in metagenomic contigs and mean contig coverage, and HDBSCAN clustering method (Campello et.al., 2013).

## Requirements
* Bash
* Perl >= 5.10
* R >= 3.2
* R packages: ggplot2, Rtsne, getopt, dbscan, RColorBrewer
* samtools
* Read mapping software: **bowtie2** (by default) or minimap or bwa.
* esl-sfetch from easel for sequence extracting

To install R packages just run R environment and execute command:

`install.packages(c("ggplot2","Rtsne","getopt","dbscan", "RColorBrewer"))`

To install samtools use your OS package manager (apt-get, yum etc.) or bioconda.

YAMB was tested on Ubuntu 16.04 and Ubuntu 18.04.

## How to use
At first, add the folder with executables to your PATH variable.

Run yamb.sh providing multifasta file with contigs, sequencing reads (forward/unpaired reads are obluigate, reverse and merged reads are optional) and output folder (optional, new folder *contigs filename*-yamb will be created by default):

`./yamb.sh -c <contigs.fasta> -f <R1.fastq[.gz]> [-r <R2.fastq[.gz]>] [-s <Merged.reads.fastq[.gz]>] [-o <output folder>] [-t <CPU threads>] [-m  <minimum contig length>]`

In output folder you will find several files `yamb-pp-<t-SNE perplexity>-hdbscan.csv` containing tab-separated values:
Contig ID, 1st tSNE component, 2nd tSNE component, cluster #, normalized coverage, contig length. 

Corresponding files `yamb-pp-<t-SNE perplexity>-hdbscan.png` are two-dimensional visualizations of contig fragments. Bins are colored by different colours.

Contig acquisition is performed by easel utility "esl-sfetch". You can find metagenomic bins in folder *bins-yamb-pp-N*, where N is perplexity parameter.

It's strongly recommended to estimate completeness and contamination of bins (e.g. [CheckM](https://github.com/Ecogenomics/CheckM>)).


## Validation

The pipeline was validated on synthetic metagenome, which resembled an AMD community and included 7 bacterial and archeal species. Reads were produced by SimSeq and assembled by SPAdes. Five bins showed completeness between 94.49% and 100%, contamination was bellow 1.82%. One archael bin with mean coverage near 3x showed completeness 71.95% and 0% contamination, one bacterial bin 21.69% and 0% respectively. Low completeness of two bins can be explained by random read generation and low initial coverage (proportionally to their shares in community) which caused messing of genome parts. Binning results were compared to taxonomy assignment using best hits of nucleotide BLAST alignment against reference genomes: precision equaled 98.56%, recall equaled 98.13%.

YAMB binnnig results are comparable to the CONCOCT binning results on the real metagenomic data in terms of completeness and contamination of metagenomics bins (Korzhenkov A., 2019).


## Realization

A. Korzhenkov, 2017-2019


## How to cite

Preprint is available on **bioRxiv**:

Korzhenkov, A. (2019). *YAMB: metagenome binning using nonlinear dimensionality reduction and density-based clustering*. BioRxiv, 521286 <https://doi.org/10.1101/521286>

## How to participate

Feel free to post an issue or clone the source and make a pull request.

## References

1. Van Der Maaten, L. (2014). Accelerating t-SNE using tree-based algorithms. The Journal of Machine Learning Research, 15(1), 3221-3245.
2. Campello, R. J., Moulavi, D., & Sander, J. (2013, April). Density-based clustering based on hierarchical density estimates. In Pacific-Asia conference on knowledge discovery and data mining (pp. 160-172). Springer, Berlin, Heidelberg.
3. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055.
