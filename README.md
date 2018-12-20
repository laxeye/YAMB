## YAMB  (Yet Another Metagenome Binner) - semi-automatic pipeline for metagenomic contigs binnig.
It's based on t-Distributed Stochastic Neighbor Embedding (t-SNE) (van der Maaten, 2014) technique for dimensionality reduction, that uses tetranucleotide frequency distribution in metagenomic contigs and mean contig coverage, and HDBSCAN clustering method (Campello et.al., 2013).

## Requirements
* Bash
* Perl >= 5.10
* R >= 3.2
* R packages: ggplot2, Rtsne, getopt, dbscan
* samtools
* Read mapping software: **bowtie2** (default) or minimap or bwa.
* esl-sfetch from easel for sequence extracting

To install R packages just run R environment and execute command:
*install.packages(c("ggplot2","Rtsne","getopt","dbscan"))*
To install samtools use your OS package manager (apt-get, yum etc.) or bioconda.

YAMB was tested on Ubuntu 16.04 and Ubuntu 18.04.

## How to use
At first, add the folder with executables to your PATH variable.

Cut metagenome contigs using cut-contigs.pl:
./cut-contigs.pl assembly.fna > cutted.assembly.fna
By default window size equals to 20 000 nucleotides.

Map reads to cutted contigs using Your prefered mapping software

Sort and index resulting mapping file using samtools. For example if You have a mapping file in SAM format called "mapping.sam" execute next commands:
samtools view -b -F 4 mapping.sam | samtools sort - > mapping.bam
samtools index mapping.bam

Run yamb.sh providing multifasta file with contigs, sequencing reads and output folder (it will be cleaned):
./yamb.sh contigs.fasta R1.fastq.gz R2.fastq.gz binning-output


In output folder you will find several files 
yamb-pp-\*-hdbscan.csv
where contig number, t-SNE perplexity and cluster count presented in the filenames.
These files are tab-separated sheets with columns:
Contig ID, 1st tSNE component, 2nd tSNE component, cluster #, normalised coverage, contig length
Corersponding files
yamb-pp-\*-hdbscan.png
are visualization of contigs in two-dimensional space, where size of point depends on contig size, transparence - on mean contig coverage. Bins are coloured by different colours and labeled with a number. Please notice that K-means clustering alghorhytm may resilt in bin overdividing or bins merging, thats why you need to check bin size and its graphical imaging.


To get a list of contigs in any bin you can use next command, where You should substitute X with bin number:
awk -F"\t" '{if($4 == X) print $1}' yamb-NK-pp-NN-cl-NN.csv | sed 's/"//g' > binX.txt

Contig acquisition performed by easel utility "esl-sfetch" according it's instructions.

It's strongly recommended to estimate completeness and contamination of bins.


## Validation

The pipeline was validates using synthetic methagenome, which ressembled an acidic mine drainage community of 7 
bacterial and archeal species. Five bins showed completeness between 94.49% and 100%, contamination was bellow 1.82%. One archael bin with coverage near 3x showed completeness 71.95% and 0% contamination, one bacterial bin 21.69% and 0% respectively. Bad results can be explained by random read generation and low initial coverage (in proportion to share in community) which caused genome parts missing.


## Realisation

A. Korzhenkov, 2017-2018
With kind help of:
S. Toshchakov, O. Golyshina and A. Tepliuk.


## References

1. van der Maaten L.J.P. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15(Oct):3221-3245, 2014.
2. Campello, Ricardo JGB, Davoud Moulavi, and Jörg Sander. "Density-based clustering based on hierarchical density estimates." Pacific-Asia conference on knowledge discovery and data mining. Springer, Berlin, Heidelberg, 2013.
3. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043–1055.
