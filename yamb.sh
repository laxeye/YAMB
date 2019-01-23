#!/usr/bin/env bash
if ( [[ -z "$1" ]] || [[ $1 == -h* ]] || [[ $1 == -v* ]] )
then
	echo "Yet Another Metagenome Binner"
	echo "Code: A. Korzhenkov. 2017-2018"
	echo "Please run as: $0 <contig.fasta> <Read1.fastq[.gz]> <Read2.fastq[.gz]> [output folder]"
	exit
fi

if ! ( which samtools > /dev/null ); then echo "You should install samtools.";exit 1; fi


if ! ( [[ -f $1 ]] && [[ -f $2 ]] && [[ -f $3 ]] )
then
	echo "Input files not found!"
	exit 1
fi

THREADS=2
contigs=$1
R1=$2
R2=$3
if [ $# -gt 3 ]; then output=$4; else output=$1-yamb; fi

if [[ ! -e $output ]]; then mkdir $output; fi

./cut-contigs.pl $contigs 10000 > "$output/cut.contigs.fna"
contigs="$output/cut.contigs.fna"

if (which bowtie2 > /dev/null);then
	echo "Mapping reads with Bowtie2"
	bowtie2-build -q $contigs "$output/idx"
	bowtie2 -p $THREADS -q -x "$output/idx" -1 $R1 -2 $R2 -S "$output/mapping.sam"
elif (which minimap2 > /dev/null); then
	echo "Mapping reads with minimap2"
	minimap2 -ax sr $contigs $R1 $R2 > "$output/mapping.sam"
elif (which bwa > /dev/null); then
	echo "Mapping reads with bwa"
	bwa index $contigs
	bwa mem $contigs $R1 $R2 > "$output/mapping.sam"
else
	echo "Read mappers not found! Exitting."
	exit 1
fi

echo "Converting SAM to BAM, sorting and indexing."
mapping=$output"/mapping.bam"
samtools view -F 4 -1 -b $output"/mapping.sam" | samtools sort - > $mapping && rm $output"/mapping.sam"
samtools index $mapping

echo "Indexing fasta."
samtools faidx $contigs

echo "Coverage calculation."
awk -v OFS="\t" '{print $1,"0",$2}' $contigs.fai > $output/contigs.bed
samtools bedcov -Q 20 $output/contigs.bed $mapping | awk -F "\t" -v OFS="\t" '{print $1,$3,$4/$3}' > $output/coverage.csv
echo "Tetramer occurence estimation."
tetramers.pl $contigs | grep -v "^#" > $output/tetramers.csv
echo "Datafile generation."
paste $output/tetramers.csv <(awk -v FS="\t" -v OFS="\t" '{print $3/100, $2}' $output/coverage.csv) | awk '$259 > 1000' > $output/data.csv
cd $output
#Default parameters here. For more information run tsne-kmean.r
tsne-clust.r -i data.csv
if (which esl-sfetch > /dev/null)
then
	esl-sfetch --index cut.contigs.fna
	for f in *scan.csv
	do 
		N=$(cut -f 4 $f | sort -nr | head -1)
		D=${f//-hdbscan.csv/}
		D="bins-$D"
		mkdir $D
		for i in $(seq 0 $N)
		do
			esl-sfetch -f cut.contigs.fna <(awk '{if($4 == '$i') print $1}' $f) | perl -lne 'if(/>(.+)_(\d+)/){$a=$1;$i=$2;if(!($a=~/$ao/ and $b==$i-1)){$ao=$a;$b=i;print "$_"}else{$b=$i}}else{print}' > $D/bin-$i.fna
		done
	done
else
	echo "Unable to locate esl-sfetch, no fasta files may be produced"
fi

cd ..
echo "Mission complete! Check $output/ for results!"
