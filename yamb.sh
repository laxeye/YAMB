#!/usr/bin/env bash
if ( [[ -z "$1" ]] || [[ $1 == -h* ]] || [[ $1 == -v* ]] )
then
	echo "Yet Another Metagenome Binner. Code by A. Korzhenkov. 2017-2019"
	echo "Usage: $0 -c <contigs> -f <Forward reads> -r <Reverse reads> [ -o <output folder> ] [ -t <CPU threads> ] [ -m  minimum contig length ]"
	echo "Reads may be gzipped."
	echo "Multi-threading makes sense only for mapping."
	exit
fi

if ! ( which samtools > /dev/null ); then echo "You should install samtools.";exit 1; fi


while (( "$#" )); do
  case "$1" in
    -f|--forward)
      R1=$2
      shift 2
      ;;
    -r|--reverse)
      R2=$2
      shift 2
      ;;
    -s|--single)
      RS=$2
      shift 2
      ;;
    -c|--contigs)
      contigs=$2
      shift 2
      ;;
    -o|--output)
      output=$2
      shift 2
      ;;
    -t|--threads)
      THREADS=$2
      shift 2
      ;;
    -m|--minimum-length)
      MINLENGTH=$2
      shift 2
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

if ! ( ( [[ -f $R1 ]] && [[ -f $R2 ]] ) || ( [[ -f $R1 ]] && [[ -z $R2 ]] ) && [[ -f $contigs ]] )
then
	echo "Input files not found!"
	exit 1
fi

#Default values for
if [ -z $THREADS ]; then THREADS=1; fi
if [ -z $MINLENGTH ]; then MINLENGTH=1000; fi
if [ -z $output ]; then output=$R1-yamb; fi

if [[ ! -e $output ]]; then mkdir $output; fi

BIN_PATH="`dirname \"$0\"`"

$BIN_PATH/cut-contigs.pl $contigs 10000 > "$output/cut.contigs.fna"
contigs="$output/cut.contigs.fna"

if (which bowtie2 > /dev/null);then
	echo "Mapping reads with Bowtie2"
	bowtie2-build -q $contigs "$output/idx"
  if (  [[ ! -z $RS ]]  && [[ -f $RS ]] )
  then
    bowtie2 -p $THREADS -q -x "$output/idx" -1 $R1 -2 $R2 -U $RS -S "$output/mapping.sam"
  else
	  bowtie2 -p $THREADS -q -x "$output/idx" -1 $R1 -2 $R2 -S "$output/mapping.sam"
  fi
      
elif (which minimap2 > /dev/null); then
	echo "Mapping reads with minimap2"
	minimap2 -x sr -a -t $THREADS $contigs $R1 $R2 > "$output/mapping.sam"
elif (which bwa > /dev/null); then
	echo "Mapping reads with bwa"
	bwa index $contigs
	bwa mem -t $THREADS $contigs $R1 $R2 > "$output/mapping.sam"
else
	echo "Read mappers not found! Exitting."
	exit 1
fi

echo "Converting SAM to BAM, sorting and indexing."
mapping=$output"/mapping.bam"
samtools view -F 4 -u -b $output"/mapping.sam" | samtools sort - > $mapping && rm $output"/mapping.sam"
samtools index $mapping

echo "Indexing fasta."
samtools faidx $contigs

echo "Coverage calculation."
awk -v "OFS=\t" '{print $1,"0",$2}' $contigs.fai > $output/contigs.bed
samtools bedcov -Q 20 $output/contigs.bed $mapping | awk -F "\t" -v OFS="\t" '{print $1,$3,$4/$3}' > $output/coverage.csv
echo "Tetramer occurence estimation."
$BIN_PATH/tetramers.pl $contigs | grep -v "^#" > $output/tetramers.csv
echo "Datafile generation."
paste $output/tetramers.csv <(awk -v FS="\t" -v OFS="\t" '{print $3/100, $2}' $output/coverage.csv) | awk '$259 > 1000' > $output/data.csv
cd $output
#Default parameters here. For more information run tsne-kmean.r
$BIN_PATH/tsne-clust.r -i data.csv -m $MINLENGTH
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
			esl-sfetch -f cut.contigs.fna <(awk '{if($4 == '$i') print $1}' $f) | perl -lne 'if(/>(.+)_\[(\d+)\]/){$a=$1;$i=$2;if(!($a=~/$ao/ and $b==$i-1)){$ao=$a;print "$_"}$b=$i}else{print}' > $D/bin-$i.fna
		done
	done
else
	echo "Unable to locate esl-sfetch, no fasta files may be produced"
fi

cd ..
echo "Mission complete! Check $output/ for results!"
