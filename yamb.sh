#!/usr/bin/env bash
if ( [[ -z "$1" ]] || [[ $1 == -h* ]] || [[ $1 == -v* ]] )
then
	echo "Yet Another Metagenome Binner. Code by A. Korzhenkov. 2017-2019"
	echo "Usage: $0 -c <contigs> -f <forward reads> [ -r <reverse reads>] [ -s <merged reads> ] [ -o <output folder> ]"
	echo "[ -t <CPU threads> ] [ -l <contig fragment length> ] [ -m  <minimum contig length> ]"
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
    -s|--merged)
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
    -l|--length)
      LENGTH=$2
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

if ! ( [[ -f $contigs ]] )
then
  echo "Contigs \"$contigs\" not found"
  exit 1
fi

if ! ( [[ -f $R1 ]] ); then
  echo "Error: forward reads not found!"
elif ! ( [[ -f $R2 ]] || [[ -z $R2 ]] ); then
  echo "Error: reverse reads not found!"
elif ! ( [[ -f $RS ]] || [[ -z $RS ]] ); then
  echo "Error: merged reads not found!"
fi

#Default values for
if [ -z $THREADS ]; then THREADS=1; fi
if [ -z $MINLENGTH ]; then MINLENGTH=1000; fi
if [ -z $output ]; then output=$R1-yamb; fi
if [ -z $LENGTH ]; then LENGTH=10000; fi
logfile="$output/log.txt"

if [[ ! -e $output ]]; then mkdir $output; fi

BIN_PATH="`dirname \"$0\"`"

$BIN_PATH/cut-contigs.pl $contigs $LENGTH > "$output/cut.contigs.fna"
contigs="$output/cut.contigs.fna"

if (which bowtie2 > /dev/null); then
	echo -e "$(date +"%T")\tMapping reads with Bowtie2" | tee $logfile
	bowtie2-build -q $contigs "$output/idx"
  if ( [[ -z $RS ]] ); then
    if [[ -z $R2 ]] ; then
      bowtie2 -p $THREADS -q -x "$output/idx" -U $R1 -S "$output/mapping.sam"
    else
      bowtie2 -p $THREADS -q -x "$output/idx" -1 $R1 -2 $R2 -S "$output/mapping.sam"
    fi
  else
    bowtie2 -p $THREADS -q -x "$output/idx" -1 $R1 -2 $R2 -U $RS -S "$output/mapping.sam"
  fi

elif (which minimap2 > /dev/null); then
	echo -e "$(date +"%T")\tMapping reads with minimap2" | tee -a $logfile
	minimap2 -x sr -a -t $THREADS $contigs $R1 $R2 > "$output/mapping.sam"
elif (which bwa > /dev/null); then
	echo -e "$(date +"%T")\tMapping reads with bwa" | tee -a $logfile
	bwa index $contigs
	bwa mem -t $THREADS $contigs $R1 $R2 > "$output/mapping.sam"
else
	echo -e "$(date +"%T")\tRead mappers not found! Exitting." | tee -a $logfile
	exit 1
fi

echo -e "Converting SAM to BAM, sorting and indexing." | tee -a $logfile
mapping=$output"/mapping.bam"
samtools view -F 4 -u -b $output"/mapping.sam" | samtools sort - > $mapping && rm $output"/mapping.sam"
samtools index $mapping

echo -e "$(date +"%T")\tIndexing fasta." | tee -a $logfile
samtools faidx $contigs

echo -e "$(date +"%T")\tCoverage calculation." | tee -a $logfile
awk -v "OFS=\t" '{print $1,"0",$2}' $contigs.fai > $output/contigs.bed
samtools bedcov -Q 20 $output/contigs.bed $mapping | awk -F "\t" -v OFS="\t" '{print $1,$3,$4/$3}' > $output/coverage.csv
echo -e "$(date +"%T")\tTetramer occurence estimation." | tee -a $logfile
$BIN_PATH/tetramers.pl $contigs | grep -v "^#" > $output/tetramers.csv
echo -e "$(date +"%T")\tDatafile generation." | tee -a $logfile
paste $output/tetramers.csv <(awk -v FS="\t" -v OFS="\t" '{print $3/100, $2}' $output/coverage.csv) | awk '$NF > "'"$MINLENGTH"'"' > $output/data.csv
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
	echo -e "$(date +"%T")\tError: Unable to locate esl-sfetch, no fasta files may be produced" | tee -a $logfile
fi

cd ..
echo -e "$(date +"%T")\tBinning complete! Check $output/ for results!" | tee -a $logfile
