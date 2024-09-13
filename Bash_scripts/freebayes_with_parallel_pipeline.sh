#! /bin/bash

REF=${1}
NJ=$(( $2/2 ))     					# Split into NJ parts with equal bp. Half of threads available
BAM=$3								# Input preprocessed BAM file
sample=$4
INTERVAL=$5							# The original target bed file
outVCF=$6
functions_dir=$7

fbayes_dir=`pwd`						# Set pwd as a variable
log_dir=$fbayes_dir/vc_logs
if [ ! -d $log_dir ]; then mkdir $log_dir; fi

# Manage threads and resource allocation
parallelthreads=$(( $2/$NJ ))
javamem=$(( $(free -g | grep "Mem" | awk '{print $2}') / 2 ))		# Set Java maximum mem to 50% of total RAM available

# Get record time function from functions.sh
source $functions_dir/./functions.sh
declare -a alltimes
starttime=""
endtime=""

# Get ref prefix by removing .fasta and show inputs
suffix=".fasta"
ref_prefix=${REF%"$suffix"}
# Check if .dict is present with the reference genome
if [ ! -f ${ref_prefix}.dict ]; then
	echo ".dict file for reference genome is not present. Creating .dict file for $REF..."
	time_func
	/usr/bin/time -v -o $log_dir/refDict.stderr \
		docker run -v /export:/export broadinstitute/gatk gatk \
		--java-options "-Xmx${javamem}g -XX:ParallelGCThreads=$NJ" \
		CreateSequenceDictionary \
		-R $REF \
		-O ${ref_prefix}.dict \
		&& echo "Finish creating dictionary for reference genome" || echo "Failed"
	time_func
else
	echo ".dict file exists. Using ${ref_prefix}.dict"
fi

# We need to create a file with all the chromsome lengths. I do this from a novoindex of the reference but you could also make from the @SQ headers in a BAM file.
time_func
chrLengths_file=$fbayes_dir/chromosomeLengths_${sample}.txt
/usr/bin/time -v -o $log_dir/preparation.stderr \
	docker run -v /export:/export staphb/samtools samtools view -H $BAM |
	 grep "^@SQ" | 
	 awk 'BEGIN { OFS="\t"; } { split($2, a, ":"); split($3, b, ":"); print a[2], b[2]; }' \
	 > $chrLengths_file	# Table chromosome length for bedtools slop

# Add 100bp slop to each target and merge overlapping targets.
tmpbed=$fbayes_dir/tmp.bed
expTargets_file=$fbayes_dir/expandedTargets_${sample}.bed 
/usr/bin/time -v -a -o $log_dir/preparation.stderr \
	docker run -v /export:/export staphb/bedtools bedtools slop -b 100 -g $chrLengths_file -i $INTERVAL > $tmpbed
/usr/bin/time -v -a -o $log_dir/preparation.stderr \
	docker run -v /export:/export staphb/bedtools bedtools merge -d 50 -i $tmpbed > $expTargets_file
rm -f $tmpbed

export PREFIX=j		# Prefix of file name for split targets files. This will make files like j1.bed, j2.bed, etc
rm -f $fbayes_dir/${PREFIX}*.bed  # Remove old files
total_length=$(awk '{ l += ($3 - $2); } END { printf "%.0f", l }' < $expTargets_file)
echo $total_length
export BP=$(( 1 + total_length / NJ ))  ## Add 1 to make sure we get all
echo $BP

awk -v bp=$BP -v p=$PREFIX -v fbdir=$fbayes_dir 'BEGIN {print "BP:" bp} { l += ($3 - $2); j=int(l/bp); print $0 >> fbdir "/" p j ".bed"; }' $expTargets_file
time_func

echo "Finish spliting"
echo "***************************************\n"

tmp_dir=$fbayes_dir/tmp_vc
mkdir $tmp_dir

echo "Start Freebayes with $NJ threads"
time_func
ls $fbayes_dir/j*.bed | parallel -j $NJ "/usr/bin/time -v -o {}_checkRes.stderr \
	docker run -v /export:/export staphb/freebayes freebayes \
	-f $REF -t {} $BAM > {}.vcf \
	2>>$fbayes_dir/freebayes_${sample}.log" \
	&& echo "Done" || echo "Fail"

mergeVCF=${outVCF%".gz"}
grep "^#" $fbayes_dir/j0.bed.vcf > $mergeVCF
for C in `ls $fbayes_dir/j*.vcf | sort`
do
   awk '$6 >= 1 && !/^#/ { print $0 }' $C >> $mergeVCF
done && echo "Done" || echo "Fail"
time_func

# bgzip the merged VCF
docker run -v /export:/export staphb/htslib bgzip -i -c -@ `nproc` $mergeVCF > $outVCF && \
	echo Done || echo Failed

if [ -f $outVCF ]; then
	echo -e "Variant Calling done. \nOutput VCF: ${outVCF}"
	echo "Removing the split interval files and parallel VCF files" && \
	rm -f ${fbayes_dir}/${PREFIX}*.bed && \
	rm -f ${fbayes_dir}/${PREFIX}*.bed.vcf && \
	rm -f ${fbayes_dir}/${PREFIX}*.bed.vcf.idx && \
	rm -rf $tmp_dir && \
	echo "Generated interval files and parallel VCF files removed"
else
	echo "Output file not found. Please check the directories."
	exit 1
fi
echo -e "\n***************************************\n"

$functions_dir/./view_stderrs.sh $fbayes_dir > $fbayes_dir/report_stderrs.txt && \
	echo "Done" || echo "Fail"

echo ${alltimes[@]}
countTotalTime
echo "\t\t\t\t\t\t\t-----End of execution-----"
