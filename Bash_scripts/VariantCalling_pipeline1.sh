#! /bin/bash
# Script original name: Split_n_runParallel_v1.1.sh
# This script has seven parameters: ./scriptname <Referenc_genome.fasta> [threads/files_to_split] <Input.bam> <sample> <Region/Interval.bed> <Output.vcf> </path/to/functions.sh/directory>
#	1) Reference genome (.fasta)
#	2) No. of files to split and to run // threads;
#	3) Input preprocessed BAM file; 
#	4) Sample (eg: HG001);
#	5) Interval file (.bed);
#	6) Output - unVQSR (.vcf)
#	7) Scripts directory (to access functions.sh)
# Example:
# - ./scriptname.sh preprocessed.bam 32 /path/to/interval.bed HG001 script_dir
#
# Note: Now this can run anywhere, but have to specify interval file
# Note: Only half of the threads available/provided will be used because running variant calling in parallel can be very CPU intensive, causing high load average. During my testing in 32 cores 64 threads server, using 64 threads will lead to some of the jobs having error. Thus, I split it in half, to align to the number of cores available.

REF=${1}
NJ=$(( $2/2 ))     					# Split into NJ parts with equal bp. Half of threads available
BAM=$3								# Input preprocessed BAM file
sample=$4
INTERVAL=$5							# The original target bed file
outVCF=$6
functions_dir=$7

vc_dir=`pwd`						# Set pwd as a variable
log_dir=$vc_dir/vc_logs
rm -rf $log_dir
mkdir $log_dir

# Manage threads and resource allocation
parallelthreads=$(( $2/$NJ ))
javamem=$(( $(free -g | grep "Mem" | awk '{print $2}') / 2 ))		# Set Java maximum mem to 50% of total RAM available
echo "[Notice]: Using ${javamem}Gb of RAM"

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
	echo ".dict file for reference genome is not present."
	echo "[Process]: Creating .dict file for $REF..."
	time_func
	/usr/bin/time -v -o $log_dir/refDict.stderr \
		gatk --java-options "-Xmx${javamem}g -XX:ParallelGCThreads=$NJ" \
		CreateSequenceDictionary \
		-R $REF \
		-O ${ref_prefix}.dict \
		&& echo "[Status]: Finish creating dictionary for reference genome" || echo "Failed"
	time_func
else
	echo ".dict file exists. Using ${ref_prefix}.dict"
fi

# We need to create a file with all the chromsome lengths. I do this from a novoindex of the reference but you could also make from the @SQ headers in a BAM file.
echo "[Process]: Splitting BED file for running HaplotypeCaller in parallel..."
time_func
chrLengths_file=$vc_dir/chromosomeLengths.txt
/usr/bin/time -v -o $log_dir/preparation.stderr \
	samtools view -H $BAM | grep "^@SQ" | \
	awk 'BEGIN { OFS="\t"; } { split($2, a, ":"); split($3, b, ":"); print a[2], b[2]; }' \
	> $chrLengths_file	# Table chromosome length for bedtools slop

tmpbed=$vc_dir/tmp.bed
/usr/bin/time -v -a -o $log_dir/preparation.stderr \
	bedtools slop -b 100 -g $chrLengths_file -i $INTERVAL > $tmpbed

expTargets_file=$vc_dir/expandedTargets.bed 
/usr/bin/time -v -a -o $log_dir/preparation.stderr bedtools merge -d 50 -i $tmpbed > $expTargets_file
# Add 100bp slop to each target and merge overlapping targets.
rm -f $tmpbed

export PREFIX=j                # Prefix of file name for split targets files. This will make files like j1.bed, j2.bed, etc
rm -f ${vc_dir}/${PREFIX}*.bed  	# Remove old files
total_length=$(awk '{ l += ($3 - $2); } END { printf "%.0f", l }' < $expTargets_file)
echo $total_length
export BP=$(( 1 + total_length / NJ ))  ## Add 1 to make sure we get all
echo $BP

awk -v bp=$BP -v p=$PREFIX -v vcdir=$vc_dir 'BEGIN {print "BP:" bp} { l += ($3 - $2); j=int(l/bp); print $0 >> vcdir "/" p j ".bed"; }' $expTargets_file
time_func

echo "[Status]: Finish spliting"
echo -e "\n*************************\n"

tmp_dir=$vc_dir/tmp_vc
mkdir $tmp_dir

echo "[Process]: Start HaplotypeCaller with $NJ threads"
time_func
ls $vc_dir/j*.bed | parallel -j $NJ "/usr/bin/time -v -o {}_checkRes.stderr \
	gatk --java-options '-Xmx${parallelthreads}g -XX:ParallelGCThreads=$parallelthreads' \
	HaplotypeCaller -L {} -R $REF -I $BAM -O {}.vcf --native-pair-hmm-threads $parallelthreads \
	--tmp-dir $tmp_dir --smith-waterman FASTEST_AVAILABLE \
	2>>$vc_dir/haploCall.log" && echo "Done" || echo "Fail"

# Merging the VCFs into a final VCF
mergeVCF=${outVCF%".gz"}
grep "^#" $vc_dir/j0.bed.vcf > $mergeVCF
for C in `ls $vc_dir/j*.vcf | sort`
do
   awk '$6 >= 1 && !/^#/ { print $0 }' $C >> $mergeVCF
done && echo "Done" || echo "Fail"
time_func

# bgzip the merged VCF
bgzip -i -c -@ `nproc` $mergeVCF > $outVCF && \
	echo Done || echo Failed

if [ -f $outVCF ]; then
	echo -e "[Status]: Variant Calling done. \nOutput VCF: ${outVCF}"
	echo "Removing the split interval files and parallel VCF files" && \
	rm -f ${vc_dir}/${PREFIX}*.bed && \
	rm -f ${vc_dir}/${PREFIX}*.bed.vcf && \
	rm -f ${vc_dir}/${PREFIX}*.bed.vcf.idx && \
	rm -rf $tmp_dir && \
	echo "Generated interval files and parallel VCF files removed"
else
	echo "[Status]: Output file not found. Please check the directories."
	exit 1
fi
echo -e "\n***************************************\n"

$functions_dir/./view_stderrs.sh $vc_dir > $vc_dir/report_stderrs.txt && \
	echo -e "[Notice]: Variant calling stderr report generated\n" || echo "Fail"

echo ${alltimes[@]}
countTotalTime
echo -e "\t\t\t\t\t\t\t-----End of execution-----"
