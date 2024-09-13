#! /bin/bash

MAPbam=$1
ref=$2
sample=$3
sprefix="HG"
sample_no=${sample#"$sprefix"}
outputfile=$4
functions_dir=$5
res_dir=$6

# Output/Intermediate files
markedDups="`pwd`/marked_dups.bam"
baseRecalTable="`pwd`/BaseRecal_out_data.table"

# Create temporary directories for processes - attempt to solve bus error (core dumped)
tempdir=`pwd`/tmpdir
mkdir $tempdir
# Set Java maximum mem to 50% of total RAM available
javamem=$(( $(free -g | grep "Mem" | awk '{print $2}') / 2 ))

# Create directory for log
log_dir=`pwd`/Preprocess_job_log
if [ -d $log_dir ]; then echo "[Notice]: Making log dir: $log_dir" && mkdir $log_dir; else echo "[Notice]: Directory already present. Storing log files here."; fi

echo "[Notice]: Directory changed to `pwd`"		# Directory: data_preprocessing

# Check if $MAPbam is mapped.bam or marked_dups.bam
if [[ $MAPbam =~ mapped ]]; then
	markdup="start"
elif [[ $MAPbam =~ marked_dup ]]; then
	echo "[Notice]: Marked Duplicates BAM file found. Skipping mark duplicates..."
	markedDups=$MAPbam
	markdup="skip"
fi

# Get record time function from functions.sh
source $functions_dir/./functions.sh
declare -a alltimes
starttime=""
endtime=""

# Preprocessing
if [ -f $MAPbam ] && [ $markdup == "start" ]; then
	echo -e "\n[Notice]: Using mapped BAM for marking duplicates."
	echo "[Process]: Start marking duplicates"
	time_func
	/usr/bin/time -v -o $log_dir/markDups.stderr \
		gatk --java-options "-Xmx${javamem}g" \
		MarkDuplicatesSpark \
		-I $MAPbam \
		-O $markedDups \
		--tmp-dir $tempdir \
		2>`pwd`/markDups.log \
		&& echo "[Status]: Finish marking duplicates" || echo "[Status]: Failed"
	time_func
	if [ -f $markedDups ]; then 
		echo "[Status]: Marked duplicates BAM file created successfully"
		echo "Removing mapped BAM file..."
		rm -rf $MAPbam && echo -e "Mapped BAM file removed successfully.\n"
		rm -rf $tempdir/*
	else
		echo "[ERROR]: Can't find marked duplicates BAM file. Exiting..."
		exit 1
	fi
elif [ -f $markedDups ] && [ $markdup == "skip" ]; then
	echo "[Notice]: Skipped mark duplicates"
else
	echo "[ERROR]: Can't find markDuplicates.bam file. Exiting..."
	exit 1
fi

if [ -f $markedDups ]; then
	echo "[Process]: Start generating recalibration table"
	time_func
	/usr/bin/time -v -o $log_dir/gatk_BQSR.stderr \
		gatk --java-options "-Xmx${javamem}g" \
		BaseRecalibratorSpark \
		-I $markedDups -R $ref \
		--known-sites $res_dir/Homo_sapiens_assembly38.dbsnp138.vcf \
		--known-sites $res_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
		-O $baseRecalTable \
		--tmp-dir $tempdir \
		2>`pwd`/baseRecalibrator.log \
		&& echo "[Status]: Finish generating recalibration table" || echo "[Status]: Failed"
	time_func
fi

if [ -f $baseRecalTable ]; then
	echo "[Process]: Start Applying BQSR"
	time_func
	/usr/bin/time -v -a -o $log_dir/gatk_BQSR.stderr \
		gatk --java-options "-Xmx${javamem}g" \
		ApplyBQSRSpark \
		-R $ref -I $markedDups \
		--bqsr-recal-file $baseRecalTable \
		-O $outputfile \
		--tmp-dir $tempdir \
		2>`pwd`/applyBQSR.log \
		&& echo "[Status]: DONE!" || echo "[Status]: Failed"
	time_func
fi
if [ -f $outputfile ]; then 
	echo "[Status]: Preprocessed BAM file created successfully. Data preprocessing completed."
	echo "Removing marked duplicates BAM file and mapped BAM file..."
	rm -rf $markedDups $markedDups.bai $markedDups.sbi && echo "Marked duplicates BAM file removed successfully."
else
	echo "[ERROR]: Can't find preprocessed BAM file. Exiting..."
	exit 1
fi
rm -rf $tempdir

echo ${alltimes[@]}
countTotalTime
