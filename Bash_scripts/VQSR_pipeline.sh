#! /bin/bash

# ./VQSR_pipeline.sh <Reference_genome.fasta> <Input.vcf> <sample(HG001)> <Region/Interval.bed> <Filtered_output.vcf> </path/to/functions.sh/directory>

REF=$1
INP=$2
sample=$3
REGs=$4
OUT=$5
functions_dir=$6
threads=$7
res_dir=$8

# Intermediate files
TRANCHES1="`pwd`/VRoutput_SNP.tranches"
RECAL1="`pwd`/VRoutput_SNP.recal"
TRANCHES2="`pwd`/VRoutput_INDEL.tranches"
RECAL2="`pwd`/VRoutput_INDEL.recal"
MID="`pwd`/appliedVQSR_SNP.vcf.gz"

# Set Java maximum mem to 50% of total RAM available
javamem=$(( $(free -g | grep "Mem" | awk '{print $2}') / 2 ))

# Create directory for log
log_dir=`pwd`/VQSR_job_log
if [ ! -d $log_dir ]; then echo "Creating $log_dir for log files" && mkdir $log_dir; fi

echo Directory changed to `pwd` 			# Directory: gatk_variant_calling

# Get record time function from functions.sh
source $functions_dir/./functions.sh
declare -a alltimes
starttime=""
endtime=""

sts="Done"
# gzip file if not gzip - VQSR need block compressed format
case "$INP" in
*.gz | *.tgz )
	echo "[Notice]: VCF file is already gzipped."
	;;
*)
	echo "[Notice]: VCF file is not gzipped." && echo "[Process]: Now gzipping ${INP}..."
	bgzip -c -@ $threads $INP > ${INP}.gz && echo "[Status]: Gzipped" || sts="Failed"
	INP=${INP}.gz
	;;
esac
if [ ! -f ${INP}.tbi ]; then
	echo -e "[Notice]: No .tbi index file found. \n[Process]: Indexing..."
	tabix -p vcf $INP && echo "[Status]: Indexed" || sts="Failed"
	tbi_size=$(stat -c%s ${INP}.tbi)
	echo "[Notice]: ${INP}.tbi file size = $tbi_size"
	if [[ $tbi_size -lt 100 ]]; then sts="Failed" && echo "[WARNING]: File size too small"; else echo "[Status]: Valid file size."; fi
fi

if [ $sts == "Failed" ]; then
	echo "[Status]: ${INP} file can't be indexed"
	echo "[Notice]: Retrying after resorting." 
	echo "[Process]: Start resorting this VCF.gz file..."
	inp_prefix=${INP%".vcf.gz"}
	zcat $INP | grep "^#" > $inp_prefix.header
	zcat $INP | grep -v "^#" > $inp_prefix.record
	sort -k 1,1V -k2,2n $inp_prefix.record > $inp_prefix.record.sorted
	
	INP_sorted=${inp_prefix}_sorted.vcf.gz
	cat $inp_prefix.header $inp_prefix.record.sorted > ${INP_sorted%".gz"}	
	bgzip -c ${INP_sorted%".gz"} > $INP_sorted
	if [ -f $INP_sorted ]; then
		echo "[Status]: Sorted the VCF file"
		echo "[Process]: Indexing to .tbi and .gzi..."
		tabix -p vcf $INP_sorted && \
		echo "[Status]: Done .tbi" || echo "[Status]: Failed .tbi"
		bcftools index $INP_sorted && \
		echo "[Status]: Done .gzi" || echo "[Status]: Failed .gzi"
		echo "Removing intermediate files"
		rm -f $inp_prefix.header $inp_prefix.record $inp_prefix.record.sorted
		INP=$INP_sorted
	else
		echo "[ERROR]: Something is wrong. Either vcf file not found, improperly formatted, or not sorted properly."
		exit 1
	fi
fi

# SNP - VariantRecalibrator
time_func
echo "[Process]: Start VariantRecalibrator - SNP"
/usr/bin/time -v -o $log_dir/VQSR.stderr \
	gatk --java-options "-Xmx${javamem}g" \
	VariantRecalibrator \
	-R $REF \
	-V $INP \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $res_dir/hapmap_3.3.hg38.vcf.gz \
	--resource:omni,known=false,training=true,truth=false,prior=12.0 $res_dir/1000G_omni2.5.hg38.vcf.gz \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 $res_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $res_dir/Homo_sapiens_assembly38.dbsnp138.vcf \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	-L $REGs \
	-O $RECAL1 \
	--tranches-file $TRANCHES1 \
	--rscript-file `pwd`/VRoutput_SNP.plots.R \
	2>`pwd`/VariantRecal_SNP.log && echo "[Status]: DONE!" || echo "[Status]: Failed"
time_func

# SNP - ApplyVQSR
time_func
echo "[Process]: Start Applying VQSR - SNP"
/usr/bin/time -v -a -o $log_dir/VQSR.stderr \
	gatk --java-options "-Xmx${javamem}g" \
	ApplyVQSR \
	-R $REF \
	-V $INP \
	-L $REGs \
	-O $MID \
	--truth-sensitivity-filter-level 99.0 \
	--tranches-file $TRANCHES1 \
	--recal-file $RECAL1 \
	-mode SNP \
	2>`pwd`/ApplyVQSR_SNP.log && echo "[Status]: DONE!" || echo "[Status]: Failed"
time_func

# INDEL - VariantRecalibrator
time_func
echo "[Process]: Start VariantRecalibrator - INDEL"
/usr/bin/time -v -a -o $log_dir/VQSR.stderr \
	gatk --java-options "-Xmx${javamem}g" \
	VariantRecalibrator \
	-R $REF \
	-V $MID \
	--resource:mills,known=false,training=true,truth=true,prior=12.0 $res_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--resource:axiom,known=false,training=true,truth=false,prior=10.0 $res_dir/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $res_dir/Homo_sapiens_assembly38.dbsnp138.vcf \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode INDEL \
	-L $REGs \
	-O $RECAL2 \
	--tranches-file $TRANCHES2 \
	--rscript-file `pwd`/VRoutput_INDEL.plots.R \
	2>`pwd`/VariantRecal_INDEL.log && echo "[Status]: DONE!" || echo "[Status]: Failed"
time_func

# INDEL - ApplyVQSR
time_func
echo "[Process]: Start Applying VQSR - INDEL"
/usr/bin/time -v -a -o $log_dir/VQSR.stderr \
	gatk --java-options "-Xmx${javamem}g" \
	ApplyVQSR \
	-R $REF \
	-V $MID \
	-L $REGs \
	-O $OUT \
	--truth-sensitivity-filter-level 99.0 \
	--tranches-file $TRANCHES2 \
	--recal-file $RECAL2 \
	-mode INDEL \
	2>`pwd`/ApplyVQSR_INDEL.log && echo "[Status]: DONE!" || echo "[Status]: Failed"
time_func

# Check if filtered VCF file exists or not
if [ -f $OUT ]; then 
	echo "[Status]: Variant called and filtered. Variant calling and filtering completed."
	echo "Removing intermediate files..."
	rm -f $TRANCHES1 $RECAL1 $RECAL1.idx $TRANCHES2 $RECAL2 $RECAL2.idx $MID $MID.tbi && \
		echo "Finish removing intermediate files"
else
	echo "[ERROR]: Filtered variants file not found. Exiting..."
	exit 1
fi

echo ${alltimes[@]}
countTotalTime
