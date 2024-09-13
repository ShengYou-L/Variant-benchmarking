#! /bin/bash

# This script has one function:
# 1.
# 	Split input VCFs into SNPs and INDELs, intersecting each with high confidence regions BED file,
#	and merge the high confidence SNPs and INDELs together into a VCF file.

# Note: 
# -	Run this script outside the output directory on the local pc.
# - This script is modified from the script for confidence region in which this aims to filter 
#	variants using the GIAB and Exome target kit high confidence regions created from the original 
#	script.

usage () {
	echo "Variant filtering"
	echo "./Variant_filtering.sh [Input VCF] [Regions <WES|WGS>] [sample name] [/path/to/giab/dir]"
}

INP="`realpath ${1}`"
inpdir=$(realpath $(dirname $INP))
ngs=$2
sample=$3
giab_dir=$4
functions_dir=$5

case "$ngs" in
WES )
	# Generated BED file that intersects GIAB benchmark regions and Agilent's WES kit regions
	regions=$giab_dir/${sample}/${sample}_Agilent_GIABConfidentRegions.bed
	OUT=$giab_dir/${sample}/${sample}_${ngs}_Agilent_GIABConfidentBOTH.vcf.gz
	;;
WGS )
	# GIAB BED file that intersects GIAB benchmark regions
	regions=$giab_dir/$sample/${sample}_GRCh38_1_22_v4.2.1_benchmark.bed
	OUT=$giab_dir/${sample}/${sample}_${ngs}_GIAB_benchmark_filtered.vcf.gz
	;;
*)
	usage
	exit 1
esac

# Generating stats using bcftools for input file
bcftools stats $INP > ${INP}_stats.txt

# Intermediate files for SNPs and INDELs extracted from freebayes VCF file
SNPs_tmp="${inpdir}/SNPs.vcf"
INDELs_tmp="${inpdir}/INDELs.vcf"

# Intermediate files for SNPs and INDELs intersected with High Confidence regions BED file
SNPs_intersected="${inpdir}/HighConfidentSNPs.vcf"
INDELs_intersected="${inpdir}/HighConfidentINDELs.vcf"

if [ ! -d $inpdir/vf_logs ]; then mkdir $inpdir/vf_logs; fi

usage
# Get record time function from functions.sh
source $functions_dir/./functions.sh
declare -a alltimes
starttime=""
endtime=""

time_func
echo "[Process]: Split Input VCF"
# Split freebayes vcf into SNP.vcf and INDEL.vcf
/usr/bin/time -v -o $inpdir/vf_logs/vcfFiltering.stderr \
	gatk SplitVcfs \
	-I $INP \
	-SNP_OUTPUT $SNPs_tmp \
	-INDEL_OUTPUT $INDELs_tmp \
	-STRICT false \
	2>$inpdir/vf_logs/variantfiltering_process.log && \
	echo "[Status]: Done...1" || echo "[Status]: Fail...1"

echo "[Process]: Intersect SNPs with high confidence regions"
# Intersect SNPs with high confidence regions to get Joined_ConfidenceSNPs.vcf
/usr/bin/time -v -a -o $inpdir/vf_logs/vcfFiltering.stderr \
	bedtools intersect -header -u \
	-a $SNPs_tmp \
	-b $regions \
	>$SNPs_intersected \
	2>>$inpdir/vf_logs/variantfiltering_process.log && \
	echo "[Status]: Done...2" || echo "[Status]: Fail...2"

echo "[Process]: Intersect INDELs with high confidence regions"
# Intersect INDELs with high confidence regions to get Joined_ConfidenceINDELs.vcf
/usr/bin/time -v -a -o $inpdir/vf_logs/vcfFiltering.stderr \
	bedtools intersect -header -u \
	-a $INDELs_tmp \
	-b $regions \
	>$INDELs_intersected \
	2>>$inpdir/vf_logs/variantfiltering_process.log && \
	echo "[Status]: Done...3" || echo "[Status]: Fail...3"

echo "[Process]: Merge intersected SNPs and INDELs to output VCF"
# Join and sort both Joined_Confidence(SNPs/INDELs).vcf to get a final HighConfidentVariants.vcf
/usr/bin/time -v -a -o $inpdir/vf_logs/vcfFiltering.stderr \
	gatk SortVcf \
	-I $SNPs_intersected \
	-I $INDELs_intersected \
	-O $OUT \
	2>>$inpdir/vf_logs/variantfiltering_process.log && \
	echo "[Status]: Done...4" || echo "[Status]: Fail...4"

# Generating stats using bcftools for output file
bcftools stats $OUT > ${OUT}_stats.txt

time_func

if [ -f $OUT ]; then
	echo "[Status]: GIAB benchmark VCF file filtered using its interval BED file"
	echo "Filtered GIAB benchmark VCF: ${OUT}"
	echo "Removing intermediate files"; \
	rm -f $SNPs_tmp $INDELs_tmp ${SNPs_tmp}.idx ${INDELs_tmp}.idx \
		$SNPs_intersected $INDELs_intersected && \
		echo "Finish removing" || echo "Fail to remove files"
else
	echo "[ERROR]: Filtered GIAB benchmark VCF not found. Exiting..." && exit 1
fi

