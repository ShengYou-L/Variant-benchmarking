#! /bin/bash

vcf=$1
sample=$2
ngs=$3
REGs=$4
REF=$5
happyOut=$6
thr=$7
functions_dir=$8
giab_dir=$9

# Seperate the header and records, sort the records, and merge back
# gzip file if not gzip - VQSR need block compressed format
case "$vcf" in
*.gz | *.tgz )
	echo "[Notice]: File is already gzipped."
	query_vcf=$vcf
	;;
*)
	echo "[Notice]: File is not gzipped." 
	echo "[Process]: Now gzipping ${vcf}..."
	vcfgz=${vcf}.gz
	bcftools view -Oz -o $vcfgz $vcf && echo "[Status]: Done" || echo "[Status]: Failed"
	query_vcf=$vcfgz
	;;
esac

# Activate python2 conda environment for hap.py
if [[ -z `echo $CONDA_DEFAULT_ENV` ]]; then
	conda activate python2 && echo "${CONDA_DEFAULT_ENV}: Activated python2 conda env"
elif [[ `echo $CONDA_DEFAULT_ENV` == "base" ]]; then
	conda deactivate && conda activate python2 && echo "${CONDA_DEFAULT_ENV}: Activated python2 conda env"
else
	echo "${CONDA_DEFAULT_ENV}: Activated python2 conda env"
fi

# Check if .sdf folder is present
suffix=".fasta"
ref_prefix=${REF%"$suffix"}
if [ ! -d ${ref_prefix}.sdf ]; then
	echo "[Notice]: .sdf folder for reference genome not found." 
	echo "[Process]: Generating sdf from reference genome with RTG tools"
	/usr/bin/time -v -o `pwd`/happy.stderr \
		/home/hap.py/libexec/rtg-tools-install/rtg format -o ${ref_prefix}.sdf $REF \
		2>`pwd`/sdf_ref.log && \
		echo "[Status]: Done" || echo "[Status]: Fail"
else
	echo "[Notice]: .sdf folder for reference genome found"
fi

# Threads
echo "[Notice]: Using ${thr} threads"

# Get record time function from functions.sh
source $functions_dir/./functions.sh
declare -a alltimes
starttime=""
endtime=""

# Get GIAB benchmark VCF
case "${ngs}" in
WES )
	truth_vcf=$giab_dir/$sample/${sample}_${ngs}_Agilent_GIABConfidentBOTH.vcf.gz
	exome_reg=$giab_dir/S33266436_hg38/S33266436_Regions.bed	# Get the WES regions as well
	if [ ! -f ${truth_vcf} ]; then
		echo "[Notice]: GIAB benchmark VCF not found. Is it not filtered yet? Please check the directories."
	else
		echo "[Notice]: Using GIAB benchmark VCF as truth: ${truth_vcf}"
	fi
	# Running hap.py
	time_func
	echo "[Process]: Running hap.py"
	echo "  - Truth: $truth_vcf"
	echo "  - Query: $query_vcf"
	/usr/bin/time -v -a -o `pwd`/happy.stderr \
		hap.py \
		$truth_vcf \
		$query_vcf \
		-f $REGs \
		-T $exome_reg \
		-r $REF \
		-o $happyOut \
		--roc QUAL --roc-filter LowQual --threads $thr \
		--engine=vcfeval \
		--engine-vcfeval-template ${ref_prefix}.sdf \
		2>`pwd`/happy.log \
		&& status="Done" || status="Fail"
	time_func
	;;
WGS )
	truth_vcf=$giab_dir/$sample/${sample}_${ngs}_GIAB_benchmark_filtered.vcf.gz
	if [ ! -f ${truth_vcf} ]; then
		echo "GIAB benchmark VCF not found. Is it not filtered yet? Please check the directories."
	else
		echo "Using GIAB benchmark VCF as truth: ${truth_vcf}"
	fi
	# Running hap.py
	time_func
	echo "Running hap.py"
	echo " - Truth: $truth_vcf"
	echo " - Query: $query_vcf"
	/usr/bin/time -v -a -o `pwd`/happy.stderr \
		hap.py \
		$truth_vcf \
		$query_vcf \
		-f $REGs \
		-r $REF \
		-o $happyOut \
		--roc QUAL --roc-filter LowQual --threads $thr \
		--engine=vcfeval \
		--engine-vcfeval-template ${ref_prefix}.sdf \
		2>`pwd`/happy.log \
		&& status="Done" || status="Fail"
	time_func
	;;
esac

echo "[Status]: $status"
if [ $status == "Done" ] && [ -f ${happyOut}.summary.csv ]; then
	echo "[Status]: hap.py benchmarking done successfully"
	echo "[Notice]: Summary file: ${happyOut}.summary.csv"
elif [ $status == "Fail" ]; then
	echo "[Status]: hap.py failed and summary file not found."
	echo "[Notice]: Retry after resorting the query VCF file."
	echo "[Process]: Start resorting query VCF file..."
	query_prefix=${query_vcf%".vcf.gz"}
	zcat $query_vcf | grep "^#" > $query_prefix.header
	zcat $query_vcf | grep -v "^#" > $query_prefix.record
	sort -k 1,1V -k2,2n $query_prefix.record > $query_prefix.record.sorted
	
	query_vcf_sorted=${query_prefix}_sorted.vcf.gz
	cat $query_prefix.header $query_prefix.record.sorted > ${query_vcf_sorted%".gz"}
	bgzip -c ${query_vcf_sorted%".gz"} > $query_vcf_sorted
	if [ -f $query_vcf_sorted ] && [ ! -z $query_vcf_sorted ]; then
		echo "[Status]: Sorted the VCF file"
		bcftools index $query_vcf_sorted
		echo "Removing intermediate files"
		rm -f $query_prefix.header $query_prefix.record $query_prefix.record.sorted
	else
		echo "[Status]: Something is wrong. Either vcf file not found or not sorted properly."
		exit 1
	fi
	
	case "${ngs}" in
	WES )
		truth_vcf=$giab_dir/$sample/${sample}_${ngs}_Agilent_GIABConfidentBOTH.vcf.gz
		exome_reg=$giab_dir/S33266436_hg38/S33266436_Regions.bed	# Get the WES regions as well
		# Running hap.py
		time_func
		echo "[Process]: Running hap.py"
		echo "  - Truth: $truth_vcf"
		echo "  - Query: $query_vcf_sorted"
		/usr/bin/time -v -a -o `pwd`/happy.stderr \
			hap.py \
			$truth_vcf \
			$query_vcf_sorted \
			-f $REGs \
			-T $exome_reg \
			-r $REF \
			-o $happyOut \
			--roc QUAL --roc-filter LowQual --threads $thr \
			--engine=vcfeval \
			--engine-vcfeval-template ${ref_prefix}.sdf \
			2>>`pwd`/happy.log \
			&& status="Done" || status="Fail"
		time_func
		;;
	WGS )
		truth_vcf=$giab_dir/$sample/${sample}_${ngs}_GIAB_benchmark_filtered.vcf.gz
		# Running hap.py
		time_func
		echo "[Process]: Running hap.py"
		echo "  - Truth: $truth_vcf"
		echo "  - Query: $query_vcf_sorted"
		/usr/bin/time -v -a -o `pwd`/happy.stderr \
			hap.py \
			$truth_vcf \
			$query_vcf_sorted \
			-f $REGs \
			-r $REF \
			-o $happyOut \
			--roc QUAL --roc-filter LowQual --threads $thr \
			--engine=vcfeval \
			--engine-vcfeval-template ${ref_prefix}.sdf \
			2>>`pwd`/happy.log \
			&& status="Done" || status="Fail"
		time_func
		;;
	esac
	echo "[Status]: $status"
	if [ $status == "Done" ] && [ -f ${happyOut}.summary.csv ]; then
		echo "[Status]: hap.py benchmarking done successfully"
		echo "[Notice]: Summary file: ${happyOut}.summary.csv"
	else
		echo "[ERROR]: hap.py failed and summary file not found."
		echo "Exiting..." && exit 1
	fi
fi
conda deactivate && conda activate && echo "Deactivated python2 env. Back to base environmnent"
echo -e "\t\t\t\t\t\t\t-----End of execution-----"
