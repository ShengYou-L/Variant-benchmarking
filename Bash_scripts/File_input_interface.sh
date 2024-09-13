#! /bin/bash

# getopts Syntax:
# - getopts optstring name [args]

usage() { 
	echo -e "\nUsage: $0 "
	echo -e " -S <HG001..HG004> \t\tSample. Available samples: HG001, HG002, HG003, HG004"
	echo -e " -n <WES|WGS> \t\t\tType of sequencing reads. Either whole exome sequencing (WES) "
	echo -e "\t\t\t\tor whole genome sequencing (WGS)"
	echo -e " -m <bwa|novoalign> \t\tUsing BWA-mem2 or Novoalign for read mapping\n" 
	echo -e " Other options (optional):" 
	echo -e " -p <gatk|novosort> \t\tUsing GATK or Novosort for data preprocessing"
	echo -e " -v <gatk|freebayes> \t\tUsing GATK or Freebayes for variant calling. Default: GATK" 
	echo -e " -q <TRUE|FALSE> \t\tPerform VQSR or not. Default: TRUE"
	echo -e " -t <INT> \t\t\tNumber of threads. Default: maximum threads available using nproc."
	echo -e " -c </path/to/directory> \tContinue with previous main directory. Only use this if you are skipping any of the tools."
	echo -e " -h <INT> \t\t\tPrint this help message.\n"
	1>&2; exit 1; 
}

echo -e "\nVariant Calling Pipeline docker test 1\n"

# Main Directory for the scripts
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#echo Main directory for this script: $SCRIPT_DIR

echo "Options selected: "
# OPTSTRING below defines and stores the available options
# options with : after (eg: m:) means that they require an argument
OPTSTRING=":S:n:m:p:v:t:c:q:h"

# For each option encounters, "getopts" assigns the option to "opt" variable
# which then executes the body of the loop.
while getopts ${OPTSTRING} opt; do
	case ${opt} in								# Switch case, check case
		S)
			sample=${OPTARG}					# Can input lower or uppercase
			sample=${sample^^}					# Convert lower to uppercase
			if [ $sample != "HG001" ] && \
				[ $sample != "HG002" ] && \
				[ $sample != "HG003" ] && \
				[ $sample != "HG004" ]; then
				echo -e "[\033[0;31mERROR\e[0m] Invalid argument → -S: ${sample}" 
				echo "Choose from ' HG001 | HG002 | HG003 | HG004 ' only" && usage
			fi
			echo -e "\t-S\tSample: ${sample}"
			;;
		n)
			ngs=${OPTARG}						# Can input lower or uppercase
			ngs=${ngs^^}						# Convert lower to uppercase
			if [ $ngs != "WES" ] && [ $ngs != "WGS" ]; then
				echo -e "[\033[0;31mERROR\e[0m] Invalid argument → -n: ${ngs}"
				echo "WES or WGS only" && usage
			fi
			echo -e "\t-n\tType: ${ngs}"
			;;
		m)										
			map_tool=${OPTARG}
			map_tool=${map_tool,,}
			if [ $map_tool != "bwa" ] && [ $map_tool != "novoalign" ] && [ $map_tool != "skip" ]; then
				echo -e "[\033[0;31mERROR\e[0m] Invalid argument for -m: ${map_tool}"
				echo "Choose from ' bwa | novoalign | skip ' only" && usage
			fi
			echo -e "\t-m\tMap tool: ${map_tool}"
			;;
		p)										
			preprocess_tool=${OPTARG}
			preprocess_tool=${preprocess_tool,,}
			if [ $preprocess_tool != "gatk" ] && [ $preprocess_tool != "novosort" ] && [ $preprocess_tool != "skip" ]; then
				echo -e "[\033[0;31mERROR\e[0m] Invalid argument for -p: ${preprocess_tool}"
				echo "Choose from ' gatk | novosort | skip ' only" && usage
			fi
			echo -e "\t-p\tPreprocess tool: ${preprocess_tool}"
			;;
		v)										
			vc_tool=${OPTARG}
			vc_tool=${vc_tool,,}
			if [ $vc_tool != "gatk" ] && [ $vc_tool != "freebayes" ] && [ $vc_tool != "skip" ]; then
				echo -e "[\033[0;31mERROR\e[0m] Invalid argument for -v: ${vc_tool}"
				echo "gatk or skip only" && usage
			fi
			echo -e "\t-v\tVariant Calling tool: ${vc_tool}"
			;;
		q)										
			vf=${OPTARG}
			vf=${vf^^}
			if [ $vf != "TRUE" ] && [ $vf != "FALSE" ]; then
				echo -e "[\033[0;31mERROR\e[0m] Invalid argument for -q: ${vf}"
				echo "Choose from ' TRUE | FALSE ' only" && usage
			fi
			echo -e "\t-q\tVQSR: ${vf}"
			;;
		t)
			threads=${OPTARG}
			if ! [[ $threads =~ ^[0-9]+$ ]]; then
				echo -e "[\033[0;31mERROR\e[0m] Not a number" && usage
			fi
			echo -e "\t-t\tThreads: ${threads}"
			;;
		c)
			previous_dir=`realpath ${OPTARG}`
			if [ ! -d ${previous_dir} ]; then
				echo -e "[\033[0;31mERROR\e[0m] Directory not found: ${previous_dir}" && usage
			else
				echo -e "\t-c\tContinue from: ${previous_dir}"
			fi
			;;
		h)
			usage
			;;
		:)				# When "opt" is passed without an argument while expecting to have one
			echo -e "[\033[0;31mERROR\e[0m] Option -${OPTARG} requires an argument." 
			exit 1
			;;
		?)				# When "opt" is an invalid character
			echo -e "[\033[0;31mERROR\e[0m] Invalid option: -${OPTARG}."  
			exit 1									
			;;
		*)				# Any other input that is not in the manual
			usage
			;;
	esac
done
shift $((OPTIND-1))

# Check if sample and type of sequencing reads are given as input or not
if [ -z "${sample}" ] || [ -z "${ngs}" ] || [ -z "${map_tool}" ]; then
	echo -e "[\033[0;31mERROR\e[0m] -S, -n, and -m must be present" && usage
fi

echo -e "\nAuto-selected: "
# Assign preprocess tool if it is not mentioned in input
# Note: Novoalign goes with novosort (for now)
if [ -z "${preprocess_tool}" ]; then
	if [ $map_tool == "bwa" ]; then
		preprocess_tool="gatk"
		echo -e "\t-p\tPreprocess tool: ${preprocess_tool}"
	elif [ $map_tool == "novoalign" ]; then
		preprocess_tool="novosort"
		echo -e "\t-p\tPreprocess tool: ${preprocess_tool}"
	else
		usage
	fi
fi

# Assign variant calling tool if it is not mentioned in input
# Note: all GATK, no Freebayes yet
if [ -z "${vc_tool}" ]; then
	vc_tool="gatk"
	echo -e "\t-v\tVariant Calling tool: ${vc_tool}"
fi

# Assign variant filtering as true (default)
if [ -z ${vf} ]; then vf=TRUE && echo -e "\t-f\tVQSR: ${vf}"; fi

# Assign threads if not provided
if [ -z ${threads} ]; then
	threads=`nproc` && echo -e "\t-t\tThreads: ${threads} (maximum available)"
fi

# Get reads based on type of sequencing and sample
main_dir=$HOST_DIR && echo -e "\nMain directory: $main_dir"	# <-- NEED TO CHANGE THIS 
echo All resources will be downloaded here. Your pipeline directory and results will also be here.
F1=$main_dir/reads/${sample}/${ngs}/${sample}_${ngs}_R1_forward_paired.fastq.gz
F2=$main_dir/reads/${sample}/${ngs}/${sample}_${ngs}_R2_reverse_paired.fastq.gz
ref=$main_dir/Ref_genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

source $SCRIPT_DIR/./functions.sh

echo -e "\nChecking and listing resources..."
# Check if input reads, reference, and mapping tool are given as input or not
if [ -f "${F1}" ] && [ -f "${F2}" ] && [ -f "${ref}" ]; then
	echo -e "\nF1 reads: $F1"
	echo "F2 reads: $F2"
	echo "Reference genome: ${ref}"
else
	create_check_reads_dir $sample $ngs $F1 $F2
	create_check_ref_dir $ref
	if [ -f "${F1}" ] && [ -f "${F2}" ] && [ -f "${ref}" ]; then
		echo -e "\nF1 reads: $F1"
		echo "F2 reads: $F2"
		echo "Reference genome: ${ref}"
	elif [ ! -f "${F1}" ] || [ ! -f "${F2}" ]; then
		echo -e "[\033[0;31mERROR\e[0m] Either F1 or F2 reads not found. Please check the directories." && echo "Exit 2" && exit 2
	elif [ ! -f "${ref}" ]; then
		echo -e "[\033[0;31mERROR\e[0m] Reference genome not found. Please check the directories." && echo "Exit 2" && exit 2
	else
		echo -e "[\033[0;31mERROR\e[0m] Reads and reference are not found. Please check the directories." && echo "Exit 2" && exit 2
	fi
fi

# Assign region file 
# 	- WES use GIAB benchmark regions + Agilent exome kit regions
# 	- WGS use GIAB benchmark regions only
case "${ngs}" in
WES )
	regs=$main_dir/GIAB_benchmark/$sample/${sample}_Agilent_GIABConfidentRegions.bed
	;;
WGS )
	regs=$main_dir/GIAB_benchmark/$sample/${sample}_GRCh38_1_22_v4.2.1_benchmark.bed
	;;
esac
if [ -f $regs ]; then
	echo -e "\nRegion/Interval BED file: ${regs}\n"
else
	echo -e "\n[Notice]: Region/Interval BED file not found."
	echo "[Process]: Checking..."
	check_bed_files_and_giab_benchmark_vcf $sample $ngs
	if [ -f $regs ]; then
		echo -e "\n[Notice]: Region/Interval BED file: ${regs}\n"
	else
		echo -e "\n[ERROR]: Region/Interval BED file not found." && exit 1
	fi
fi

giab_dir=$main_dir/GIAB_benchmark
case "${ngs}" in
WES )
	if [ ! -f $giab_dir/${sample}/${sample}_${ngs}_Agilent_GIABConfidentBOTH.vcf.gz ]; then
		echo "[Process]: Start filtering the GIAB benchmark VCF using its BED file..."
		$SCRIPT_DIR/Variant_filtering_pipeline.sh $main_dir/GIAB_benchmark/$sample/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz $ngs $sample $giab_dir $SCRIPT_DIR
	else
		echo "[Notice]: Filtered GIAB benchmark VCF for $sample $ngs found."
		echo "        : $giab_dir/${sample}/${sample}_${ngs}_Agilent_GIABConfidentBOTH.vcf.gz"
		echo "        : This VCF file will be used as truth VCF file in hap.py benchmarking for this pipeline."
	fi
	;;
WGS )
	if [ ! -f $giab_dir/${sample}/${sample}_${ngs}_GIAB_benchmark_filtered.vcf.gz ]; then
		echo "[Process]: Start filtering the GIAB benchmark VCF using its BED file..."
		$SCRIPT_DIR/Variant_filtering_pipeline.sh $main_dir/GIAB_benchmark/$sample/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz $ngs $sample $giab_dir $SCRIPT_DIR
	else
		echo "[Notice]: Filtered GIAB benchmark VCF for $sample $ngs found."
		echo "        : $giab_dir/${sample}/${sample}_${ngs}_GIAB_benchmark_filtered.vcf.gz"
		echo "        : This VCF file will be used as truth VCF file in hap.py benchmarking for this pipeline."
	fi
	;;
esac

if [ $vc_tool == "gatk" ] || [ $preprocess_tool == "gatk" ]; then
	check_get_gatk_resources $main_dir
fi

#echo "---Stopped: For testing---" && exit 69
# Create output directories
source $SCRIPT_DIR/./functions.sh
create_work_dir $map_tool $vc_tool $sample $ngs $previous_dir		# Now in created main directory

# Read mapping
if [ $map_tool == "bwa" ]; then
	cd ${map_tool}_read_mapping 				# Now in read_mapping directory
	mapbam="`pwd`/${sample}_${ngs}_mapped.bam"
	echo -e "\n-----Calling bwa script for read mapping-----"			# 7 parameters
	$SCRIPT_DIR/./bwa_mem2_pipeline.sh $F1 $F2 $ref $mapbam $threads $sample $SCRIPT_DIR 1>>`pwd`/bwamem2.log 2>&1
	cd ..		 								# Back to created main directory
	if [ ! -f ${mapbam} ]; then echo "[ERROR]: Can't find mapped BAM file" && exit 1; fi	
elif [ $map_tool == "novoalign" ]; then
	cd ${map_tool}_read_mapping 				# Now in read_mapping directory
	mapbam="`pwd`/${sample}_${ngs}_mapped.bam"
	echo -e "\n-----Calling novoalign script for read mapping-----"		# 7 parameters
	$SCRIPT_DIR/./novoalign_pipeline.sh $F1 $F2 $ref $threads $mapbam $sample $SCRIPT_DIR 1>>`pwd`/novoalign.log 2>&1
	cd ..										# Back to created main directory
	if [ ! -f ${mapbam} ]; then echo "[ERROR]: Can't find mapped BAM file" && exit 1; fi	
elif [ $map_tool == "skip" ] && [ $preprocess_tool == "skip" ]; then
	echo "[Notice]: Read mapping skipped based on user request."
elif [ $map_tool == "skip" ]; then
	echo "[Notice]: Read mapping skipped based on user request."
	echo "[Process]: Looking for BAM file..."
	cd $(ls | grep 'read_mapping')
	mapbam=$(realpath $(ls *.bam | grep ".bam$") || echo "None")
	if [ ${mapbam} == "None" ]; then
		cd ../data_preprocessing
		mapbam=$(realpath $(ls *.bam | grep ".bam$") || echo "None")
		if [ ${mapbam} == "None" ]; then
			echo "[ERROR]: There is no BAM file found."
			exit 1
		else
			echo "[Notice]: BAM file found: ${mapbam}. Going to data preprocessing"
		fi
	else
		echo "[Notice]: BAM file found: ${mapbam}. Going to data preprocessing"
	fi
	cd ..
else
	echo [ERROR]: unknown error && usage
fi

# Data preprocessing
if [ $preprocess_tool == "gatk" ]; then
	cd data_preprocessing						# Now in data_preprocessing directory
	preprocessedbam="`pwd`/${sample}_${ngs}_preprocessed.bam"
	echo -e "\n-----Calling gatk script for data preprocessing-----"	# 6 parameters
	$SCRIPT_DIR/./gatk_preprocess_pipeline1.sh $mapbam $ref $sample $preprocessedbam $SCRIPT_DIR $main_dir/resource 1>>`pwd`/gatk_preprocessing.log 2>&1
	vc_inputbam=$preprocessedbam
	cd ..								# Back to created main directory
elif [ $preprocess_tool == "novosort" ]; then
	cd data_preprocessing 						# Now in data_preprocessing directory
	sortedbam="`pwd`/${sample}_${ngs}_sorted.bam"
	echo -e "\n-----Calling novosort for data preprocessing-----"
	mkdir `pwd`/tmp_dir
	source ${SCRIPT_DIR}/./functions.sh
	starttime=""
	endtime=""
	echo -e "Command: \n" >`pwd`/novosort.log
	echo "novosort \
		-t tmp_dir -c $threads -s -i --md -o $sortedbam $mapbam --stats `pwd`/novosort_stats.csv \
		&& echo 'Finish sorting' || echo 'Failed'" >>`pwd`/novosort.log
	time_func >>`pwd`/novosort.log
	echo "[Process]: Start sorting" >>`pwd`/novosort.log
	novosort \
		-t tmp_dir -c $threads -s -i --md -o $sortedbam $mapbam --stats `pwd`/novosort_stats.csv \
		1>>`pwd`/novosort.log 2>&1 && echo '[Status]: Finish sorting' >>`pwd`/novosort.log || echo '[Status]: Failed' >>`pwd`/novosort.log
	time_func >>`pwd`/novosort.log
	# Check if output exists or not
	if [ -f ${sortedbam} ]; then
		echo "[Notice]: Sorted BAM file exists. File: ${sortedbam}" >>`pwd`/novosort.log
		rm -rf `pwd`/tmp_dir
	else
		echo "[ERROR]: Sorted BAM does not exists. Exiting..." >>`pwd`/novosort.log 
		exit 1
	fi
	vc_inputbam=$sortedbam
	cd ..								# Back to created main directory
elif [ $preprocess_tool == "skip" ] && [ $vc_tool == "skip" ]; then
	echo "[Notice]: Data preprocessing skipped based on user request."
elif [ $preprocess_tool == "skip" ]; then
	echo "[Notice]: Data Preprocessing skipped based on user request."
	echo "[Process]: Looking for preprocessed/sorted BAM file..."
	cd data_preprocessing
	vc_inputbam=$(realpath $(ls | grep -E "(sorted|preprocessed).bam$") || echo "None")
	if [ ${vc_inputbam} == "None" ]; then
		echo "[ERROR]: There is no preprocessed/sorted BAM file. Please check the directories."
		echo "Exiting..." && exit 1
	else
		echo -e "[Notice]: BAM file found: ${vc_inputbam}. Going to data preprocessing\n"
	fi
	cd ..
else
	echo [ERROR]: unknown error && usage
fi

# Variant calling
if [ $vc_tool == "gatk" ]; then
	# REMEMBER TO COMMENT/REMOVE THE FOLLOWING LINE IF NOT USING DBSNP
	# ref=/export/home/shengyou/Downloads/Project1/Ref_genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
	cd ${vc_tool}_variant_calling				# Now in variant_calling directory
	output_vcf="`pwd`/${sample}_${ngs}_HaplotypeCaller_variants.vcf.gz"
	echo -e "\n-----Calling gatk script for variant calling-----"				# 7 parameters
	$SCRIPT_DIR/./VariantCalling_pipeline1.sh $ref $threads $vc_inputbam $sample $regs $output_vcf $SCRIPT_DIR 1>>`pwd`/gatk_variant_calling.log 2>&1
	unVQSR_vcf=$output_vcf
	cd ..
elif [ $vc_tool == "freebayes" ]; then
	echo "Not yet developed. Please try again next time." && exit 2
	cd ${vc_tool}_variant_calling
	echo -e "\n-----Calling freebayes script for variant calling-----"
	# $SCRIPT_DIR/gatk_preprocess_pipeline.sh $mapsam $ref $f1 $f2 $preprocessedbam $SCRIPT_DIR
	cd ..
elif [ $vc_tool == "skip" ]; then
	echo "[Notice]: Variant calling skipped based on user request."
	echo "[Process]: Looking for merged VCF file..."
	vc_dir=$(ls | grep 'variant_calling') && echo $vc_dir && cd $vc_dir
	vc_tool=${vc_dir%"_variant_calling"} && echo $vc_tool
	unVQSR_vcf=$(realpath $(ls | grep 'sorted_variants.vcf.gz$') || realpath $(ls | grep 'variants.vcf.gz$') || echo "None")
	if [ ${unVQSR_vcf} == "None" ]; then
		echo "[ERROR]: There is no merged VCF file. Please check the directories."
		echo "Exiting..." && exit 1
	else
		echo -e "[Notice]: Merged VCF file found: ${unVQSR_vcf}. Going to VQSR\n"
	fi
	cd ..
else
	echo [ERROR]: unknown error && usage
fi
if [ ! -f ${unVQSR_vcf} ]; then 
	echo "[ERROR]: Can't find VCF file." && exit 1
fi

# VQSR
if [ ${vf} == TRUE ]; then
	cd $(ls | grep 'variant_calling')				# Now in variant_calling directory
	filteredOutput_vcf="`pwd`/${sample}_${ngs}_filtered_variants.vcf.gz"
	echo -e "\n-----Calling gatk VQSR script-----"					# 8 parameters
	$SCRIPT_DIR/./VQSR_pipeline.sh $ref $unVQSR_vcf $sample $regs $filteredOutput_vcf $SCRIPT_DIR $threads $main_dir/resources 1>>`pwd`/vqsr.log 2>&1
	cd ..
	if [ -f ${filteredOutput_vcf} ]; then
		echo "[Notice]: Final VCF: ${filteredOutput_vcf}"
	else
		echo "[ERROR]: Final VCF not found. Please check the directories."
	fi
	happyInput=$filteredOutput_vcf
else
	echo "[Notice]: No VQSR done on called-variants VCF file"
	happyInput=$unVQSR_vcf
fi

# hap.py benchmarking
cd happy_output
happyOutput=`pwd`/happy_output_${ngs}_${sample}vsGIAB
echo -e "\n-----Calling hap.py script for benchmarking-----"				# 8 parameters
$SCRIPT_DIR/./hap.py_run_pipeline.sh $happyInput $sample $ngs $regs $ref $happyOutput $threads $SCRIPT_DIR $giab_dir 1>>`pwd`/happy.log 2>&1
happyGraph=`pwd`/${sample}_${ngs}_ROC_PR_graph.pdf
happyAUC=`pwd`/${sample}_${ngs}_${map_tool}_${vc_tool}.AUC.txt
if [[ $ngs == "WES" ]]; then
	Rscript $SCRIPT_DIR/happyPlot_n_AUC_pipeline.R --sample $sample --seq-type $ngs --Ymin 0.7 --graph-output $happyGraph --auc-output $happyAUC ${happyOutput}.summary.csv 1>>`pwd`/plot.log 2>&1
elif [[ $ngs == "WGS" ]]; then
	Rscript $SCRIPT_DIR/happyPlot_n_AUC_pipeline.R --sample $sample --seq-type $ngs --Ymin 0.97 --graph-output $happyGraph --auc-output $happyAUC ${happyOutput}.summary.csv 1>>`pwd`/plot.log 2>&1
fi
cd ../ && pwd
if [ -f ${happyOutput}.summary.csv ]; then
	echo -e "\nhap.py results generated at `pwd`/happy_output".
	echo "Pipeline completed. Thanks for using." && exit 0
else
	echo -e "\nIf file not found after manual checking, there might be something wrong with the pipeline."
	echo "Sorry for any inconvenience caused. Thanks for using." && exit 1
fi
