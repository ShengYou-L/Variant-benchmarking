#! /bin/bash

# Record start, end time, and time used
time_func () {    
    if [ -z $starttime ]; then
        starttime=`date +%s`
        echo `date`
    else
        endtime=`date +%s`
        echo `date`
        timeused=$(( ${endtime} - ${starttime} ))
        echo Total time used: ${timeused}s
        echo -n "Time used in minutes: "; echo "scale=2; ${timeused}/60" | bc -l
        echo "************************************************"
        alltimes+=($timeused)
        starttime=""
        endtime=""
    fi
}

# Count total time used for preprocessing
countTotalTime () {
    declare -i TTime
    for i in ${alltimes[@]}
    do 
    	TTime=$(( TTime + i ))
    done
    echo Total time used: ${TTime}s
    echo -n "Total time used in minutes: "; echo "scale=2; ${TTime}/60" | bc -l
}

# Create working directory for the pipeline
create_work_dir () {
	maptool=$1
	vc_tool=$2
	sample=$3
	ngs=$4
	previous_dir=$5
	echo Current directory: $HOST_DIR			# Directory where the script is called
	dir="$HOST_DIR/Pipeline_${sample}_${ngs}_${maptool}_${vc_tool}"
	map_dir="$dir/${maptool}_read_mapping"
	dp_dir="$dir/data_preprocessing"
	vc_dir="$dir/${vc_tool}_variant_calling"
	happy_dir="$dir/happy_output"
	if [ $maptool == "skip" ]; then 
		echo "No read mapping. No creating new directories. Using ${previous_dir}"
		cd $previous_dir && pwd
	else
		echo "Creating directory: "
		echo "└───> ${dir}"
		echo "      ├──────> ${map_dir}"
		echo "      ├──────> ${dp_dir}"
		echo "      ├──────> ${vc_dir}"
		echo "      └──────> ${happy_dir}"
		mkdir $dir; cd $dir									# Now in created main directory
		mkdir ${map_dir}; mkdir ${dp_dir}; mkdir ${vc_dir}; mkdir ${happy_dir}
	fi
}

# Create directory for reads and check if reads are present. If reads are not present, then it will download the reads and preprocess them properly.
create_check_reads_dir () {
	sample=$1
	ngs=$2
	f1=$3
	f2=$4
	echo Current directory: $HOST_DIR			# Directory where the script is called
	read_dir="$HOST_DIR/reads"
	sample_dir="${read_dir}/${sample}"
	seqtype_dir="${sample_dir}/${ngs}"
	echo "Reads directory: "
	echo "└───> ${read_dir}"
	echo "      └──────> ${sample_dir}"
	echo "      		 └──────> ${seqtype_dir}"
	if [ ! -d $read_dir ]; then mkdir $read_dir; else echo "Reads directory exists. Skipped"; fi
	if [ ! -d $sample_dir ]; then mkdir $sample_dir; else echo "$sample directory exists. Skipped"; fi
	if [ ! -d $seqtype_dir ]; then mkdir $seqtype_dir; else echo "$ngs directory in $sample directory exists. Skipped"; fi
	
	# Check if reads are present
	if [ -f $f1 ] && [ -f $f2 ]; then
		echo "Trimmed $sample $ngs F1 and F2 reads found."
	else
		echo "Trimmed F1 reads not found. Assumming QC and trimming are not done."
		echo "Searching for original reads..."
		if [ -d $seqtype_dir/Ori_reads ]; then
			F1_ori=$(ls $seqtype_dir/Ori_reads 2>/dev/null | grep "R1" 2>/dev/null || echo "none")
			F2_ori=$(ls $seqtype_dir/Ori_reads 2>/dev/null | grep "R2" 2>/dev/null || echo "none")
			F1_ori=$seqtype_dir/Ori_reads/$F1_ori
			F2_ori=$seqtype_dir/Ori_reads/$F2_ori
			echo F1: $F1_ori
			echo F2: $F2_ori
			if [ $F1_ori == "none" ] && [ $F2_ori == "none" ]; then
				echo "Original $sample $ngs F1 and F2 reads not found"
				download_reads $sample $ngs $seqtype_dir/Ori_reads
				if [ $ngs == "WES" ]; then
					F1_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.wes_agilent.100x.R1.fastq.gz
					F2_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.wes_agilent.100x.R2.fastq.gz
				elif [ $ngs == "WGS" ]; then
					F1_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.pcr-free.50x.R1.fastq.gz
					F2_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.pcr-free.50x.R2.fastq.gz
				fi
				echo "Start QC and trimming..."
				qc_trimming $F1_ori $F2_ori $sample $ngs
			elif [ -f $F1_ori ] && [ -f $F2_ori ]; then
				echo "Original $sample $ngs F1 reads found: $F1_ori"
				echo "original $sample $ngs F2 reads found: $F2_ori"
				echo "Start QC and trimming..."
				qc_trimming $F1_ori $F2_ori $sample $ngs
			fi
		else
			echo "Original reads directory not found. Creating directory..."
			mkdir $seqtype_dir/Ori_reads && echo "Directory created."
			echo "Original $sample $ngs F1 and F2 reads not found"
			if [ $ngs == "WES" ]; then
				F1_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.wes_agilent.100x.R1.fastq.gz
				F2_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.wes_agilent.100x.R2.fastq.gz
			elif [ $ngs == "WGS" ]; then
				F1_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.pcr-free.50x.R1.fastq.gz
				F2_ori=$seqtype_dir/Ori_reads/${sample}.novaseq.pcr-free.50x.R2.fastq.gz
			fi
			download_reads $sample $ngs $seqtype_dir/Ori_reads
			echo "Start QC and trimming..."
			qc_trimming $F1_ori $F2_ori $sample $ngs
		fi
	fi
}

# Download reads if not present
download_reads () {
	sample=$1
	ngs=$2
	target_dir=$3
	echo "[Downloading]: ${sample} ${ngs} F1 & F2 reads"
	case "${ngs}" in
	WES )
		wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wes_agilent/100x/${sample}.novaseq.wes_agilent.100x.R1.fastq.gz -P $target_dir && \
			echo "F1 reads downloaded"
		wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wes_agilent/100x/${sample}.novaseq.wes_agilent.100x.R2.fastq.gz -P $target_dir && \
			echo "F2 reads downloaded"
		echo "Original reads at $target_dir"
		;;
	WGS )
		wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/50x/${sample}.novaseq.pcr-free.50x.R1.fastq.gz -P $target_dir && \
			echo "F1 reads downloaded"
		wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/50x/${sample}.novaseq.pcr-free.50x.R2.fastq.gz -P $target_dir && \
			echo "F2 reads downloaded"
		echo "Original reads at $target_dir"
		;;
	esac
}

# QC and Trimming with FastQC and Trimmomatic for originally downloaded reads
# Parameters are all the same because the QC results for all of the reads are almost the same.
qc_trimming () {
	f1_ori=$1
	f2_ori=$2
	sample=$3
	ngs=$4
	ori_dir=$(dirname $f1_ori)
	out_dir=${ori_dir%"/Ori_reads"}
	starttime=""
	endtime=""
	time_func
	echo "[Process]: FastQC on $sample $ngs original reads"
	fastqc -t 2 -q -o $ori_dir $f1_ori $f2_ori && \
		echo "1st QC Done. Report generated here: $ori_dir" || echo "Fail"
	echo "[Process]: Trimmomatic trimming adapters using Trimmomatic's adapter file: TruSeq3-PE-2.fa"
	echo "Parameters: ILLUMINACLIP:/home/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:80"
	trimmomatic PE -threads `nproc` \
		$f1_ori $f2_ori \
		$out_dir/${sample}_${ngs}_R1_forward_paired.fastq.gz \
		$out_dir/${sample}_${ngs}_R1_forward_unpaired.fastq.gz \
		$out_dir/${sample}_${ngs}_R2_reverse_paired.fastq.gz \
		$out_dir/${sample}_${ngs}_R2_reverse_unpaired.fastq.gz \
		ILLUMINACLIP:/home/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
		SLIDINGWINDOW:5:20 MINLEN:80 \
		1>>$ori_dir/trimmomatic_${sample}_${ngs}.log \
		2>$ori_dir/trim.err && \
		echo "Trimmomatic Done" || echo "Fail"
	echo "[Process]: FastQC on $sample $ngs trimmed reads"
	fastqc -t 2 -q -o $out_dir \
		$out_dir/${sample}_${ngs}_R1_forward_paired.fastq.gz \
		$out_dir/${sample}_${ngs}_R2_reverse_paired.fastq.gz && \
		echo "2nd QC Done. Report generated here: $out_dir" || echo "Fail"
	time_func
	echo "[Notice]: Removing original reads: $f1_ori & $f2_ori"
	rm -rf $f1_ori $f2_ori && echo "[Notice]: Removed original reads."
}

create_check_ref_dir () {
	ref=$1
	ref_dir=$(dirname $ref)
	if [ -f $ref ]; then
		echo "Reference genome fasta file for hg38 from GIAB found. $ref"
	elif [ -f ${ref}.gz ]; then
		echo "Gzipped Reference genome fasta file for hg38 from GIAB found. $ref"
		echo "Better to unzip the file. hap.py need unzipped fasta file."
		echo "[Process]: Unzipping reference genome..."
		gzip -d ${ref}.gz
	else
		echo "Reference genome fasta file for hg38 from GIAB not found."
		echo "[Downloading]: GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz -P $ref_dir && \
			echo "Reference genome fasta file for hg38 from GIAB downloaded at $ref_dir"
		gzip -d ${ref}.gz
	fi
	# Check if .fai file is present
	if [ -f ${ref}.fai ]; then
		echo ".fai file found"
	else
		echo ".fai file not found."
		echo "[Process]: Indexing with samtools..."
		samtools faidx $ref && echo "Done" || echo "Fail"
	fi
}

check_get_gatk_resources () {
	main_dir=$1
	# Check resources for BQSR and VQSR
	if [ ! -d $main_dir/resources ]; then mkdir $main_dir/resources; fi
	rm $main_dir/resources/current_resources.txt $main_dir/resources/req_resources.txt $main_dir/resources/need_download_res.txt
	ls $main_dir/resources | grep -v '.tbi\|.fai\|.idx' > $main_dir/resources/current_resources.txt
	echo -e "\nList of resources required by GATK found: "
	cat $main_dir/resources/current_resources.txt
	
	# All required resources
	resources=( 1000G_omni2.5.hg38.vcf.gz 
		1000G_phase1.snps.high_confidence.hg38.vcf.gz 
		Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz 
		hapmap_3.3.hg38.vcf.gz 
		Homo_sapiens_assembly38.dbsnp138.vcf 
		Homo_sapiens_assembly38.known_indels.vcf.gz 
		Mills_and_1000G_gold_standard.indels.hg38.vcf.gz )
	for res in ${resources[@]}; do echo $res >> $main_dir/resources/req_resources.txt; done
	
	# Find difference between two files using comm and output the ones that are not present in the current_resources.txt
	comm -13 $main_dir/resources/current_resources.txt $main_dir/resources/req_resources.txt > $main_dir/resources/need_download_res.txt
	if [ ! -s $main_dir/resources/need_download_res.txt ]; then 
		echo "***All resources are present. No need to download anything.***"
	else
		echo -e "\nNeed to download resources: " && cat $main_dir/resources/need_download_res.txt
		for line in $(cat $main_dir/resources/need_download_res.txt); do
			case "${line}" in
				1000G_phase1.snps.high_confidence.hg38.vcf.gz )
					echo "[Downloading]: 1000G high confidence snps"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz -P $main_dir/resources && \
						echo "1000G high confidence snps downloaded" || echo "1000G high confidence snps download fail"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi -P $main_dir/resources && \
						echo "1000G high confidence snps tbi downloaded" || echo "1000G high confidence snps tbi download fail"
					;;
				1000G_omni2.5.hg38.vcf.gz )
					echo "[Downloading]: 1000G omni2.5"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz -P $main_dir/resources && \
						echo "1000G omni2.5 downloaded" || echo "1000G omni2.5 download fail"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi -P $main_dir/resources && \
						echo "1000G omni2.5 tbi downloaded" || echo "1000G omni2.5 tbi download fail"
					;;
				Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz )
					echo "[Downloading]: Axiom Exome Plus genotypes..."
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz -P $main_dir/resources && \
						echo "Axiom Exome Plus genotypes downloaded" || echo "Axiom Exome Plus genotypes download fail"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi -P $main_dir/resources && \
						echo "Axiom Exome Plus genotypes tbi downloaded" || echo "Axiom Exome Plus genotypes tbi download fail"
					;;
				hapmap_3.3.hg38.vcf.gz )
					echo "[Downloading]: hapmap_3.3"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz -P $main_dir/resources && \
						echo "hapmap_3.3 downloaded" || echo "hapmap_3.3 download fail"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi -P $main_dir/resources && \
						echo "hapmap_3.3 tbi downloaded" || echo "hapmap_3.3 tbi download fail"
					;;
				Homo_sapiens_assembly38.dbsnp138.vcf )
					echo "[Downloading]: dbsnp138"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf -P $main_dir/resources && \
						echo "dbsnp138 downloaded" || echo "dbsnp138 download fail"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx -P $main_dir/resources && \
						echo "dbsnp138 idx downloaded" || echo "dbsnp138 idx download fail"
					;;
				Homo_sapiens_assembly38.known_indels.vcf.gz )
					echo "[Downloading]: known_indels"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz -P $main_dir/resources && \
						echo "known_indels downloaded" || echo "known_indels download fail"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi -P $main_dir/resources && \
						echo "known_indels tbi downloaded" || echo "known_indels tbi download fail"
					;;
				Mills_and_1000G_gold_standard.indels.hg38.vcf.gz )
					echo "[Downloading]: 1000G gold standard indels"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -P $main_dir/resources && \
						echo "1000G gold standard indels downloaded" || echo "1000G gold standard indels download fail"
					wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi -P $main_dir/resources && \
						echo "1000G gold standard indels tbi downloaded" || echo "1000G gold standard indels tbi download fail"
					;;
			esac
		done
	fi
}

check_bed_files_and_giab_benchmark_vcf () {
	sample=$1
	ngs=$2
	giab_dir=$HOST_DIR/GIAB_benchmark
	echo "GIAB resources directory: "
	echo "└───> ${giab_dir}"
	echo "      └──────> $giab_dir/$sample"
	if [ ! -d $giab_dir ]; then 
		mkdir $giab_dir && echo Created directory: $giab_dir
	else 
		echo GIAB directory present.
	fi
	if [ ! -d $giab_dir/$sample ]; then 
		mkdir $giab_dir/$sample && echo Created directory: $giab_dir/$sample
	else 
		echo $sample directory in GIAB directory present.
	fi
	case "$ngs" in
	WES )
		if [ ! -f $giab_dir/${sample}/${sample}_Agilent_GIABConfidentRegions.bed ]; then
			echo "No merged interval file found."
			echo "Searching for GIAB benchmark interval BED file..."
			if [ ! -f $giab_dir/${sample}/${sample}_GRCh38_1_22_v4.2.1_benchmark.bed ]; then
				echo "GIAB benchmark interval BED file for hg38-$sample not found"
				echo "[Downloading]: ${sample}_GRCh38_1_22_v4.2.1_benchmark.bed"
				download_giab_benchmark $sample $giab_dir/$sample
			fi
			echo "[Process]: Joining Agilent WES kit regions and GIAB benchmark regions with bedtools"
			# Joins two bed files: WES regions and GIAB confidence -> Joined_confidenceRegions.bed
			# The Agilent WES kit regions can't be install with wget or link from their official website so need to put at the GitHub repository. The below -b need to change.
			/usr/bin/time -v -o $giab_dir/$sample/mergeRegions.stderr \
				bedtools intersect -loj \
				-a $giab_dir/${sample}/${sample}_GRCh38_1_22_v4.2.1_benchmark.bed \
				-b $giab_dir/S33266436_hg38/S33266436_Regions.bed \
				> $giab_dir/${sample}/${sample}_Agilent_GIABConfidentRegions.bed\
				&& echo "Done" || echo "Fail"
		fi
		;;
	WGS )
		if [ ! -f $giab_dir/${sample}/${sample}_GRCh38_1_22_v4.2.1_benchmark.bed ]; then
			echo "GIAB benchmark interval BED file for hg38-${sample} not found."
			echo "[Downloading]: ${sample}_GRCh38_1_22_v4.2.1_benchmark.bed"
			download_giab_benchmark $sample $giab_dir/$sample
		fi
		;;
	esac
}

download_giab_benchmark () {
	sample=$1
	targetDir=$2
	echo "[Downloading]: "
	case "$sample" in
	HG001 )
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed -P $targetDir && \
			echo "BED file downloaded" || echo "Fail to download BED file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -P $targetDir && \
			echo "VCF file downloaded" || echo "Fail to download VCF file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi -P $targetDir && \
			echo "VCF index file downloaded" || echo "Fail to download VCF index file"
		;;
	HG002 )
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -P $targetDir && \
			echo "BED file downloaded" || echo "Fail to download BED file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -P $targetDir && \
			echo "VCF file downloaded" || echo "Fail to download VCF file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi -P $targetDir && \
			echo "VCF index file downloaded" || echo "Fail to download VCF index file"
		mv $targetDir/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed $targetDir/HG002_GRCh38_1_22_v4.2.1_benchmark.bed && \
			echo "Renamed to HG002_GRCh38_1_22_v4.2.1_benchmark.bed"
		;;
	HG003 )
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -P $targetDir && \
			echo "BED file downloaded" || echo "Fail to download BED file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -P $targetDir && \
			echo "VCF file downloaded" || echo "Fail to download VCF file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi -P $targetDir && \
			echo "VCF index file downloaded" || echo "Fail to download VCF index file"
		mv $targetDir/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed $targetDir/HG003_GRCh38_1_22_v4.2.1_benchmark.bed && \
			echo "Renamed to HG003_GRCh38_1_22_v4.2.1_benchmark.bed"
		;;
	HG004 )
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -P $targetDir && \
			echo "BED file downloaded" || echo "Fail to download BED file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -P $targetDir && \
			echo "VCF file downloaded" || echo "Fail to download VCF file"
		wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi -P $targetDir && \
			echo "VCF index file downloaded" || echo "Fail to download VCF index file"
		mv $targetDir/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed $targetDir/HG004_GRCh38_1_22_v4.2.1_benchmark.bed && \
			echo "Renamed to HG004_GRCh38_1_22_v4.2.1_benchmark.bed"
		;;
	esac
}
