#! /bin/bash

# Now in read_mapping directory: 2 directories forward than original directory
# Input files and output files
f1=$1
f2=$2
ref=$3
MAPbam=$4
threads=$5
sample=$6
functions_dir=$7

echo "[Notice]: Directory changed to `pwd`"

# Get ref prefix by removing .fasta and show inputs
suffix=".fasta"
ref_prefix=${ref%"$suffix"}
ref_index="${ref_prefix}.bwt.2bit.64"
echo BWA-mem2 index output file: ${ref_index}

sampleNum=${sample#"HG"}

# Read group information
header=$(zcat $f1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+\+[ATGCN]+$")
rginfo="@RG\tID:rg$sampleNum\tSM:$sample\tPU:$id"_"$sm\tLB:lib-$sampleNum\tPL:ILLUMINA"
echo "Read Groups Information from FastQ header: $rginfo"

# Get record time function from functions.sh
source $functions_dir/./functions.sh
starttime=""
endtime=""
declare -a alltimes

# Run novoindex if .nix file not found
if [ ! -e ${ref_index} ]; then
	echo "[Process]: Start indexing with BWA index"
	time_func
	/usr/bin/time -v -o `pwd`/bwa_index.stderr bwa-mem2 index $ref 2>`pwd`/indexing.log && \
		echo "[Status]: Done" || echo "[Status]: Failed"
	time_func
else
	echo "[Notice]: Index files found. Skipping bwa-mem2 index."
	echo "Index files: ${ref_index}"
fi

# Run bwa
#mapsam="`pwd`/${sample}_mapped.sam"
#bwa_cmd="/usr/bin/time -v -o bwa_mapp.stderr bwa-mem2 mem -t ${threads} ${ref_prefix} ${f1} ${f2} -o ${mapsam} 2>`pwd`/bwa_mapping.log && echo Done || echo Failed"
#echo "[Process]: Start bwa with command:"
#echo $bwa_cmd
time_func
/usr/bin/time -v -o `pwd`/bwa_mapp.stderr \
	bwa-mem2 mem \
	-t ${threads} \
	${ref_prefix} \
	${f1} ${f2} \
	-R $rginfo \
	-o ${mapsam} \
	2>`pwd`/bwa_mapping.log && \
	echo "[Status]: Done" || echo "[Status]: Failed"
time_func

# Create temporary directories for processes - attempt to solve bus error (core dumped)
tempdir=`pwd`/tmpdir
mkdir $tempdir

echo "[Process]: Start converting SAM to BAM"
time_func
/usr/bin/time -v -o `pwd`/samtools_samtobam.stderr \
	gatk --java-options "-Xmx${threads}g" \
	SamFormatConverter \
	-I $mapsam -O $MAPbam \
	--TMP_DIR $tempdir \
	2>`pwd`/samtobam_conversion.log && \
	echo "[Status]: Finish converting SAM to BAM" || echo "[Status]: Failed"
time_func
# Remove SAM if BAM exists
if [ -f $MAPbam ]; then
	echo "[Status]: mapped.bam file found. Removing mapped.sam file..."
	rm -rf $mapsam && echo -e "SAM file removed successfully\n"
	rm -rf $tempdir
else
	echo "[ERROR]: Mapped.bam file not found."
	echo "Exiting..." && exit 1
fi

echo ${alltimes[@]}
countTotalTime
