#! /bin/bash

# Input
R1=$1
R2=$2
ref=$3
threads=$4
MAPbam=$5
sample=$6
functions_dir=$7

# Filename and path for novoindex .nix file
suffix=".fasta"
ref_prefix=${ref%"$suffix"}
ref_index="${ref_prefix}.nix"
echo "Novoindex output file: ${ref_index}"

sampleNum=${sample#"HG"}

# Read group information
header=$(zcat $1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+\+[ATGCN]+$")
rginfo="@RG\tID:rg${sampleNum}\tSM:${sample}\tPU:${id}_${sm}\tLB:lib-${sampleNum}\tPL:ILLUMINA"
echo "Read Groups information from FastQ: $rginfo"

# Get record time function from functions.sh
source ${functions_dir}/./functions.sh
starttime=""
endtime=""
declare -a alltimes

# Run novoindex
if [ ! -e ${ref_index} ]; then
	echo "[Process]: Start indexing reference genome with Novoindex"
	time_func
	novoindex -k 15 -s 2 $ref_index $ref 2>`pwd`/novoindex.log && \
		echo "[Status]: DONE!" || echo "[Status]: Failed"
	time_func
else
	echo "[Notice]: Novoindex .nix file found. Skipping Novoindex."
	echo "Novoindex file: ${ref_index}"
fi

# Run novoalign
echo Command:
echo "novoalign -f $R1 $R2 -d $ref_index -c $threads --tune NOVASEQ -o BAM 1 $rginfo \> $MAPbam 2>`pwd`/novoalign.log && echo 'Finish novoalign' || echo 'Failed' "
time_func
echo "[Process]: Start read mapping with novoalign"
novoalign \
    -f $R1 $R2 \
    -d $ref_index \
    -c $threads \
    --tune NOVASEQ \
    -o BAM 1 $rginfo > $MAPbam \
    2>`pwd`/novoalign_mapping.log \
    && echo "[Status]: Finish novoalign" || echo "[Status]: Failed" 
time_func
if [ -e ${MAPbam} ]; then
	echo "[Notice]: Mapped BAM found."
else
	echo "[ERROR]: Mapped.bam file not found."
	echo "Exiting..." && exit 1
fi

echo ${alltimes[@]}
countTotalTime
