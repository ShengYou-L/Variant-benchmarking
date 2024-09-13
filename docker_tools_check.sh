#! /bin/bash

cd /home

# Check for conda environment
echo "Checking conda environment"
conda info --envs

if [ -z `echo $CONDA_DEFAULT_ENV` ]; then			# Current activated environment
	. ~/.bashrc && conda activate && echo $CONDA_DEFAULT_ENV
fi
for i in {0..1}
do
	if [ `echo $CONDA_DEFAULT_ENV` == "base" ]; then
		echo "Checking base environment..."
		echo `which java`: `java --version`
		echo -e "`which python`: `python --version`\n"
		. ~/.bashrc && conda deactivate && conda activate python2
	elif [ `echo $CONDA_DEFAULT_ENV` == "python2" ]; then
		echo "Checking python2 environment..."
		echo `which java`: `java --version`
		echo -e "`which python`: `python2 --version`\n"
		. ~/.bashrc && conda deactivate && conda activate
	fi
done

# Check for novocraft tools
echo "Novoalign: " && echo `which novoalign`: `novoalign --version` || echo "Novoalign not found. Check if tool is installed properly and path is set correctly."
echo "Novosort: " && echo -e "`which novosort`: `novosort --version`\n" || echo "Novosort not found. Check if tool is installed properly and path is set correctly."
echo "Novoindex: " && echo -e "`which novoindex`: `novoindex --version`\n" || echo "Novoindex not found. Check if tool is installed properly and path is set correctly."

# Check for BWA-mem2
echo "BWA-mem2: " && echo -e "`which bwa-mem2`: `bwa-mem2 version`\n" || echo "BWA-mem2 not found. Check if tool is installed properly and path is set correctly"

# Check for GATK
echo "GATK: " && echo -e `which gatk`: "\n"`gatk --version`"\n" || echo "GATK not found. Check if tool is installed properly and path is set correctly"

# Check for samtools
echo "samtools: " && echo -e `which samtools`: "\n"`samtools --version | grep "samtools\|Using htslib"`"\n" || echo "samtools not found. Check if tool is installed properly and path is set correctly"

# Check for bcftools
echo "bcftools: " && echo -e `which bcftools`: "\n"`bcftools --version | grep "bcftools\|Using htslib"`"\n" || echo "bcftools not found. Check if tool is installed properly and path is set correctly"

# Check for bedtools
echo "bedtools: " && echo -e `which bedtools`: "\n"`bedtools --version`"\n" || echo "bedtools not found. Check if tool is installed properly and path is set correctly"

# Check for hap.py
. ~/.bashrc && conda deactivate && conda activate python2
echo "hap.py: " && echo -e `which hap.py`: "\n"`hap.py --version 2> >(grep "Hap.py")`"\n" || echo "hap.py not found. Check if tool is installed properly and path is set correctly"

# Check fastqc
. ~/.bashrc && conda deactivate && conda activate
echo "fastqc: "
echo `which fastqc`
echo `fastqc --version` || echo "fastqc not found. Check if fastqc is installed properly"

# Check trimmomatic
echo "Trimmomatic: "
echo `which trimmomatic`
echo `trimmomatic -version` || echo "Trimmomatic not found. Check if trimmomatic is installed properly"

# Check /usr/bin/parallel
echo "parallel: "
echo `which parallel`
echo `parallel --version` || echo "Parallel not found. Check if parallel is installed properly"

# Check for paths
echo PATHS
echo "ANT_HOME: " && echo $ANT_HOME
echo "JAVA_HOME: " && echo $JAVA_HOME
echo "PATH: " && echo $PATH
conda deactivate

# Check for R
echo "R: "
echo `which R`: `R --version`
echo -e `which Rscript`: `Rscript --version`"\n"

echo "DONE!"
