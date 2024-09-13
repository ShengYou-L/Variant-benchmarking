#! /bin/bash

cd /home

# Check for conda environment
conda info --envs

if [ -z `echo $CONDA_DEFAULT_ENV` ]; then			# Current activated environment
	conda activate && echo $CONDA_DEFAULT_ENV
fi
for i in {0..1}
do
	if [ `echo $CONDA_DEFAULT_ENV` == "base" ]; then
		echo "Checking base environment..."
		echo -e `which java`: `java --version`
		echo -e `which python`: `python --version`
		conda deactivate && conda activate python2
	elif [ `echo $CONDA_DEFAULT_ENV` == "python2" ]; then
		echo "Checking base environment..."
		echo -e `which java`: `java --version`
		echo -e `which python`: `python --version`
		conda deactivate && conda activate
	fi
done

# Check for novocraft tools
echo `which novoalign`: `novoalign --version` || echo "Novoalign not found. Check if tool is installed properly and path is set correctly."
echo `which novosort`: `novosort --version` || echo "Novosort not found. Check if tool is installed properly and path is set correctly."

# Check for BWA-mem2
echo `which bwa-mem2`: `bwa-mem2 version` || echo "BWA-mem2 not found. Check if tool is installed properly and path is set correctly"

# Check for GATK
echo -e `which gatk`: "\n"`gatk --version` || echo "GATK not found. Check if tool is installed properly and path is set correctly"

# Check for samtools
echo -e `which samtools`: "\n"`samtools --version | grep "samtools\|Using htslib"` || echo "samtools not found. Check if tool is installed properly and path is set correctly"

# Check for hap.py
conda deactivate && conda activate python2
echo -e `which hap.py`: "\n"`hap.py --version 2> >(grep "Hap.py")` || echo "hap.py not found. Check if tool is installed properly and path is set correctly"
echo $ANT_HOME
echo $JAVA_HOME
echo $PATH
conda deactivate

# Check for R
echo `which R`: `R --version`
echo `which Rscript`: `Rscript --version`

echo "DONE!"
