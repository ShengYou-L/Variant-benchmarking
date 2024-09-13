############## Novocraft Dockerfile #########
#Novocraft Module made for Variant Calling pipeline
#
# The Dockerfiles used are meant for Novocraft Technologies Sdn Bhd staff ONLY
# FOR INTERNAL USE ONLY. DO NOT SHARE WITH NON STAFF MEMBERS
#
# Creator:      ShengYou - 12th July 2024
#               FarhanT - 29th June 2024
#		A. Malik - 24th July 2024
# Revision:
#
# Source: (GITHUB)

FROM r-base:4.4.0
LABEL maintainer="shengyou@novocraft.com"

## Prepare conda env ##
WORKDIR /home
RUN apt-get update && apt-get install -y --no-install-recommends \
	autoconf automake \
	bc build-essential bzip2 \
	curl cmake \
	gcc git gpg-agent \
	jq \
	libbz2-dev libblas-dev libboost-all-dev libc6-dev libcurl4-gnutls-dev \
	libicu-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev \
	libtool libfontconfig1-dev libharfbuzz-dev libfribidi-dev libxml2-dev \
	libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
	make m4 \
	openjdk-17-jdk \
	perl pkg-config \
	software-properties-common \
	tabix time \
	unzip \
	vim \
	wget \
	zlib1g-dev && \
	apt -y clean && apt -y autoclean && \
	apt -y autoremove && rm -rf /var/lib/apt/lists/*

## Install required R packages ##
COPY installRpackages.R /home
RUN wget https://r-lib.r-universe.dev/bin/linux/noble/4.5/src/contrib/ragg_1.3.2.9000.tar.gz && \
	Rscript installRpackages.R && \
	rm installRpackages.R && \
	rm -rf ragg_1.3.2.9000.tar.gz

## Install miniconda for different python versions and java ##
RUN mkdir -p /home/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/miniconda3/miniconda.sh
RUN bash /home/miniconda3/miniconda.sh -b -u -p /home/miniconda3 && \
        rm -rf /home/miniconda3/miniconda.sh
ENV PATH="$PATH:/home/miniconda3/bin"
RUN conda init bash && \
        . ~/.bashrc && \
        conda activate && \
	conda install -y -c conda-forge openjdk && \
	conda install -y --solver=classic conda-forge::conda-libmamba-solver conda-forge::libmamba conda-forge::libmambapy conda-forge::libarchive

# Create python2 environment
RUN conda create -y -n python2 python=2.7.18 && \
        . ~/.bashrc && \
        conda activate python2 && conda deactivate

# Install require tools in python2 env #
RUN . ~/.bashrc && conda activate python2 && \
	pip install bx-python distribute scipy cython pandas && \
	conda install -y \
		conda-forge::openjdk=11.0.9.1 \
		conda-forge::boost-cpp \
		bioconda::pysam

# Install ANT and set environment paths (based on official website)
RUN wget https://archive.apache.org/dist/ant/binaries/apache-ant-1.9.16-bin.tar.gz && \
    	tar -xzf apache-ant-1.9.16-bin.tar.gz && \
    	rm apache-ant-1.9.16-bin.tar.gz
ENV ANT_HOME=/home/apache-ant-1.9.16
ENV PATH="$PATH:/home/apache-ant-1.9.16/bin"

# Copy hap.py from current directory into dockerfile 
WORKDIR /home
RUN mkdir -p /home/hap.py-source
COPY hap.py /home/hap.py-source
RUN mkdir -p /home/hap.py

## Require modification in order to let installation complete successfully due to error in CMake
## Solution from https://github.com/Illumina/hap.py/issues/166
## Add the following lines to hap.py/src/c++/lib/tools/Roc.cpp
##	#include <limits>
RUN sed -i '1i#include <limits>' /home/hap.py-source/src/c++/lib/tools/Roc.cpp

## Install and build hap.py with rtg tools ##
RUN . ~/.bashrc && conda activate python2 && \
	python /home/hap.py-source/install.py /home/hap.py --with-rtgtools --no-tests
ENV PATH="$PATH:/home/hap.py/bin"

## BWA-mem2 Installation ##
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
  | tar jxf -
ENV PATH="$PATH:/home/bwa-mem2-2.2.1_x64-linux"

## Samtools Installation ##
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
        tar -vxjf samtools-1.20.tar.bz2 && \
        cd samtools-1.20 && make && make install && \
	rm /home/samtools-1.20.tar.bz2 
ENV PATH="$PATH:/home/samtools-1.20"

## Bcftools Installation ##
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
        tar -vxjf bcftools-1.20.tar.bz2 && \
        cd bcftools-1.20 && make && make install && \
	rm /home/bcftools-1.20.tar.bz2
ENV PATH="$PATH:/home/bcftools-1.20"

## Bedtools Installation ##
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
        tar -zxvf bedtools-2.31.1.tar.gz && \
        cd bedtools2 && make && \
	rm /home/bedtools-2.31.1.tar.gz
ENV PATH="$PATH:/home/bedtools2/bin"

## Install novocraft ##
WORKDIR /home
COPY novocraft /home/novocraft
ENV PATH="$PATH:/home/novocraft"

## GATK Installation ##
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
RUN unzip gatk-4.5.0.0.zip && rm gatk-4.5.0.0.zip
ENV PATH="$PATH:/home/gatk-4.5.0.0"

## Install FastQC ##
RUN conda install -c bioconda fastqc

## Install Trimmomatic ##
RUN conda install -c bioconda trimmomatic

## Install parallel ##
RUN wget https://ftpmirror.gnu.org/parallel/parallel-latest.tar.bz2 && \
	tar -xjf parallel-latest.tar.bz2 && \
	rm -rf parallel-latest.tar.bz2 && \
	cd $(ls | grep 'parallel') && \
	./configure && \
	make && \
	make install

## Check all tools ##
COPY docker_tools_check.sh /home
RUN ./docker_tools_check.sh > check_tools.txt

## Import bash scripts ##
RUN mkdir -p /home/bash_scripts
COPY Bash_scripts/* /home/bash_scripts
RUN cd /home/bash_scripts && \
	chmod +x File_input_interface.sh && \
	chmod +x functions.sh && \
	chmod +x view_stderrs.sh && \
	chmod +x bwa_mem2_pipeline.sh && \
	chmod +x freebayes_with_parallel_pipeline.sh && \
	chmod +x gatk_preprocess_pipeline1.sh && \
	chmod +x VariantCalling_pipeline1.sh && \
	chmod +x VQSR_pipeline.sh && \
	chmod +x hap.py_run_pipeline.sh && \
	chmod +x novoalign_pipeline.sh && \
	chmod +x Variant_filtering_pipeline.sh

ENTRYPOINT ["/bin/bash", "/home/bash_scripts/File_input_interface.sh"]

