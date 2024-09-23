# Variant-benchmarking-pipelines

## Benchmark Variant Calls from different combinations of Read Mapping and Variant Calling tools to GIAB High Confidence Variant Calls

This project compares the GIAB high confidence variant calls to the variant calls from different pipelines using hg38 reference genome and hg38 HG001 to HG004 sample reads. Illumina/hap.py was used for the benchmarking process.

If you are interested in reproducing the result, then you can either try run the pipelines using [docker]() or you can download this repository and build the docker image using the Dockerfile and resources provided.

## Pipelines used in this project

1. **Novoalign + GATK**
	- Read Mapping: [Novoalign](https://www.novocraft.com/products/novoalign/)
	- Preprocessing and Sorting: [Novosort](https://www.novocraft.com/products/novosort/)
	- Variant Calling: GATK's HaplotypeCaller
	- Variant Filtering: GATK's VQSR 
2. **BWA-mem2 + GATK**
	- Read Mapping: [BWA-mem2](https://github.com/bwa-mem2/bwa-mem2)
	- Preprocessing and Sorting: [GATK's pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
	- Variant Calling: GATK's HaplotypeCaller
	- Variant Filtering: GATK's VQSR

*Note: Variant Calling and Filtering using GATK is done by following [GATK's best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels)

## Results
### ROC/PR curves and AUC values
#### WGS

![Figure 1: Graph of ROC/PC curves for WGS HG001 SNP_SEL](/results/WGS/HG001/happy_HG001_WGS_SNP_SEL.svg)
![Figure 2: Graph of ROC/PC curves for WGS HG001 INDEL_SEL](/results/WGS/HG001/happy_HG001_WGS_INDEL_SEL.svg)

![Figure 3: Graph of ROC/PC curves for WGS HG002 SNP_SEL](/results/WGS/HG002/happy_HG002_WGS_SNP_SEL_a2.svg)
![Figure 4: Graph of ROC/PC curves for WGS HG002 INDEL_SEL](/results/WGS/HG002/happy_HG002_WGS_INDEL_SEL_a1.svg)

![Figure 5: Graph of ROC/PC curves for WGS HG003 SNP_SEL](/results/WGS/HG003/happy_HG003_WGS_SNP_SEL_a2.svg)
![Figure 6: Graph of ROC/PC curves for WGS HG003 INDEL_SEL](/results/WGS/HG003/happy_HG003_WGS_INDEL_SEL_a1.svg)

![Figure 7: Graph of ROC/PC curves for WGS HG004 SNP_SEL](/results/WGS/HG004/happy_HG004_WGS_SNP_SEL_a2.svg)
![Figure 8: Graph of ROC/PC curves for WGS HG004 INDEL_SEL](/results/WGS/HG004/happy_HG004_WGS_INDEL_SEL_a1.svg)

**Table 1: AUC values for Selectively-filtered SNPs (SNP_SEL) in WGS HG001**
|SNP_SEL|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novoalign+GATK|0.937|0.934|0.935|0.934|
|BWA-mem2+GATK|0.908|0.901|0.901|0.9|

**Table 2: AUC values for Selectively-filtered INDELs (INDEL_SEL) in WGS HG001**
|INDEL_SEL|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novoalign+GATK|0.968|0.967|0.966|0.966|
|BWA-mem2+GATK|0.952|0.958|0.961|0.96|

#### WES

![Figure 1: Graph of ROC/PC curves for WES sample SNP_SEL](/results/WES/HG001/happy_HG001_WES_SNP_SEL_a1.svg)
![Figure 2: Graph of ROC/PC curves for WES sample INDEL_SEL](/results/WES/HG001/happy_HG001_WES_INDEL_SEL_a1.svg)

![Figure 3: Graph of ROC/PC curves for WES sample SNP_SEL](/results/WES/HG002/happy_HG002_WES_SNP_SEL_a1.svg)
![Figure 4: Graph of ROC/PC curves for WES sample INDEL_SEL](/results/WES/HG002/happy_HG002_WES_INDEL_SEL_a1.svg)

![Figure 5: Graph of ROC/PC curves for WES sample SNP_SEL](/results/WES/HG003/happy_HG003_WES_SNP_SEL_a1.svg)
![Figure 6: Graph of ROC/PC curves for WES sample INDEL_SEL](/results/WES/HG003/happy_HG003_WES_INDEL_SEL_a1.svg)

![Figure 7: Graph of ROC/PC curves for WES sample SNP_SEL](/results/WES/HG004/happy_HG004_WES_SNP_SEL_a2.svg)
![Figure 8: Graph of ROC/PC curves for WES sample INDEL_SEL](/results/WES/HG004/happy_HG004_WES_INDEL_SEL_a1.svg)

**Table 1: AUC values for Selectively-filtered SNPs (SNP_SEL) in WES HG001**
|SNP_SEL|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novoalign+GATK|0.918|0.916|0.916|0.913|
|BWA-mem2+GATK|0.913|0.906|0.912|0.910|

**Table 2: AUC values for Selectively-filtered INDELs (INDEL_SEL) in WES HG001**
|INDEL_SEL|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novoalign+GATK|0.810|0.796|0.842|0.847|
|BWA-mem2+GATK|0.805|0.828|0.840|0.827|

>[!IMPORTANT]
>Conclusion: The AUC values, and overall Precision and Recall values using the Novoalign+GATK pipeline are higher than that in the BWA-mem2+GATK pipeline especially in WGS reads.

### Time used

#### WGS

**Read mapping**:

|Tool|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novoalign|20.96 hrs|27.23 hrs|24.04 hrs|24.94 hrs|
|BWA-mem2|4.35 hrs|8.01 hrs|6.64 hrs|6.69 hrs|


**Data preprocessing**:

|Tool|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novosort|4.52 hrs|0.87 hrs|0.88 hrs|0.86 hrs|
|GATK(MarkDuplicate Spark + BQSR Spark)|47.37 hrs|24.27 hrs|18.94 hrs|14.09 hrs|

**Variant Calling**: Both pipelines used GATK's HaplotypeCaller with GNU's parallel for variant calling. The time taken on average is about **1.5 hrs**.

**Variant Filtering**: Both pipelines used GATK's VQSR for variant filtering. The time taken on average is about **1 hr**.

**Benchmarking with hap.py**: The time taken is between **35 to 45 minutes**.

#### WES

**Read mapping**:

|Tool|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novoalign|277.36 mins|304.25 mins|269.13 mins|273.86 mins|
|BWA-mem2|30.2 mins|53.23 mins|69.55 mins|57.15 mins|


**Data preprocessing**:

|Tool|HG001|HG002|HG003|HG004|
|---|---|---|---|---|
|Novosort|7.08 mins|6.61 mins|7.05 mins|6.6 mins|
|GATK(MarkDuplicate Spark + BQSR Spark)|45.74 mins|72.48 mins|86.6 mins|52.62 mins|

**Variant Calling**: Both pipelines used GATK's HaplotypeCaller with GNU's parallel for variant calling. The time taken on average is about **32.5 minutes**.

**Variant Filtering**: Both pipelines also used GATK's VQSR for variant filtering. The time taken on average is also about **32.4 minutes**.

**Benchmarking with hap.py**: The time taken is between **6 to 7 minutes**.


Note: All analysis were done on a server with 32 cores and 64 threads.

## Installation

You can docker pull the [docker image for this pipeline]() from docker hub.

    docker pull

Note: If you pull from docker hub, you might want to download the "GIAB_benchmark" directory from this repository that contains the WES Agilent Kit regions/intervals file which is not easily downloadable with a link.

Alternatively, you can download this repository with the Dockerfile and build the docker yourself.

    docker build -t <image_name>:<image_tag> /path/to/Dockerfile

To build the docker yourself with the Dockerfile, you will have to download [Novoalign+Novosort](https://www.novocraft.com/support/download/) and [hap.py](https://github.com/Illumina/hap.py) into the directory containing the Dockerfile.

After that, if you try to build the Docker image, it might fail due to a bootstrap error when installing hap.py. 

So in this repository, there is a directory "hap.py_fixes" which contains the boost_subset_1_58_0.tar.bz2 and make_dependencies.sh. You can replace these two files in the hap.py/external directory. 

Or you can apply the fixes to the boost library yourself from [here](https://github.com/boostorg/build/commit/48e9017139dd94446633480661e5447c7e0d8b1b#diff-cfcd688bf771e8e7d55826fd1f3916623e47543ddfad1d3428c81437046bcec6).

*Note: 
1. The docker image is 10.3Gb in memory.
2. You can run the pipeline in the downloaded directory. If you want to run in another directory, please move the "GIAB_benchmark" directory to your directory as well.

## Novoalign License

To run Novoalign, a trial license is required. Please request one from Novocraft Technologies (https://www.novocraft.com/contact-us/).

## Instructions
**TL;DR**

Command template:
```
docker run -v $(pwd):$(pwd) -v </path/to/novoalign.lic>:/home/novocraft/novoalign.lic -e HOST_DIR=$(pwd) variant-benchmarking:latest -S <HG001|..|HG004> -n <WES|WGS> -m <novoalign|bwa>
```

**Details**

Basic command structure:
```
docker run -v <your/volume>:<your/volume> -v </path/to/novoalign.lic>:/home/novocraft/novoalign.lic -e HOST_DIR=<path/to/desired/directory> variant-benchmarking:latest -S <HG001|..|HG004> -n <WES|WGS> -m <novoalign|bwa>
```

The -e HOST_DIR parameter sets the docker environment to point to your directory and will be used to as a main diretory for all the resources to be downloaded. 
The path should be the same or be under the directory in the first option in the -v mount.

Both the -v <your/volume>:<your/volume> and -e HOST_DIR=<path/to/your/dir> are required if you are wondering.

Example:
```
docker run -v /home:/home -v /export/home/usr1/novocraft/novoalign.lic:/home/novocraft/novoalign.lic -e HOST_DIR=$(pwd) variant-benchmarking:latest -S hg004 -n WES -m novoalign
```

Note: If the directory that contains the GIAB_benchmark is somewhere else, then please change the path in -e HOST_DIR to your path (-e HOST_DIR=/path/to/your/dir).

### Options
```
Required options:
 -S <HG001..HG004>  Sample. Available samples: HG001, HG002, HG003, HG004
 -n <WES|WGS>       Type of sequencing reads.
 -m <bwa|novoalign> Read mapping tool. BWA-mem2 or Novoalign.

Optional options:
 -p <gatk|novosort> Data preprocessing tool.
 -v <gatk>          Variant Calling tool.
 -q <TRUE|FALSE>    Perform variant filtering with VQSR or not.
 -t <INT>           Number of threads.
 -c </path/to/dir>  Continue with previous main directory.
		    Ex: /home/Downloads/Pipeline_HG002_WES_novoalign_gatk
		    You can use this if one of the processes failed.
                    You can specify <skip> for -m, -p, -v options if they are completed.
                    The pipeline with find the file required for the next process.
 -h                 Print this help message.
```

You can run ``` docker run variant-benchmarking:latest ``` to view the parameters required and available.

## Directories and resources
### Directories
Below are the directories that will be created when you run the pipeline (some names are examples):
```
your_current_working_directory (HOST_DIR)
	│
	├───> Pipeline_{sample}_{seqtype}_{maptool}_{variantcallingtool} (ex: Pipeline_HG001_WGS_novoalign_gatk)
	│     ├──────> {maptool}_read_mapping (ex: novoalign_read_mapping)
	│     ├──────> data_preprocessing
	│     ├──────> {variantcallingtool}_variant_calling (ex: gatk_variant_calling)
	│     └──────> happy_output
	│
	├───> reads 
	│     ├──────> HG001
	│     │        ├──────> WES
	│     │        └──────> WGS
	│     ├──────> HG002
	│     ├──────> HG003
	│     └──────> HG004
	│
	├───> GIAB_benchmark
	│     ├──────> HG001
	│     ├──────> HG002
	│     ├──────> HG003
	│     └──────> HG004
	│
	└───> resources
```
*Note:
1. If you run one sample, let's say HG001 WGS, only the directories related to this will be created like HG001 and WGS directories.
2. It is best that you run the pipeline in a designated directory because all required resources (see Resources) will be downloaded here.

### Resources

The resources will be downloaded automatically when you run the pipeline. 
Reads will be downloaded and pre-processed automatically using Trimmomatic. FastQC reports of before and after the trimming will be generated for you. Only related resources will be downloaded. 

- reads are downloaded from [brain-genomics-public](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/sequencing/fastq?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false)
- Reference genome used is downloaded from GIAB [GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/)
- GIAB benchmark VCF and BED files are downloaded from [GIAB](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/)
- resources for GATK BQSR and VQSR are downloaded from [genomics-public-data](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false)

Let's say you run a pipeline with HG003 WES, the WES HG003 reads, and the HG003 GIAB benchmark VCF and BED files will be downloaded. The reference genome will be downloaded during the first time you run the pipeline (only once, unless you delete the reference). If GATK is part of the pipeline, then the resources required by GATK will be downloaded accordingly.

## All Tools Used in this Project
List of tools and version installed in the docker:
1. Novoalign v4.04.01
2. Novosort v3.02.03
3. GATK v4.5.0
4. BWA-mem2 v2.2.1
5. Samtools v1.20
6. Bedtools v2.31.1
7. hap.py v0.3.15
8. R v4.4.0
9. Bcftools v1.4.1
10. GNU parallel 20240722

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.
