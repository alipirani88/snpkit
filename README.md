# Variant Calling and core SNP diagnostics Pipeline

## Synopsis

This is a highly customisable, automated variant detection pipeline that can be easily deployed for infectious disease outbreak investigations and other clinical microbiology projects.  

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Input](#input)

## Installation

The pipeline can be set up in two easy steps:

> 1. Clone the github directory onto your system.

```
git clone https://github.com/alipirani88/variant_calling_pipeline.git

```

> 2. Use variant_calling_pipeline/environment.yml and variant_calling_pipeline/environment_gubbins.yml files to create conda environment.

Create two new environments - varcall and varcall_gubbins
```
conda env create -f variant_calling_pipeline/environment.yml -n varcall
conda env create -f variant_calling_pipeline/environment_gubbins.yml -n varcall_gubbins
```

Check installation

```
conda activate varcall

python variant_calling_pipeline/variant_call.py -h
```

## Quick Start

Assuming you want to generate core snps for more than a few hundred samples and run the analysis in parallel on cluster (time and memory efficient). The default resources used for parallel jobs are: 

```
nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00

OR

--nodes=1 --ntasks=1 --cpus-per-task=1 --mem=5g --time=125:00:00
```

See option resources in the scheduler section of the [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file. Detailed information is in the section [Customizing config file](#customizing-the-config-file)

- Run variant calling step (All) on a set of PE reads with default parameters

```
python variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index reference.fasta -steps All -cluster cluster -scheduler SLURM -clean

```

The above command will run the variant calling (step 1) pipeline on a set of PE reads residing in test_readsdir. The results will be saved in the output directory test_output_core. The config file contains options for some frequently used reference genomes. To know which reference genomes are included in the config file, look up the [config]() file or check the help menu of the pipeline.

The results of variant calling will be placed in an individual folder generated for each sample in the output directory. A log file for each sample will be generated and can be found in each sample folder inside the output directory. A single log file of this step will be generated in the main output directory. For more information on log file prefix and convention, please refer to the [log](#log) section below.

- Generate core SNPs/Matrices from variant calling results 

```
python variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index reference.fasta -steps core_All -cluster cluster -gubbins yes -scheduler SLURM

```

## Input

The pipeline requires three main inputs - 

**1. readsdir:** Place your Illumina SE/PE reads in a folder and give path to this folder with -readsdir argument. Apart from the standard Miseq/Hiseq fastq naming convention (R1_001_final.fastq.gz), other acceptable fastq extensions are: 

```

- R1.fastq.gz/_R1.fastq.gz, 
- 1_combine.fastq.gz, 
- 1_sequence.fastq.gz, 
- _forward.fastq.gz, 
- _1.fastq.gz/.1.fastq.gz.

```

**2. config:** A high level easy to write YAML format configuration file that lets you configure your system wide runs and specify analysis parameters, path to the installed tools, data and system wide information.

- This config file will contain High level information such as locations of installed programs like GATK, cores and memory usage for running on HPC compute cluster, path to a reference genome, various parameters used by different tools. These settings will apply across multiple runs and samples. 

- The config file stores data in KEY: VALUE pair. 

- An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters is included with the installation folder. You can customize this config file and provide it with the -config argument or edit this config file based on your requirements. 

- Parameters for each of the tools can be customised under the 'tool_parameter' attribute of each tool in config file. 

- If you wish to run pipeline in hpc compute environment such as PBS or SLURM, change the number of nodes/cores memory reuirements based on your needs else the pipeline will run with default settings.


**3. index:** a reference genome index name as specified in a config file. For example; if you have set the reference genome path in config file as shown below, then the required value for command line argument -index would be -index KPNIH1

```
[KPNIH1]
# path to the reference genome fasta file.
Ref_Path: /nfs/esnitkin/bin_group/variant_calling_bin/reference/KPNIH1/
# Name of reference genome fasta file.
Ref_Name: KPNIH1.fasta
```

Here, Ref_Name is the reference genome fasta file located in Ref_Path. Similarly, if you want to use a different version of KPNIH reference genome, you can create a new section in your config file with a different index name.

```
[KPNIH1_new]
# path to the reference genome fasta file.
Ref_Path: /nfs/esnitkin/bin_group/variant_calling_bin/reference/KPNIH1_new/
# Name of reference genome fasta file.
Ref_Name: KPNIH1_new.fasta
```

THe pipeline also requires Phaster results of your reference genome to mask phage region. To enable this place the phaster results files in the reference genome folder.



