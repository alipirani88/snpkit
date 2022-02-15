# SNPKIT

## Synopsis

SNPKIT is a microbial variant calling pipeline/toolkit that can be used for outbreak investigations and other clinical microbiology projects.

## Author
[Ali Pirani](https://twitter.com/alipirani88)

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Input](#input)

## Installation

The pipeline can be set up in two steps:

> 1. Clone the snpkit github directory onto your system.

```
git clone https://github.com/alipirani88/snpkit.git

```

> 2. Use snpkit/environment.yml and snpkit/environment_gubbins.yml files to set up the conda environments.

```
conda env create -f snpkit/envs/environment.yml -n snpkit
conda env create -f snpkit/envs/environment_gubbins.yml -n gubbins
```

Check installation

```
conda activate snpkit

python snpkit/snpkit.py -h
```

## Quick Start

Assuming you want to call variants for more than a few samples against a reference genome KPNIH1 and run the analysis in parallel with SLURM HPC manager. 


- You can run the first step of pipeline with option "-steps All" that will call variants for samples placed in test_readsdir against a reference genome KPNIH1.

```

python snpkit/snpkit.py \
-type PE \
-readsdir /Path-To-Your/test_readsdir/ \
-outdir /Path/test_output_core/ \
-analysis output_prefix \
-index KPNIH1 \
-steps All \
-cluster cluster \
-scheduler SLURM \
-clean

```

- The above command will run variant calling steps of the pipeline on a set of PE reads residing in test_readsdir. 
- The results will be saved in output directory test_output_core. 
- The reference genome path will be read from YAML config file.


The results of variant calling will be placed in an individual folder generated for each sample in the output directory. 

- Run the second part of the pipeline to generate SNP and Indel Matrices and various multiple sequence alignments outputs.

```
python snpkit/snpkit.py \
-type PE \
-readsdir /Path-To-Your/test_readsdir/ \
-outdir /Path/test_output_core/ \
-analysis output_prefix \
-index reference.fasta \
-steps core_All \
-cluster cluster \
-gubbins yes \
-scheduler SLURM

```

This step will gather all the variant call results from the first step, generate SNP-Indel Matrices, qc reports and core/non-core sequence alignments that can be used as an input for downstream phylogenetic analysis such as gubbins, iqtree and beast.

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

**2. config:** A high level YAML format configuration file that lets you configure your system wide runs and specify analysis parameters, path to the installed tools, data and system wide information.

- The config file stores data in KEY: VALUE pair. An example [config](https://github.com/alipirani88/snpkit/blob/master/config) file with default parameters is included with the installation folder. You can customize this config file and provide it with the -config argument or edit this config file based on your requirements. 

- Parameters for each of the tools can be customised under the 'tool_parameter' attribute of each tool in config file. 

- If you wish to run pipeline in hpc compute environment such as PBS or SLURM, change the number of nodes/cores memory reuirements based on your needs.

**3. index:** a reference genome index name as specified in a config file. 

For example; if you have set the reference genome path in config file as shown below, then the required value for command line argument -index would be -index KPNIH1

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


### For detailed information, please refer to the [wiki](https://github.com/alipirani88/snpkit/wiki) page.

