# Variant Calling and core SNP diagnostics Pipeline

## Synopsis


This pipeline calls variants on PE/SE reads provided in a directory and generates core SNP/consensus fasta files that can be used to to build a phylogeny or as an input for Gubbins/Beast analysis

## Contents

- [Installation](#installation)
- [Input](#input)
- [Steps](#steps)
- [Command line options](#command-line-options)
- [Run pipeline on Compute cluster](#run-pipeline-on-compute-cluster)
- [Quick Start](#quick-start)
- [Output Files](#output-files)
- [Customizing Config file](#customizing-config-file)
- [Log](#log)
- [Bonus Ducks](#bonus-ducks)

## Installation

Pending. 

Ignore this if you are in snitkin lab. The dependencies are already installed in lab bin_group folder: 

/nfs/esnitkin/bin_group/variant_calling_bin/. 

Use the python version installed in:

/nfs/esnitkin/bin_group/anaconda2/bin/python


## Input


Input is a directory(-readsdir) containing SE/PE reads and a config file where all the configuration settings for the pipeline are set. This config file settings will be used universally on all samples available in readsdir. An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters are included in the pipeline folder. You can customize this config file and provide it with the -config argument. 

Detailed information in section [Customizing Config file](#customizing-config-file)

Note: Apart from standard Miseq/Hiseq fastq naming extensions (R1_001_final.fastq.gz), other acceptable fastq extensions are: R1.fastq.gz/_R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, _1.fastq.gz/.1.fastq.gz. 


## Steps


![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/pipeline.png)

There are three main steps to generate core SNPs which can be provided with -steps argument and should be run in sequential order.

**1. Variant Calling:** This step will run all the standard variant calling steps on sample read files residing in input reads directory. 
The possible options are:

Option ***All*** :  This will run all the variant calling steps starting from cleaning the reads to variant calling and generating intermediate files that will be used in later steps. 

Option ***clean,align,post-align,varcall,filter,stats*** :  This option will run variant calling steps starting from cleaning to variant calling, variant filtering and generating stats.

You can modify this argument to run a part of the pipeline. For example, if you want to restart the pipeline from alignment step, you can do it by supplying the following option: "align,post-align,varcall,filter,stats" which will skip the Trimmomatic cleaning part.


Note: The order of variant calling steps needs to be sequential to avoid any errors. While skipping any of the steps, make sure the results for skipped steps are already present in your output folder.


**2. Preparing files for Core SNP extraction and diagnostics purposes:**


Option ***core_prep*** : Run this step before running the last core steps. This will prepare all the intermediate data required for the last core step i.e generating the consensus fasta file out of core snps. This can be a time consuming step depending on how closely related the samples are to the reference genome.  


**3. Generate Core SNP consensus and data matrix for diagnostics plots:**

Option ***core*** : This step will generate core SNP consensus fasta file and a consensus fasta of only core variant positions. Various data matrices will generated at this step that can be used later for diagnostics purposes. 

**4. Generate report for the pipeline:**

Option ***report*** : This step will generate final core results directory and various reports that will summarize the alignment and core SNP results. All the final results will be saved into date_time_core_results directory inside the output folder. 

```

2018_01_13_13_18_03_core_results
├── core_snp_consensus
└── data_matrix

```

core_snp_consensus directory contains the core consensus fasta and vcf files.
data_matrix contains all matrices and reports generated during report step.  

As the name suggests, data_matrix will contain various matrices that can be queried or plotted for further diagnosing the variant call results. Alternatively, you can run a R script provided inside the data_matrix folder to generate the plots. 

Require: ggplot2 and heatmap.3

```

module load R

Rscript generate_diagnostics_plots.R 

```

| Extension | Description |
| --------- | ----------- |
| . barplot.pdf |  Distribution of filter-pass variant positions(not just core) in each sample. colors represents the filter criteria that caused them to get filtered out in that particular sample.|
| . barplot_DP.pdf |  |
| . temp_Only_filtered_positions_for_closely_matrix_FQ.pdf |  |
| . DP_position_analysis.pdf |  |
| . temp_Only_filtered_positions_for_closely_matrix_DP.pdf |  |


## Command line options

```

usage: variant_call.py [-h] -type TYPE -readsdir DIR -outdir OUTPUT_FOLDER
                       -index INDEX [-steps STEPS] -analysis ANALYSIS_NAME
                       [-config CONFIG] [-suffix SUFFIX]
                       [-filenames FILENAMES] [-cluster CLUSTER]

Variant Calling pipeline for Illumina PE/SE data.

Required arguments:

  -type         Type of reads: SE or PE
  -readsdir     Path to Sequencing Reads Data directory. Requires full/absolute path.
  -outdir       Output Folder Path ending with output directory name to save the results. Requires full/absolute path.
  -index        Reference Index Name. Most Frequently used reference genomes index options: KPNIH1 | MRSA_USA_300 | MRSA_USA_100 | CDIFF_630 | paris
  -steps        Variant Calling Steps in sequential order.
                1.   All: This will run all the steps starting from cleaning the reads to variant calling;
                2.   clean,align,post-align,varcall,filter,stats : This will also run all steps starting from cleaning to variant calling.
                3.   coverage_depth_stats: Run Only Depth of Coverage Stats module after cleaning and read mapping steps
                4.   core_prep: Run this step before running the core steps. This will prepare the data required for generating core SNPs
                5.   core: extract core snps and generate diagnostics plot data matrices to explore filtered snps.
  -analysis     Unique analysis name that will be used as prefix to saving results and log files.

Optional arguments:
  -config       Path to Config file, Make sure to check config settings before running pipeline
  -suffix       Fastq reads suffix such as fastq, fastq.gz, fq.gz, fq; Default: fastq.gz
  -filenames    fastq filenames with one single-end filename per line. if the type is set to PE, it will detect the second paired-end filename with the suffix from first filename.
  -cluster      Run variant calling pipeline in one of the four modes. Default: local. The possible modes are: cluster/parallel-cluster/parallel-local/local

```

## Run pipeline on Compute cluster


The variant calling can be run in parallel on each sample using the -cluster argument. Set the pbs resources such as resources(-l), email(-M), queue(-q), flux_account(-A) and notification(-m) under [scheduler] section in config file. For more details, check out UMICH [flux](http://arc-ts.umich.edu/systems-and-services/flux/) website.


Possible options for -cluster option(Supported system: pbs):

***local*** : This is the default option. When this option is set, the pipeline will analyze each sample one after the another. This will not make use of multiple cores present in your system or clusters. Use this option for few number of samples or for testing purposes.

***parallel-local*** :  This option will run the pipeline and analyze samples in parallel but on a local system. This is preferred for less than 20 sample size input and multiple core local system.

***cluster*** : This option will run the pipeline on a single cluster. This option is similar to local but will rather run on a cluster node. It will not make use of multiple cores present on a cluster. Use this option for few number of samples or for testing purposes or for few number of large size samples which requires multiple cores to analyze individual samples.

***parallel-cluster*** : The pipeline is optimized for this option. When this option is set, the pipeline will run everything on a cluster making use of all the cores available. This option will also submit individual jobs for each sample seperately and whenever needed. When you set this option for each step in the pipeline, make sure all the jobs submitted by the pipeline at each step is completed before proceeding to another step. You can check the status of the job with: 

```
qstat -u USERNAME  
```

## Quick Start

Assuming you want to generate core snps for more than a few hundred samples and run the analysis in parallel on cluster(Time and memory efficient). The default pbs resources used for parallel jobs are: 

```
nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00
```
See option resources in scheduler section of [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file. Detailed information in section [Customizing Config file](#customizing-config-file)

- Run variant calling step (All) on a set of PE reads with default parameters

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps All -cluster parallel-cluster

```

The above command will run variant calling (step 1) pipeline on a set of PE reads residing in test_readsdir. The results will be saved in output directory test_output_core. The config file contains options for some frequently used reference genome. To know which reference genomes are included in config file, look up the [config]() file or check the help menu of the pipeline.

The results of variant calling will be placed in an individual folder generated for each sample in output directory. A log file for each sample will be generated and can be found in each sample folder inside the out directory. A single log file of this step will be generated in main output directory. For more information on log file prefix and convention, please refer [log](#log) section below.

- Run core_prep step to generate files for core SNP calling.

Run this steps to generate various intermediate files that will be used for generating core SNPs.

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core_prep -cluster parallel-cluster

```

- Run core step to generate final core SNP consensus fasta files.

Since this step compares multiple files simultaneously and involves multiple I/O operations, It is recommended to provide higher memory compute resources. 

example:

```
nodes=1:ppn=4,mem=47000mb,walltime=24:00:00
```

Replace the resources option in scheduler section of config file with the above line before running the command.

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core -cluster parallel-cluster

```

## Output Files

***Variant Calling***: Each sample folder under output directory contains standard variant calling outputs such as clean reads, aligned BAM files, filtered non-core vcf and various stats file.

Pending...

| Extension | Description |
| --------- | ----------- |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |




***Core Variants***: The final core step results will be moved to \*_core_results directory under output directory. There are two results folder, core_snp_consensus and data_matrix. The core SNP vcf and consensus fasta files are stored in core_snp_consensus while the data matrices for useful for variant diagnostics can be found in data_matrix.

Pending...

| Extension | Description |
| --------- | ----------- |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |


## Bonus Ducks

- Run only a part of variant calling step

OPTIONAL: In case, you want to rerun the above analysis with different variant call/filter parameters (i.e skip the read cleaning, alignment and post-alignment steps), you can do it as follows:

NOTE: OPTIONAL
```

python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps varcall,filter,stats -cluster parallel-cluster

```

- Run variant calling step on selected samples provided in filenames.txt. filenames.txt should contain fastq filenames with one single-end filename per line.

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps varcall,filter,stats -cluster parallel-cluster -filenames filenames.txt
```

- Run core_prep and core steps on selected files

If you decide to exclude some samples from core snp analysis step (step 3), there is a way to do that without running the entire pipeline. 


```

python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core -cluster parallel-cluster -filenames filenames_custom

```

## Customizing Config file:

By default, the pipeline uses config file that comes with the pipeline. Make sure to edit this config file or copy it to your local system, edit it and provide path of this edited config file with -config argument.

```

cp variant_calling_pipeline/config /Path-to-local/config_edit

```

The pipeline implements customisable variant calling configurations using config file. Config file can be customised to use your choice of aligner and variant caller by changing two parameters under the section [pipeline]
Currently, The pipeline supports BWA aligner(mem algorithm) for aligning reads to the reference genome and samtools for variant calling.

```
# Set which tools to use in pipeline:
[pipeline]
# Options for Aligner:bwa / smalt / bowtie
aligner: bwa
# Options for variant_caller:  gatkhaplotypecaller /samtools
variant_caller: samtools

```

Make sure you have downloaded all the dependencies for the pipeline in a folder path provided in binbase option of [bin_path] section.

```
# Set bin folder path. Please make sure all the executables are placed in bin folder. Also make sure the path for individual tools are correct.
[bin_path]
binbase: /nfs/esnitkin/bin_group/variant_calling_bin/
```

NOTE: Add the required perl libraries(such as in the case of vcftools) PERL5LIB environment variable. For flux users, you can do that by loading perl-modules

```
module load perl-modules
```

If you wish to run the jobs on clsuter, make sure you cahnge the necessary schedular parameters in scheduler section shown below: for more information, visit [flux](http://arc-ts.umich.edu/systems-and-services/flux/) homepage.

```

[scheduler]
resources: nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00
email: username@umich.edu
queue: XXX
flux_account: XXX
notification: a

```

Every tool has its own *_bin option where you can set the folder name in which the tool resides. For example, in the below Trimmomatic section example, the Trimmomatic tool resides in /Trimmomatic/ folder that is set with trimmomatic_bin option which in itself resides in /nfs/esnitkin/bin_group/variant_calling_bin/ folder that was set in binbase option above.

```
[Trimmomatic]
trimmomatic_bin: /Trimmomatic/
adaptor_filepath: adapters/TruSeq3-Nextera_PE_combined.fa
seed_mismatches: 2
palindrome_clipthreshold: 30
simple_clipthreshold: 10
minadapterlength: 8
keep_both_reads: true
window_size: 4
window_size_quality: 20
minlength: 40
headcrop_length: 0
colon: :
targetlength: 125
crop_length: 40
f_p: forward_paired.fq.gz
f_up: forward_unpaired.fq.gz
r_p: reverse_paired.fq.gz
r_up: reverse_unpaired.fq.gz
```

Parameters for each tools can be customised under the 'tool_parameter' attribute of each tool in config file.


For example, to change the minadapterlength parameter of Trimmomatic from 8 to 10, replace minadapterlength of 8 with suppose 10 and restart the pipeline.

## Log:

The pipeline generates a log file following the naming convention: yyyy_mm_dd_hrs_mins_secs_analysisname.log.txt and tracks each event/command. The log file sections follow standard [Python logging conventions](https://docs.python.org/2/howto/logging.html): 

***INFO*** to print STDOUT messages; 

***DEBUG*** to print commands ran by pipeline, 

***ERROR*** to print STDERR messages and 

***EXCEPTION*** to print an exception that occured while the pipeline was running.





<!---
#### The variant calling pipeline runs sequentially as follows:
***

>1. Pre-Processing Raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
>2. Read Alignment using [BWA](http://bio-bwa.sourceforge.net/)
>3. Post-Alignment steps using [SAMTOOLS](http://samtools.sourceforge.net/), [GATK](https://software.broadinstitute.org/gatk/), [PICARD](https://broadinstitute.github.io/picard/), [Bedtools](http://bedtools.readthedocs.io/en/latest/) etc
>4. Variant Calling using [SAMTOOLS](http://samtools.sourceforge.net/)
>5. Variant Filtering and generating Consensus using [GATK](https://software.broadinstitute.org/gatk/), [Bedtools](http://bedtools.readthedocs.io/en/latest/), [vcftools](http://vcftools.sourceforge.net/), in-house scripts() etc.


**Usage:**
***

```
python pipeline.py [-h] -PE1 path-to-forward-PE-read -PE2 path-to-reverse-PE-read -o path-to-OUTPUT_FOLDER -analysis ANALYSIS_NAME -index INDEX_NAME_as_per_config_file -config path-to-config-file
```
**Note:**
***

- ***Input file format and extension***: Either .fastq or .fastq.gz
- ***Output Directory***: Pipeline creates output folder by output directory name mentioned at the end of the path. e.g: -o /any-path-followed-by/output_directory-name/
- ***ANALYSIS_NAME***: Name by which pipeline saves the results inside output directory. e.g: -analysis first_analysis
- ***INDEX***: Reference Index Name as mentioned under the section 'Reference Genome to be used for pipeline'. 
  e.g: In config file; under the section [KPNIH1], mention the attributes REF_NAME and REF_PATH for reference fasta filename and path to the reference fasta file resp. The section name KPNIH1 is required by the INDEX argument.
- config: Path to your customized config file. Make sure the section names are similar to the default config file.


**Variant Calling Output**:
***

The pipeline generates various output files from different tools at different steps. The most notable ones are:
- ***Clean reads***: *.fq.gz files from trimmomatic.

- ***Alignment files***: analysisname_aln.sam and analysisname_aln.bam from BWA, analysisname_aln_marked.bam from GATK MarkDuplicates, and finally a sorted BAM from marked bam file analysisname_aln_sort.bam. Also including *bai index files.
- ***Bed file***: analysisname_unmapped.bed and analysisname_unmapped.bed_positions with positions that were unmapped. Bedcoverage file analysisname_.bedcov

- ***VCF file***: Various vcf files are generated removing different types of variants at different steps.
>1. ***analysisname_aln_mpileup_raw.vcf***: The raw variant calls without any variant filtering
>2. ***analysisname_aln_mpileup_raw.vcf_5bp_indel_removed.vcf.gz***: variants that proximate to an indel by 5 bp are filtered out
>3. ***analysisname_filter2_gatk.vcf***: variants that does not pass the GATK variant filter parameters are filtered out(parameters can be changed in config file) + variants that are proximate to an indel by 5 bp are filtered out
>4. ***analysisname_final.vcf_no_proximate_snp.vcf.gz***: variants that does not pass the GATK variant filter parameters are filtered out(parameters can be changed in config file) + variants that are proximate to an indel by 5 bp are filtered out + variants that are proximate to each other by 5 bp 

- ***Statistics Reports***:
***
>1. ***analysisname_alignment_stats***: Alignment stats file generated using SAMTOOLS flagstat.
>2. ***analysisname_vcf_stats***: vcf stats(raw) generated using vcftools
>3. ***analysisname_depth_of_coverage***: Depth of Coverage generated using GATK Depth of Coverage.
>4. ***analysisname_markduplicates_metrics***: Mark Duplicates metrics generated during Picard Mark Duplicates step.
--->
