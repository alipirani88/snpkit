# Variant Calling and core SNP diagnostics Pipeline

## Synopsis

This is a highly customisable, automated variant detection pipeline that can be easily deployed for infectious disease outbreak investigations and other clinical microbiology projects.  

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Steps](#steps)
- [Input](#input)
- [Command line options](#command-line-options)
- [Run pipeline on compute cluster](#run-pipeline-on-compute-cluster)
- [Output files](#output-files)
- [Customizing the config file](#customizing-the-config-file)
- [Log](#log)
- [Bonus Ducks](https://knowyourmeme.com/memes/bonus-ducks)

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

**2. index:** a reference genome index name as specified in a config file. For example; if you have set the reference genome path in config file as shown below, then the required value for command line argument -index would be -index KPNIH1

```
[KPNIH1]
# path to the reference genome fasta file.
Ref_Path: /nfs/esnitkin/bin_group/variant_calling_bin/reference/KPNIH1/
# Name of reference genome fasta file.
Ref_Name: KPNIH1.fasta
```

Here, Ref_Name is the reference genome fasta file located in Ref_Path. Similarly, if you want to use a different version of KPNIH reference genome, you can create a new section with a different index name.

```
[KPNIH1_new]
# path to the reference genome fasta file.
Ref_Path: /nfs/esnitkin/bin_group/variant_calling_bin/reference/KPNIH1_new/
# Name of reference genome fasta file.
Ref_Name: KPNIH1_new.fasta
```

For more information, refer to [Customizing the config file](#customizing-the-config-file).

**3. config:** a config file to set pipeline configuration settings such as setting up environment path for various tools, path to reference genomes and filter parameters. 

The config file is a YAML format file that stores data in KEY: VALUE pair. This settings will be applied globally on all variant call jobs. An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters is included in code folder. You can customize this config file and provide it with the -config argument. An example parameter setting is shown below where we are setting the bin directory path. This is another way of telling the pipeline that all the tools required for variant calling are located in "binbase" directory of bin_path section.

```
[bin_path]
binbase: /nfs/esnitkin/bin_group/variant_calling_bin/
```



## Steps

The piepline calls variants on Illumina paired end(PE) / single end (SE) reads provided in a directory and generate a phylogenetic tree from recombinant filtered high quality variants found against the reference genome.


The pipeline is divided into two individual steps (-steps option) which should be run in sequential order (see below for command line arguments). 

**1. Variant Calling:** This step will run all standard variant calling steps on sample reads residing in the input reads directory.

The possible options are:

Option ***All***:  This will run all variant calling steps from read trimming to variant calling. 

Option ***clean,align,post-align,varcall,filter,stats***:  Use these options, if you want to run individual steps involved in variant calling steps. This will run variant calling steps starting from read cleaning/trimming, alignment to reference genome, post-alignment sam/bam format conversion, variant calling/filtering and finally generates mapping/variant statistics.

You can run a part of the pipeline by customizing the order of the -steps argument. For example, to skip the Trimmomatic cleaning part and instead run the pipeline from alignment step, run it with the following options: 

"align,post-align,varcall,filter,stats"

Note: The order of variant calling steps needs to be sequential. If skipping any of the steps, make sure the previous steps finished without any errors and their results are already present in output folder.

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/pipeline_All.png)

**2. Generate various Core/Non-core SNP consensus, SNP/Indel Matrices and recombination filtering with Gubbins:**

Option ***core_All***: This will run all the core SNP consensus/matrix generating steps. Once the core variants and matrices are generated, the command will run gubbins and RAxML/iqtree jobs. 

Note: Samples not meeting minimum Depth of Coverage threshold (10X) will be removed before running this step. 

Some of the types of results generated by this step are:

- prefix_core_results directory under the output directory which will be the final results folder.
- intermediate data files required for generating core SNP matrix/consensus in vcf/fasta format.
- Various data matrices will be generated during this step that can be used for diagnosing variant filter criterias and their impact on overall distribution of core variants.
- Gubbins recombination filtered consensus fasta files and RaxML/Iqtree trees generated from this recombination filtered consensus.

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/pipeline_core_All.png)

## Command line options

```

usage: variant_call.py [-h] -type TYPE -readsdir DIR -outdir OUTPUT_FOLDER
                       -index INDEX [-steps STEPS] -analysis ANALYSIS_NAME
                       [-gubbins_env GUBBINS_ENV] [-config CONFIG]
                       [-suffix SUFFIX] [-filenames FILENAMES]
                       [-cluster CLUSTER] [-clean]
                       [-extract_unmapped EXTRACT_UNMAPPED] [-datadir DATADIR]
                       [-snpeff_db SNPEFF_DB] [-debug_mode DEBUG_MODE]
                       [-gubbins GUBBINS] [-outgroup OUTGROUP]
                       [-downsample DOWNSAMPLE]
                       [-coverage_depth COVERAGE_DEPTH] [-scheduler SCHEDULER]

Variant Calling pipeline for Illumina PE/SE data.

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -type TYPE            Type of reads: SE or PE
  -readsdir DIR         Path to Sequencing Reads Data directory. NOTE: Provide full/absolute path.
  -outdir OUTPUT_FOLDER
                        Output Folder Path ending with output directory name to save the results. Creates a new output directory path if it doesn't exist. NOTE: Provide full/absolute path.
  -index INDEX          Reference Index Name. Most Frequently used reference genomes index options: KPNIH1 | MRSA_USA_300 | MRSA_USA_100 | CDIFF_630 | paris Make sure the paths are properly set in config file
  -steps STEPS          Variant Calling Steps in sequential order.
                        1.   All: Run all variant calling  steps starting from trimming the reads, mapping, post-processing the alignments and calling variants;
                        2.   core_All: Extract core snps and generate different types of alignments, SNP/Indel Matrices and diagnostics plots.
  -analysis ANALYSIS_NAME
                        Unique analysis name that will be used as prefix to saving results and log files.

Optional arguments:
  -gubbins_env GUBBINS_ENV
                        Name of the Gubbins Raxml Iqtree environment to load for Phylogenetic analysis
  -config CONFIG        Path to Config file, Make sure to check config settings before running pipeline
  -suffix SUFFIX        Fastq reads suffix such as fastq, fastq.gz, fq.gz, fq; Default: fastq.gz
  -filenames FILENAMES  fastq filenames with one single-end filename per line. 
                        If the type is set to PE, it will detect the second paired-end filename with the suffix from first filename. 
                        Useful for running variant calling pipeline on selected files in a reads directory or extracting core snps for selected samples in input reads directory. 
                        Otherwise the pipeline will consider all the samples available in reads directory.
  -cluster CLUSTER      Run variant calling pipeline in local or cluster mode.
                        Default: local.
                        Set your specific hpc cluster parameters in config file under the [scheduler] section. Supports PBS/SLURM scheduling system.
  -clean                clean up intermediate files. Default: ON
  -extract_unmapped EXTRACT_UNMAPPED
                        Extract unmapped reads, assemble it and detect AMR genes using ariba
  -datadir DATADIR      Path to snpEff data directory
  -snpeff_db SNPEFF_DB  Name of pre-build snpEff database to use for Annotation
  -debug_mode DEBUG_MODE
                        yes/no for debug mode
  -gubbins GUBBINS      yes/no for running gubbins
  -outgroup OUTGROUP    outgroup sample name
  -downsample DOWNSAMPLE
                        yes/no: Downsample Reads data to default depth of 100X or user specified depth
  -coverage_depth COVERAGE_DEPTH
                        Downsample Reads to this user specified depth
  -scheduler SCHEDULER  Type of Scheduler for generating cluster jobs: PBS, SLURM, LOCAL

```



## Run pipeline on compute cluster


Variant calling can be run in parallel on each sample using the -cluster argument. Set pbs resources such as resources (-l), email (-M), queue (-q), flux_account (-A) and notification (-m) under the [scheduler] section in config file. 

For more details, refer to the UMICH [flux](http://arc-ts.umich.edu/systems-and-services/flux/) website for a detailed explanation of each pbs specification.


Possible options for the -cluster option (Supported system: pbs):

***cluster***: The pipeline is optimized for this option. When this option is set, the pipeline will run variant call jobs for each sample on individual small compute clusters or on a single large cluster for core snp generation step. 

***local***: This is the default option. When this option is set, the pipeline will analyze each sample one after the another. This will not make use of multiple cores present in your system or clusters. Use this option for small numbers of samples or for testing purposes.

***parallel-local***:  This option will run the pipeline and analyze samples in parallel but on a local system. This is preferred for an input sample size of less than 20  and a multiple core local system.





<!--- Input is a directory (-readsdir) containing SE/PE reads and a config file where all the configuration settings for the pipeline are set. These config file settings will be used universally on all samples available in readsdir. An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters is included in the pipeline folder. You can customize this config file and provide it with the -config argument. -->

## Output Files

All the final results will be saved under date_time_core_results directory under the output folder. 
```

2020_03_01_16_12_11_core_results
├── core_snp_consensus
├── data_matrix
├── gubbins
├── qc_report
└── README


```


**gubbins** contains different combinations of core/non-core multi-fasta alignments. Alignments with "gubbins.fa" in their extension will be used as an input for Gubbins. Once the gubbins jobs finished, the pipeline runs RaxML and Iqtree on gubbins generated recombination filtered consensus fasta files. Raxml and iqtree results generated using these files will be placed in raxml_results and iqtree_results folder respectively.

The MSA file that is used for gubbins/and/iqtree is genome_aln_w_alt_allele_unmapped_gubbins.fa


- genome_aln_w_alt_allele_unmapped_gubbins.fa: This genome alignment cotains all the variants that were called against each sample. It contains all core variant positions, non-core variant positions will be substituted with a dash (denoting unmapped) and N’s for not meeting hard variant filter thresholds or falling under functional class filters. Non-core variants that met all filter thresholds will be substituted with variant allele. 

The different types of variants and the symbols that are used for each variant positions are described below:

1. If a position is unmapped in a sample, then it will be denoted by '-'
2. If a position did not meet a hard filter criteria in a sample then it will be denoted by 'N' in that particular sample.
3. If a position falls under any of the functional class filter such as phage region, repeat region or custom mask region, then it will be denoted by 'N'


**data_matrix** This folder contains different types of data matrices and reports that can be queried for variant diagnostics/QC plots. 

- matrices/SNP_matrix_allele.csv and Indel_matrix_allele.csv: contain allele information for each unique variant position (rownames) called in each individual sample (columns). Rownames are the position where a variant was called in one of the given samples and its associated annotations


- matrices/SNP_matrix_allele_new.csv will also contain allele information for each unique variant position (rownames) called in each individual sample (columns) except the positions that were unmapped and filtered out will be substituted with dash (-) and N's.

- matrices/SNP_matrix_code.csv and Indel_matrix_code.csv: contain one of the seven status codes for each unique variant position (rownames) called in individual sample (columns). 

- Functional_annotation_results/phage_region_positions.txt contains positions identified by Phaster as phage regions.

- Functional_annotation_results/repeat_region_positions.txt contains positions identified by nucmer as tandem repeats.

- Functional_annotation_results/mask_positions.txt contains positions set by the user to mask from the final core position list.

- Functional_annotation_results/Functional_class_filter_positions.txt is an aggregated unique list of positions that fall under repeat, mask and phage region (REPEAT_MASK_PHAGE).

- Row annotation: The rows in the matrix represent a variant position and its associated annotation divided into two parts seperated by a semi colon. An example row annotation is shown below:

```

Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand;ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product	Sample_name

Coding SNP at 31356 > T functional=NULL_NULL_NULL locus_tag=USA300HOU_0023 strand=+;T|stop_gained|HIGH|USA300HOU_0023|c.373C>T|p.Gln125*|373/2361|125/786|null or hypothetical protein|possible 5'-nucleotidase;	0
```

The first part contains information such as type of SNP, position, variant allele found, functional annotation, locus tag for the gene and strand. Second part contains snpEff annotation for the variant found and its impact on the gene. 


The different status codes are: 

| Code | Description |
| --------- | ----------- |
| -1 | unmapped |
| 0  | reference allele |
| 1  | core |
| 2  | filtered |
| 3  | non-core or true variant but filtered out due to another sample |
| -2 | Phage Region |
| -3 | FQ Region Masked |
| -4 | MQ Region Masked |


Some toy examples of how codes are arranged for different type of variants and how they would be represented in allele matrix is shown below:

- Reference Allele:

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/reference_allele.png)

- Unmapped Positions:

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/unmapped_positions.png)

- Core Variants:

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/core_variants.png)

- Filtered Positions:

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/filtered_positions.png)



**core_snp_consensus** The core_vcf folder under this directory contains annotated core vcf files that were used for generating core SNP consensus fasta results. Other folders contain different combination of core/non-core consensus fasta files for individual samples. The consensus file from these folders are concatenated to generate the multiple sequence alignment file which are then placed in gubbins folder and are used as an input for gubbins jobs.

<!--
**Bargraph matrices**

bargraph_counts.txt and bargraph_indel_counts.txt contain distributions of all the variant positions called in each individual sample. Each color represents a type of variant filter. The bar corresponds to the number of variants that got filtered out due to that particular hard filter parameter in each sample.

bargraph_indel_percentage.txt and bargraph_percentage.txt contain same distribution as the counts bargraphs but the counts are transformed into percentages.

These matrices can be used for QC checks and to understand the effects of different variant filters on the total number of core variants. 

Note: Samples with an unusual bar height should be considered outlier samples.

Run the generate_diagnostics_plots.R script created inside the data_matrix folder to generate various QC plots. 

Requires: ggplot2 and heatmap.3

```

module load R

Rscript generate_diagnostics_plots.R 

```

| File | Description |
| --------- | ----------- |
| barplot.pdf |  Distribution of filter-pass variant positions (variants observed in all the samples) in each sample. Colors represent the filter criteria that caused them to get filtered out in that particular sample.|
| barplot_DP.pdf | Distribution of filter-pass variant positions in each sample. Colors represent the read-depth range that they fall in. |
| temp_Only_filtered_positions_for_closely_matrix_FQ.pdf | Heatmap spanning the reference genome that shows positions that were filtered out due to low FQ values |
| DP_position_analysis.pdf | same information as in barplot_DP.pdf but shown in heatmap format|
| temp_Only_filtered_positions_for_closely_matrix_DP.pdf | Heatmap spanning the reference genome that shows positions that were filtered out due to low DP values |


- barplot

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/barplot.png)

- barplot_DP

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/img/barplot_DP.png)


<!-- 
barplot.pdf

![click here](https://github.com/alipirani88/variant_calling_pipeline/blob/master/modules/variant_diagnostics/R_scripts/barplot.pdf)


barplot_DP.pdf 

![click here](https://github.com/alipirani88/variant_calling_pipeline/blob/master/modules/variant_diagnostics/R_scripts/barplot_DP.pdf)

 

***Core variant positions***
Only_ref_variant_positions_for_closely_matrix.txt
Only_ref_variant_positions_for_closely_without_functional_filtered_positions





| Extension | Description |
| --------- | ----------- |
| barplot.pdf |  Distribution of filter-pass variant positions (variants observed in all the samples) in each sample. colors represents the filter criteria that caused them to get filtered out in that particular sample.|
| barplot_DP.pdf | Distribution of filter-pass variant positions in each sample. color represents the read-depth range that they fall in. |
| temp_Only_filtered_positions_for_closely_matrix_FQ.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low FQ values |
| DP_position_analysis.pdf | same information as in barplot_DP.pdf but shown in heatmap format|
| temp_Only_filtered_positions_for_closely_matrix_DP.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low DP values |


barplot.pdf

![click here](https://github.com/alipirani88/variant_calling_pipeline/blob/master/modules/variant_diagnostics/R_scripts/barplot.pdf)

barplot_DP.pdf 

![click here](https://github.com/alipirani88/variant_calling_pipeline/blob/master/modules/variant_diagnostics/R_scripts/barplot_DP.pdf)


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

-->


-->
## Tips and Tricks

- Run only part of a variant calling step

OPTIONAL: If you want to rerun the above analysis with different variant call/filter parameters (i.e skip the read cleaning, alignment and post-alignment steps), you can do it as follows:

NOTE: OPTIONAL
```

python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps varcall,filter,stats -cluster parallel-cluster

```

- Run the variant calling steps on selected samples provided in filenames.txt. filenames.txt should contain fastq filenames with one single-end filename per line.

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps varcall,filter,stats -cluster parallel-cluster -filenames filenames.txt
```

- Run core_prep and core steps on selected files

If you decide to exclude some samples from the core snp analysis step (step 3), there is a way to run core_prep and core steps on selected files without running the entire pipeline. 

```

python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core -cluster parallel-cluster -filenames filenames_custom

```

## Customizing the config file

By default, the pipeline uses the config file that comes with the pipeline. Make sure to edit this config file or copy it to your local system, edit it and provide the path of this edited config file with the -config argument.

```

cp variant_calling_pipeline/config /Path-to-local/config_edit

```

The pipeline implements customisable variant calling configurations using the config file. The config file can be customised to use your choice of aligner and variant caller by changing two parameters under the section [pipeline]
Currently, The pipeline supports BWA aligner (mem algorithm) for aligning reads to the reference genome and samtools for variant calling.

```
# Set which tools to use in pipeline:
[pipeline]
# Options for Aligner:bwa / smalt / bowtie
aligner: bwa
# Options for variant_caller:  gatkhaplotypecaller /samtools
variant_caller: samtools

```

Make sure you have downloaded all the dependencies for the pipeline in a folder path provided in the binbase option of the [bin_path] section.

```
# Set bin folder path. Please make sure all the executables are placed in the bin folder. Also make sure the path for individual tools are correct.
[bin_path]
binbase: /nfs/esnitkin/bin_group/variant_calling_bin/
```

NOTE: Add the required perl libraries (such as in the case of vcftools) PERL5LIB environment variable. For flux users, you can do that by loading perl-modules

```
module load perl-modules
```

If you wish to run the jobs on a cluster, make sure you cahnge the necessary scheduler parameters in the scheduler section shown below: for more information, visit the [flux](http://arc-ts.umich.edu/systems-and-services/flux/) homepage.

```

[scheduler]
resources: nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00
email: username@umich.edu
queue: XXX
flux_account: XXX
notification: a

```

Every tool has its own *_bin option where you can set the folder name in which the tool resides. For example, in the below Trimmomatic section example, the Trimmomatic tool resides in the /Trimmomatic/ folder that is set with the trimmomatic_bin option which in itself resides in /nfs/esnitkin/bin_group/variant_calling_bin/ folder that was set in the binbase option above.

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

Parameters for each of the tools can be customised under the 'tool_parameter' attribute of each tool in the config file.


For example, to change the minadapterlength parameter of Trimmomatic from 8 to 10, replace minadapterlength of 8 with suppose 10 and restart the pipeline.

## Log:

The pipeline generates a log file following the naming convention: yyyy_mm_dd_hrs_mins_secs_analysisname.log.txt and tracks each event/command. The log file sections follow standard [Python logging conventions](https://docs.python.org/2/howto/logging.html): 

***INFO*** to print STDOUT messages; 

***DEBUG*** to print commands run by the pipeline, 

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


