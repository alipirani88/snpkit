# Variant Calling and core SNP diagnostics Pipeline

## Synopsis


The pipeline call variants on Illumina PE/SE reads provided in a directory and generates various combination of core/non-core consesensus fasta files that can be used for phylogenetic reconstruction or as an input for Gubbins/Beast analysis.

## Contents

- [Installation](#installation)
- [Input](#input-requirements)
- [Steps](#steps)
- [Command line options](#command-line-options)
- [Run pipeline on Compute cluster](#run-pipeline-on-compute-cluster)
- [Quick Start](#quick-start)
- [Output Files](#output-files)
- [Customizing Config file](#customizing-config-file)
- [Log](#log)
- [Bonus Ducks](#bonus-ducks)

## Installation

The dependencies are already installed in Snitkin lab bin_group folder:

```
/nfs/esnitkin/bin_group/variant_calling_bin/
```

Requires Python2 version:

```
/nfs/esnitkin/bin_group/anaconda2/bin/python
```

## Input Requirements

- readsdir: folder containing SE/PE reads. Apart from standard Miseq/Hiseq fastq naming extensions (R1_001_final.fastq.gz), other acceptable fastq extensions are: R1.fastq.gz/_R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, _1.fastq.gz/.1.fastq.gz.

- config: a config file to set pipeline configuration settings such as setting up environment path for various tools, path to reference genomes and filter parameters. This settings will be applied globally on all variant call jobs. An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters is included in code folder. You can customize this config file and provide it with the -config argument.

For more information, refer [Customizing Config file](#customizing-config-file)

<!--- Input is a directory(-readsdir) containing SE/PE reads and a config file where all the configuration settings for the pipeline are set. This config file settings will be used universally on all samples available in readsdir. An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters are included in the pipeline folder. You can customize this config file and provide it with the -config argument. -->



## Steps


![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/pipeline.png)

Pipeline is divided into three individual steps(-steps option) which should be run in sequential order. 

**1. Variant Calling:** This step will run all standard variant calling steps on sample reads residing in input reads directory.

The possible options are:

Option ***All*** :  This will run all variant calling steps starting from read trimming to variant calling. 

Option ***clean,align,post-align,varcall,filter,stats*** :  This option will run variant calling steps starting from read trimming, alignment to reference genome, post-alignment format conversion, variant call/filter and finally generates mapping/variant statistics.

You can run a part of the pipeline by customizing the order of -steps argument. For example, to skip Trimmomatic cleaning part, run it with following options: "align,post-align,varcall,filter,stats" 


Note: The order of variant calling steps needs to be sequential. While skipping any of the steps, make sure those skipped steps finished without any errors.


**2. Preparing files for Core SNP extraction and diagnostics purposes:**


Option ***core_prep*** : Run this step before running the last core steps. This will prepare all the intermediate data required for generating core SNP matrix/consensus.


**3. Generate Core SNP consensus and data matrix for diagnostics plots:**

Option ***core*** : This step will generate core SNP/Indel Matrix and different types of consensus fasta files. Various data matrices will be generated during this step that can be used for diagnosing variant filter criterias. 

**4. Generate report and aggregate results for the pipeline:**

Option ***report*** : This step will aggregate the results in prefix_core_results directory under the output directory.

**5. Phylogenetic reconstruction and recombination filtering using FastTree/RAxML/Gubbins:**

Option ***tree*** : This step will generate FastTree/RAxML tree from pre-recombination filtered consensus files. If -gubbins option is set to yes, it will run gubbins on all final consensus files..

## Command line options

```

usage: variant_call.py [-h] -type TYPE -readsdir DIR -outdir OUTPUT_FOLDER
                       -index INDEX [-steps STEPS] -analysis ANALYSIS_NAME
                       [-config CONFIG] [-suffix SUFFIX]
                       [-filenames FILENAMES] [-cluster CLUSTER]
                       [-clean CLEAN] [-extract_unmapped EXTRACT_UNMAPPED]
                       [-datadir DATADIR] [-snpeff_db SNPEFF_DB]
                       [-debug_mode DEBUG_MODE] [-gubbins GUBBINS]

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
                        1.   All: This will run all the steps starting from cleaning the reads to variant calling;
                        2.   clean,align,post-align,varcall,filter,stats : This will also run all steps starting from cleaning to variant calling. 
                        You can also run part of the pipeline by giving "align,post-align,varcall,filter,stats" which will skip the cleaning part.
                        The order is required to be sequential while using this option. Also, while skipping any of the step make sure you have results already present in your output folder.
                        3.   coverage_depth_stats: Run Only Depth of Coverage Stats module after cleaning and read mapping steps
                        4.   core_prep: Run this step before running the core steps. This will prepare the data required for generating core SNPs
                        5.   core: extract core snps and generate diagnostics plot data matrices to explore filtered snps.
  -analysis ANALYSIS_NAME
                        Unique analysis name that will be used as prefix to saving results and log files.

Optional arguments:
  -config CONFIG        Path to Config file, Make sure to check config settings before running pipeline
  -suffix SUFFIX        Fastq reads suffix such as fastq, fastq.gz, fq.gz, fq; Default: fastq.gz
  -filenames FILENAMES  fastq filenames with one single-end filename per line. 
                        If the type is set to PE, it will detect the second paired-end filename with the suffix from first filename. 
                        Useful for running variant calling pipeline on selected files in a reads directory or extracting core snps for selected samples in input reads directory. 
                        Otherwise the pipeline will consider all the samples available in reads directory.
  -cluster CLUSTER      Run variant calling pipeline in one of the four modes. Default: local. Suggested mode for core snp is cluster that will run all the steps in parallel with the available cores. Make sure to provide a large memory node for this option
                        The possible modes are: cluster/parallel-cluster/parallel-local/local
                        cluster: Runs all the jobs on a single large cluster. This will mimic the local run but rather on a large compute node.
                        parallel-cluster: Submit variant call jobs for each sample in parallel on compute nodes. This mode is no available for core snp extraction step.
                        parallel-local: Run variant call jobs for each sample in parallel locally.
                        local: Run variant call jobs locally.
                        Make Sure to check if the [scheduler] section in config file is set up correctly for your cluster.
  -clean CLEAN          clean up intermediate files. Default: OFF
  -extract_unmapped EXTRACT_UNMAPPED
                        Extract unmapped reads, assemble it and detect AMR genes using ariba
  -datadir DATADIR      Path to snpEff data directory
  -snpeff_db SNPEFF_DB  Name of pre-build snpEff database to use for Annotation
  -debug_mode DEBUG_MODE
                        yes/no for debug mode
  -gubbins GUBBINS      yes/no for running gubbins

```

## Run pipeline on Compute cluster


Variant calling can be run in parallel on each sample using -cluster argument. Set pbs resources such as resources(-l), email(-M), queue(-q), flux_account(-A) and notification(-m) under [scheduler] section in config file. 

For more details, refer UMICH [flux](http://arc-ts.umich.edu/systems-and-services/flux/) website for detailed explanation of each pbs specifications.


Possible options for -cluster option(Supported system: pbs):

***local*** : This is the default option. When this option is set, the pipeline will analyze each sample one after the another. This will not make use of multiple cores present in your system or clusters. Use this option for few number of samples or for testing purposes.

***parallel-local*** :  This option will run the pipeline and analyze samples in parallel but on a local system. This is preferred for less than 20 sample size input and multiple core local system.

***cluster*** : This option will run the pipeline on a single cluster. This option is similar to local but will rather run on a cluster node. It will not make use of multiple cores present on a cluster. Use this option for few number of samples or for testing purposes or for few number of large size samples which requires multiple cores to analyze individual samples.

***parallel-cluster*** : The variant call step(All) is optimized for this option. When this option is set, the pipeline will run variant call jobs for each sample on individual compute cluster. When you set this option for each step in the pipeline, make sure all the jobs submitted by the pipeline at each step is completed before proceeding to another step. You can check the status of the job with: 

```
qstat -u USERNAME  
```

Note: Use parallel-cluster mode for All/clean,align,post-align,varcall,filter,stats steps. Use cluster/parallel-local mode for core_prep and core steps. parallel-cluster and cluster will be merged in next release.

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
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core_prep -cluster cluster

```

- Run core step to generate final core SNP consensus fasta files.

Since this step compares multiple files simultaneously and involves multiple I/O operations, It is recommended to provide higher memory compute resources. 

example:

```
nodes=1:ppn=4,mem=47000mb,walltime=24:00:00
```

Replace the resources option in scheduler section of config file with the above line before running the command.

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core -cluster cluster

```

- Run report step to aggregate results under prefix_core_results folder

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps report -cluster cluster

```

- Run tree step to generate FastTree/RAxML phylogenetic tree and recombination filterating using Gubbins

```
python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps tree -cluster cluster -gubbins yes

```

## Output Files

All the final results from report step will be saved under date_time_core_results directory under the output folder. 


```

2018_01_13_13_18_03_core_results
├── core_snp_consensus
│   ├── consensus_allele_variant_positions
│   ├── consensus_ref_allele_unmapped_variant
│   ├── consensus_ref_allele_variant_positions
│   ├── consensus_ref_variant_positions
│   ├── consensus_variant_positions
│   └── core_vcf
├── data_matrix
│   ├── Indel_matrix_allele.csv
│   ├── Indel_matrix_code.csv
│   ├── SNP_matrix_allele.csv
│   ├── SNP_matrix_allele_new.csv
│   ├── SNP_matrix_code.csv
│   ├── phage_region_positions.txt
│   ├── repeat_region_positions.txt
│   ├── mask_positions.txt
│   ├── Functional_class_filter_positions.txt
│   ├── bargraph_counts.txt
│   ├── bargraph_indel_counts.txt
│   ├── bargraph_indel_percentage.txt
│   ├── bargraph_percentage.txt
│   ├── snpEff_results
│   └── temp
├── gubbins
│   ├── 2018_07_11_14_45_01_KPNIH1_allele_var_consensus.fa
│   ├── 2018_07_11_14_45_01_KPNIH1_ref_allele_unmapped_consensus.fa
│   ├── 2018_07_11_14_45_01_KPNIH1_ref_allele_var_consensus.fa
│   ├── 2018_07_11_14_45_01_KPNIH1_ref_var_consensus.fa
│   └── 2018_07_11_14_45_01_KPNIH1_var_consensus.fa
└── trees

```

### core_snp_consensus
This folder contain different combination of core/non-core consensus fasta for individual samples. core_vcf folder contain annotated vcf files used for generating the consensus. 

### gubbins: 
This folder contain different combination of core/non-core multi-fasta consensus generated by merging individual consensus files in core_snp_consensus folder. 

**var_consensus.fa:** consensus generated from core variant positions.

**ref_var_consensus.fa** full consensus relative to the reference genome generated from core variant positions (final core variants + reference alleles + non-core variants). 

Note: non-core variant in this consensus refers to positions where a variant was observed but was filtered out for not meeting one of three criteria. Therefore, it will be replaced by a reference allele observed at this position.

- core variant criteria: Since this position was filtered out in another sample, this position will be regarded as non-core position.
- hard variant filters: If the variant called at this position didn't meet hard variant filters such as DP(depth), MQ(Mapping quality), FQ(all reads supporting one variant), QUAL(variant quality), proximate_snp (the variant was in close proximity to another variant by 10 bp range), they would be filtered out from the final core SNPs.
- functional filters: If the position fall in one of the three functional class (repetitive, masked or phage region), they would be filtered out from final core SNPs.

**allele_var_consensus.fa:** consensus generated from all unique variant positions that was called in any sample regardless of it being a core/non-core. This will contain all core variant positions, non-core variant positions will be substituted with dash (denoting unmapped) and N’s for not meeting hard variant filters or falling under functional class filter. non-core variant that meet all filter parameters will be substituted with variant allele.

**ref_allele_var_consensus.fa** full consensus relative to the reference genome generated from all unique variant positions that was called in any sample regardless of it being a core/non-core. This will contain all core variant positions, non-core variant positions will be substituted with dash (denoting unmapped) and N’s for not meeting hard variant filters or falling under functional class filter. non-core variant that meet all filter parameters will be substituted with variant allele.

**ref_allele_unmapped_consensus.fa:** full consensus relative to the reference genome generated from all unique variant positions that was called in any sample. This will contain all core variant positions, non-core variant positions will be substituted with dash (denoting unmapped) and N’s for not meeting hard variant filters or falling under functional class filter. non-core variant that meet all filter parameters will be substituted with variant allele. Positions that were unmapped relative to the reference genome will be denoted by dashes.


### data_matrix
This folder contains different type of data matrices and reports that can be queried for variant diagnostics/QC plots. 

**SNP/Indel Matrix:**

SNP_matrix_allele.csv and Indel_matrix_allele.csv: contain allele information for each unique variant position(rownames) called in individual sample(columns)

SNP_matrix_code.csv and Indel_matrix_code.csv: contain one of the five status code for each unique variant position(rownames) called in individual sample(columns). 

Different status codes are: 

unmapped = -1

reference allele = 0

core = 1

filtered = 2

non-core or True variant but filtered out due to another sample = 3


**Functional class:**

phage_region_positions.txt contain positions identified by Phaster as phage region.

repeat_region_positions.txt contain positions identified by nucmer as tandem repeats.

mask_positions.txt contain positions set by user to mask from final core position list.

Functional_class_filter_positions.txt is an aggregated unique list of positions that fall under repeat, mask and phage region.

**Bargraph matrices**

bargraph_counts.txt and bargraph_indel_counts.txt contain distribution of all the variant positions called in individual sample. Each color represents a type of variant filter. The bar corresponds to the number of variants that got filtered out due to that particular hard filter parameter in each sample.

bargraph_indel_percentage.txt and bargraph_percentage.txt contain same distribution but the counts are transformed into percentage.

These matrices can be used for QC checks and understand effects of different variant filters on the total number of core variants. 

Note: Sample with unusual bar height should be considered as outlier samples.

Run generate_diagnostics_plots.R script created inside the data_matrix folder to generate various QC plots. 

Require: ggplot2 and heatmap.3

```

module load R

Rscript generate_diagnostics_plots.R 

```

| Extension | Description |
| --------- | ----------- |
| . barplot.pdf |  Distribution of filter-pass variant positions(variants observed in all the samples) in each sample. colors represents the filter criteria that caused them to get filtered out in that particular sample.|
| . barplot_DP.pdf | Distribution of filter-pass variant positions in each sample. color represents the read-depth range that they fall in. |
| . temp_Only_filtered_positions_for_closely_matrix_FQ.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low FQ values |
| . DP_position_analysis.pdf | same information as in barplot_DP.pdf but shown in heatmap format|
| . temp_Only_filtered_positions_for_closely_matrix_DP.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low DP values |


barplot.pdf

![click here](https://github.com/alipirani88/variant_calling_pipeline/blob/master/modules/variant_diagnostics/R_scripts/barplot.pdf)

barplot_DP.pdf 

![click here](https://github.com/alipirani88/variant_calling_pipeline/blob/master/modules/variant_diagnostics/R_scripts/barplot_DP.pdf)


<!-- 

***Core variant positions***
Only_ref_variant_positions_for_closely_matrix.txt
Only_ref_variant_positions_for_closely_without_functional_filtered_positions





| Extension | Description |
| --------- | ----------- |
| barplot.pdf |  Distribution of filter-pass variant positions(variants observed in all the samples) in each sample. colors represents the filter criteria that caused them to get filtered out in that particular sample.|
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

## Tips and Tricks

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

