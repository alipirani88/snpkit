# Variant Calling and core SNP diagnostics Pipeline

## Synopsis

This is a highly customisable, automated variant detection pipeline that can be easily deployed for infectious disease outbreak investigations and other clinical microbiology projects.  

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Input](#input)
- [Command line options](#command-line-options)
- [Output files](#output-files)


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


