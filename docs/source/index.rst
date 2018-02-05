Variant Calling and core SNP diagnostics Pipeline
===============================

Synopsis
--------

**This pipeline calls variants on PE/SE reads provided in a directory and generates core SNP/consensus fasta files that can be used to to build a phylogeny or as an input for Gubbins/Beast analysis**

Installation
------------

Ignore this if you are in snitkin lab. The dependencies are already installed in lab bin_group folder: 

	/nfs/esnitkin/bin_group/variant_calling_bin/. 

Use the python version installed in:
	
	/nfs/esnitkin/bin_group/anaconda2/bin/python

Steps
-----

![alt tag](https://github.com/alipirani88/variant_calling_pipeline/blob/master/pipeline.png)

There are three main steps to generate core SNPs which can be provided with -steps argument and should be run in sequential order.

**1. Variant Calling:** This step will run all the standard variant calling steps on sample read files residing in input reads directory. 
The possible options are:

Option **All** :  This will run all the variant calling steps starting from cleaning the reads to variant calling and generating intermediate files that will be used in later steps. 

Option **clean,align,post-align,varcall,filter,stats** :  This option will run variant calling steps starting from cleaning to variant calling, variant filtering and generating stats.

You can modify this argument to run a part of the pipeline. For example, if you want to restart the pipeline from alignment step, you can do it by supplying the following option: "align,post-align,varcall,filter,stats" which will skip the Trimmomatic cleaning part.


Note: The order of variant calling steps needs to be sequential to avoid any errors. While skipping any of the steps, make sure the results for skipped steps are already present in your output folder.


**2. Preparing files for Core SNP extraction and diagnostics purposes:**


Option **core_prep** : Run this step before running the last core steps. This will prepare all the intermediate data required for the last core step i.e generating the consensus fasta file out of core snps. This can be a time consuming step depending on how closely related the samples are to the reference genome.  


**3. Generate Core SNP consensus and data matrix for diagnostics plots:**

Option **core** : This step will generate core SNP consensus fasta file and a consensus fasta of only core variant positions. Various data matrices will generated at this step that can be used later for diagnostics purposes. 

**4. Generate report for the pipeline:**

Option **report** : This step will generate final core results directory and various reports that will summarize the alignment and core SNP results. 


Quickstart
==========

Input
-----

Input is a directory(-readsdir) containing SE/PE reads and a config file where all the configuration settings for the pipeline are set. This config file settings will be used universally on all samples available in readsdir. An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters are included in the pipeline folder. You can customize this config file and provide it with the -config argument. Detailed information in section [Customizing Config file](#customizing-config-file)

.. note::

Apart from standard Miseq/Hiseq fastq naming extensions (R1_001_final.fastq.gz), other acceptable fastq extensions are: R1.fastq.gz/_R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, _1.fastq.gz/.1.fastq.gz. 

Run pipeline on Compute cluster
-------------------------------

The variant calling can be run in parallel on each sample using the -cluster argument. Set the pbs resources such as resources(-l), email(-M), queue(-q), flux_account(-A) and notification(-m) under [scheduler] section in config file. For more details, check out UMICH [flux](http://arc-ts.umich.edu/systems-and-services/flux/) website.


Possible options for -cluster option(Supported system: pbs):

**local** : This is the default option. When this option is set, the pipeline will analyze each sample one after the another. This will not make use of multiple cores present in your system or clusters. Use this option for few number of samples or for testing purposes.

**parallel-local** :  This option will run the pipeline and analyze samples in parallel but on a local system. This is preferred for less than 20 sample size input and multiple core local system.

**cluster** : This option will run the pipeline on a single cluster. This option is similar to local but will rather run on a cluster node. It will not make use of multiple cores present on a cluster. Use this option for few number of samples or for testing purposes or for few number of large size samples which requires multiple cores to analyze individual samples.

**parallel-cluster** : The pipeline is optimized for this option. When this option is set, the pipeline will run everything on a cluster making use of all the cores available. This option will also submit individual jobs for each sample seperately and whenever needed. When you set this option for each step in the pipeline, make sure all the jobs submitted by the pipeline at each step is completed before proceeding to another step.

Command line options
--------------------

::

	usage: python variant_call.py [-h] -type TYPE -readsdir DIR -outdir OUTPUT_FOLDER -index INDEX [-steps STEPS] -analysis ANALYSIS_NAME [-config CONFIG] [-suffix SUFFIX] [-filenames FILENAMES] [-cluster CLUSTER]

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


Assuming you want to generate core snps for more than a few hundred samples and run the analysis in parallel on cluster(Time and memory efficient). The default pbs resources used for parallel jobs are: 

::
	nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00


See option resources in scheduler section of [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file. Detailed information in section [Customizing Config file](#customizing-config-file)

- Run variant calling step (All) on a set of PE reads with default parameters

::
	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps All -cluster parallel-cluster


The above command will run variant calling (step 1) pipeline on a set of PE reads residing in test_readsdir. The results will be saved in output directory test_output_core. The config file contains options for some frequently used reference genome. To know which reference genomes are included in config file, look up the [config]() file or check the help menu of the pipeline.

The results of variant calling will be placed in an individual folder generated for each sample in output directory. A log file for each sample will be generated and can be found in each sample folder inside the out directory. A single log file of this step will be generated in main output directory. For more information on log file prefix and convention, please refer [log](#log) section below.

- Run core_prep step to generate files for core SNP calling.

Run this steps to generate various intermediate files that will be used for generating core SNPs.

::
	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core_prep -cluster parallel-cluster



- Run core step to generate final core SNP consensus fasta files.

Since this step compares multiple files simultaneously and involves multiple I/O operations, It is recommended to provide higher memory compute resources. 

example:

::
	nodes=1:ppn=4,mem=47000mb,walltime=24:00:00


Replace the resources option in scheduler section of config file with the above line before running the command.

::
	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core -cluster parallel-cluster


