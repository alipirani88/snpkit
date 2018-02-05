.. _quick-start:

Quickstart
==========

Input
-----

Input is a directory(-readsdir) containing SE/PE reads and a config file where all the configuration settings for the pipeline are set. This config file settings will be used universally on all samples available in readsdir. An example `config https://github.com/alipirani88/variant_calling_pipeline/blob/master/config`_ file with default parameters are included in the pipeline folder. You can customize this config file and provide it with the -config argument. Detailed information in section [Customizing Config file](#customizing-config-file)

.. note::

Apart from standard Miseq/Hiseq fastq naming extensions (R1_001_final.fastq.gz), other acceptable fastq extensions are: R1.fastq.gz/_R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, _1.fastq.gz/.1.fastq.gz. 

Run pipeline on Compute cluster
-------------------------------

The variant calling can be run in parallel on each sample using the -cluster argument. Set the pbs resources such as resources(-l), email(-M), queue(-q), flux_account(-A) and notification(-m) under [scheduler] section in config file. For more details, check out UMICH `flux http://arc-ts.umich.edu/systems-and-services/flux/`_ website.


Possible options for -cluster option(Supported system: pbs):

**local** : This is the default option. When this option is set, the pipeline will analyze each sample one after the another. This will not make use of multiple cores present in your system or clusters. Use this option for few number of samples or for testing purposes.

**parallel-local** :  This option will run the pipeline and analyze samples in parallel but on a local system. This is preferred for less than 20 sample size input and multiple core local system.

**cluster** : This option will run the pipeline on a single cluster. This option is similar to local but will rather run on a cluster node. It will not make use of multiple cores present on a cluster. Use this option for few number of samples or for testing purposes or for few number of large size samples which requires multiple cores to analyze individual samples.

**parallel-cluster** : The pipeline is optimized for this option. When this option is set, the pipeline will run everything on a cluster making use of all the cores available. This option will also submit individual jobs for each sample seperately and whenever needed. When you set this option for each step in the pipeline, make sure all the jobs submitted by the pipeline at each step is completed before proceeding to another step.

Run pipeline
------------

Assuming you want to generate core snps for more than a few hundred samples and run the analysis in parallel on cluster(Time and memory efficient). The default pbs resources used for parallel jobs are: 

::

	nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00


See option resources in scheduler section of `config https://github.com/alipirani88/variant_calling_pipeline/blob/master/config`_ file. Detailed information in section [Customizing Config file](#customizing-config-file)

1. Run variant calling step (All) on a set of PE reads with default parameters

::

	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps All -cluster parallel-cluster


The above command will run variant calling (step 1) pipeline on a set of PE reads residing in test_readsdir. The results will be saved in output directory test_output_core. The config file contains options for some frequently used reference genome. To know which reference genomes are included in config file, look up the [config]() file or check the help menu of the pipeline.

The results of variant calling will be placed in an individual folder generated for each sample in output directory. A log file for each sample will be generated and can be found in each sample folder inside the out directory. A single log file of this step will be generated in main output directory. For more information on log file prefix and convention, please refer [log](#log) section below.

2. Run core_prep step to generate files for core SNP calling.

Run this steps to generate various intermediate files that will be used for generating core SNPs.

::

	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core_prep -cluster parallel-cluster



3. Run core step to generate final core SNP consensus fasta files.

Since this step compares multiple files simultaneously and involves multiple I/O operations, It is recommended to provide higher memory compute resources. 

example:

::

	nodes=1:ppn=4,mem=47000mb,walltime=24:00:00


Replace the resources option in scheduler section of config file with the above line before running the command.

::

	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core -cluster parallel-cluster


4. Generate various pipeline reports and results folder

::
	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps report -cluster parallel-cluster

