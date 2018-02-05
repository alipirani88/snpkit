.. _command-line-options:

Command line options
====================

::

	usage: python variant_call.py [-h] -type TYPE -readsdir DIR -outdir OUTPUT_FOLDER -index INDEX [-steps STEPS] -analysis ANALYSIS_NAME [-config CONFIG] [-suffix SUFFIX] [-filenames FILENAMES] [-cluster CLUSTER]

Required arguments
------------------

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

Optional arguments
------------------

  -config       Path to Config file, Make sure to check config settings before running pipeline
  -suffix       Fastq reads suffix such as fastq, fastq.gz, fq.gz, fq; Default: fastq.gz
  -filenames    fastq filenames with one single-end filename per line. if the type is set to PE, it will detect the second paired-end filename with the suffix from first filename.
  -cluster      Run variant calling pipeline in one of the four modes. Default: local. The possible modes are: cluster/parallel-cluster/parallel-local/local

Compute cluster options
-----------------------

	**local** : This is the default option. When this option is set, the pipeline will analyze each sample one after the another. This will not make use of multiple cores present in your system or clusters. Use this option for few number of samples or for testing purposes.

	**parallel-local** :  This option will run the pipeline and analyze samples in parallel but on a local system. This is preferred for less than 20 sample size input and multiple core local system.

	**cluster** : This option will run the pipeline on a single cluster. This option is similar to local but will rather run on a cluster node. It will not make use of multiple cores present on a cluster. Use this option for few number of samples or for testing purposes or for few number of large size samples which requires multiple cores to analyze individual samples.

	**parallel-cluster** : The pipeline is optimized for this option. When this option is set, the pipeline will run everything on a cluster making use of all the cores available. This option will also submit individual jobs for each sample seperately and whenever needed. When you set this option for each step in the pipeline, make sure all the jobs submitted by the pipeline at each step is completed before proceeding to another step.
