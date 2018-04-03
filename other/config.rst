.. _config:

Customizing Config file
=======================

By default, the pipeline uses config file that comes with the pipeline. Make sure to edit this config file or copy it to your local system, edit it and provide path of this edited config file with -config argument.

::

	cp variant_calling_pipeline/config /Path-to-local/config_edit



The pipeline implements customisable variant calling configurations using config file. Config file can be customised to use your choice of aligner and variant caller by changing two parameters under the section [pipeline]
Currently, The pipeline supports BWA aligner(mem algorithm) for aligning reads to the reference genome and samtools for variant calling.


::
	# Set which tools to use in pipeline:
	[pipeline]
	# Options for Aligner:bwa / smalt / bowtie
	aligner: bwa
	# Options for variant_caller:  gatkhaplotypecaller /samtools
	variant_caller: samtools


Make sure you have downloaded all the dependencies for the pipeline in a folder path provided in binbase option of [bin_path] section.

::

	# Set bin folder path. Please make sure all the executables are placed in bin folder. Also make sure the path for individual tools are correct.
	[bin_path]
	binbase: /nfs/esnitkin/bin_group/variant_calling_bin/


NOTE: Add the required perl libraries(such as in the case of vcftools) PERL5LIB environment variable. For flux users, you can do that by loading perl-modules

::

	module load perl-modules


If you wish to run the jobs on clsuter, make sure you cahnge the necessary schedular parameters in scheduler section shown below: for more information, visit [flux](http://arc-ts.umich.edu/systems-and-services/flux/) homepage.

::

	[scheduler]
	resources: nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00
	email: username@umich.edu
	queue: XXX
	flux_account: XXX
	notification: a



Every tool has its own *_bin option where you can set the folder name in which the tool resides. For example, in the below Trimmomatic section example, the Trimmomatic tool resides in /Trimmomatic/ folder that is set with trimmomatic_bin option which in itself resides in /nfs/esnitkin/bin_group/variant_calling_bin/ folder that was set in binbase option above.

::

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


Parameters for each tools can be customised under the 'tool_parameter' attribute of each tool in config file.


For example, to change the minadapterlength parameter of Trimmomatic from 8 to 10, replace minadapterlength of 8 with suppose 10 and restart the pipeline.
