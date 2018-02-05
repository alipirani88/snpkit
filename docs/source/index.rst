Variant Calling and core SNP diagnostics Pipeline
===============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Synopsis
--------

**This pipeline calls variants on PE/SE reads provided in a directory and generates core SNP/consensus fasta files that can be used to to build a phylogeny or as an input for Gubbins/Beast analysis**

Installation
------------

Ignore this if you are in snitkin lab. The dependencies are already installed in lab bin_group folder: 

/nfs/esnitkin/bin_group/variant_calling_bin/. 

Use the python version installed in:

/nfs/esnitkin/bin_group/anaconda2/bin/python

Input
-----

Input is a directory(-readsdir) containing SE/PE reads and a config file where all the configuration settings for the pipeline are set. This config file settings will be used universally on all samples available in readsdir. An example [config](https://github.com/alipirani88/variant_calling_pipeline/blob/master/config) file with default parameters are included in the pipeline folder. You can customize this config file and provide it with the -config argument. 

Detailed information in section [Customizing Config file](#customizing-config-file)

Note: Apart from standard Miseq/Hiseq fastq naming extensions (R1_001_final.fastq.gz), other acceptable fastq extensions are: R1.fastq.gz/_R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, _1.fastq.gz/.1.fastq.gz.
