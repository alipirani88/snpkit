Quick Start
===========

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