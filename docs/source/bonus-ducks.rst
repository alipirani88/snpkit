.. _bonus-ducks:


Bonus Ducks
===========

Run only a part of variant calling step
---------------------------------------

OPTIONAL: In case, you want to rerun the above analysis with different variant call/filter parameters (i.e skip the read cleaning, alignment and post-alignment steps), you can do it as follows:

NOTE: OPTIONAL

::

	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps varcall,filter,stats -cluster parallel-cluster


Run variant calling step on selected samples provided in filenames.txt. 
-----------------------------------------------------------------------


filenames.txt should contain fastq filenames with one single-end filename per line.

::

	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps varcall,filter,stats -cluster parallel-cluster -filenames filenames.txt

Run core_prep and core steps on selected files
----------------------------------------------


If you decide to exclude some samples from core snp analysis step (step 3), there is a way to do that without running the entire pipeline. 


::

	python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/variant_call.py -type PE -readsdir /Path-To-Your/test_readsdir/ -outdir /Path/test_output_core/ -analysis output_prefix -index MRSA_USA_300 -steps core -cluster parallel-cluster -filenames filenames_custom



