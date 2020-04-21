[Main]
# Request to run pipeline on these samples were submitted by:
submitter: Ali
# Samples belongs to which project.
project_name: Pipeline Demo
# Date when the pipeline was run on this samples.
date: 2020-04-18
# Version of pipeline used
version: 1.2.7

[Description]
# Describe in few words why these variant sets are being generated; Why is this work being requested?
comments: Fill up this form to refelct these details in a README file.
# Any special request of changing the variant calling/filtering/other parameters. Also describe why these parameters were changed.
parameters: No changes to the default parameters. The parameters mentioned below in snitkin_filters were used for this run.
# Name of the filters applied
filters: snitkin_filters

# Default snitkin filters
[snitkin_filters]
avg_depth: no
# If AVG_DEPTH is yes, the below DP threshold will be ignored. Instead, LOW_DEPTH and HIGH_DEPTH filter parameter will be used.
# Filter variants with Depth less than the below threshold
dp: 9
# A value of 2 means that regions with less than half of the average coverage of the entire genome will fail
low_depth: 2
# A value of 5 means that regions with 5x depth greater than the average coverage will fail
high_depth: 5
# Filter variants with FQ(Consensus Quality) greater than the below threshold
fq: 0.025
fq2: 0.025
# Filter variants with MQ(Root Mean Square Quality) less than the below threshold
mq: 50
# Filter variants with Variant QUAL less than the below threshold
qual: 100
# Filter variants with GATK QualbyDepth QD parameter; filter less than the below threshold. Currently, being used for Indel SNPS only.
qd: 2.00
# Filter variants with AF1 less than the below threshold
af: 0.900
# Filter Variants that are proximate to each other within this number of range. To turn this off, use 0(zero).
prox: 0


# Reference Genome Path
[reference_genome]
name: KPNIH1
path: /scratch/esnitkin_root/esnitkin/apirani/Testing_pipelines/reference/KPNIH1/KPNIH1.fasta

