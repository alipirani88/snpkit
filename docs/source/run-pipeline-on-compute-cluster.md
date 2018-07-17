Run pipeline on Compute cluster
===============================

The variant calling can be run in parallel on each sample using the -cluster argument. Set the pbs resources such as resources(-l), email(-M), queue(-q), flux_account(-A) and notification(-m) under [scheduler] section in config file. For more details, check out UMICH [flux](http://arc-ts.umich.edu/systems-and-services/flux/) website.

Possible options for -cluster option(Supported system: pbs):

***local*** : This is the default option. When this option is set, the pipeline will analyze each sample one after the another. This will not make use of multiple cores present in your system or clusters. Use this option for few number of samples or for testing purposes.

***parallel-local*** :  This option will run the pipeline and analyze samples in parallel but on a local system. This is preferred for less than 20 sample size input and multiple core local system.

***cluster*** : This option will run the pipeline on a single cluster. This option is similar to local but will rather run on a cluster node. It will not make use of multiple cores present on a cluster. Use this option for few number of samples or for testing purposes or for few number of large size samples which requires multiple cores to analyze individual samples.

***parallel-cluster*** : The pipeline is optimized for this option. When this option is set, the pipeline will run everything on a cluster making use of all the cores available. This option will also submit individual jobs for each sample seperately and whenever needed. When you set this option for each step in the pipeline, make sure all the jobs submitted by the pipeline at each step is completed before proceeding to another step. You can check the status of the job with:

```
qstat -u USERNAME
```