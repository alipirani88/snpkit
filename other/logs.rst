.. _logs:

Log files
=========

The pipeline generates a log file following the naming convention: yyyy_mm_dd_hrs_mins_secs_analysisname.log.txt and tracks each event/command. The log file sections follow standard [Python logging conventions](https://docs.python.org/2/howto/logging.html): 

**INFO** to print STDOUT messages; 

**DEBUG** to print commands ran by pipeline, 

**ERROR** to print STDERR messages and 

**EXCEPTION** to print an exception that occured while the pipeline was running.
