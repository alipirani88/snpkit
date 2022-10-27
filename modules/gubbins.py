__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *

def gubbins(gubbins_dir, input_fasta, jobrun, logger, Config):
    keep_logging('\nRunning Gubbins on input: %s\n' % input_fasta, '\nRunning Gubbins on input: %s\n' % input_fasta,
                 logger,
                 'info')
    call("module load bioperl python-anaconda2/201607 biopython dendropy reportlab fasttree RAxML fastml/gub gubbins",logger)
    gubbins_cmd = "%s --threads 6 --prefix %s/%s %s" % (
    ConfigSectionMap("gubbins", Config)['base_cmd'], gubbins_dir,
    (os.path.basename(input_fasta)).replace('.fa', ''), input_fasta)
    keep_logging('\nRunning Gubbins on: %s' % input_fasta, '\nRunning Gubbins: %s\n' % input_fasta,
                 logger,
                 'info')
    keep_logging('Running: %s' % gubbins_cmd, '%s' % gubbins_cmd, logger, 'info')
    if jobrun == "parallel-local" or jobrun == "local":
        call("cd %s" % gubbins_dir, logger)
        call(gubbins_cmd, logger)
    elif jobrun == "cluster":
        call("cd %s" % gubbins_dir, logger)
        call(gubbins_cmd, logger)
    elif jobrun == "parallel-cluster":
        job_file_name = "%s/gubbins_%s.pbs" % (gubbins_dir, os.path.basename(input_fasta))
        job_name = os.path.basename(job_file_name)
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=250:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], gubbins_dir, gubbins_cmd)
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        call("qsub %s" % job_file_name, logger)
