__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *


def iqtree(tree_dir, input_fasta, jobrun, logger, Config):
    keep_logging('Running iqtree on input: %s' % input_fasta, 'Running iqtree on input: %s' % input_fasta, logger, 'info')
    iqtree_cmd = "%s/%s -s %s %s -pre %s" % (
        ConfigSectionMap("iqtree", Config)['iqtree_bin'],
        ConfigSectionMap("iqtree", Config)['base_cmd'],
        input_fasta, ConfigSectionMap("iqtree", Config)['parameters'],
        (os.path.basename(input_fasta)).replace('.fa', ''))

    keep_logging('%s' % iqtree_cmd, '%s' % iqtree_cmd, logger, 'info')
    if jobrun == "parallel-local" or jobrun == "local":
        call("cd %s" % tree_dir, logger)
        call(iqtree_cmd, logger)
    elif jobrun == "cluster":
        call("cd %s" % tree_dir, logger)
        call(iqtree_cmd, logger)
    elif jobrun == "parallel-cluster":
        job_file_name = "%s/iqtree_%s.pbs" % (tree_dir, os.path.basename(input_fasta))
        job_name = os.path.basename(job_file_name)
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=12,mem=47000mb,walltime=250:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s" % (
        job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'],
        ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], tree_dir,
        iqtree_cmd)
        f1 = open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        # os.system("qsub %s" % job_file_name)
        keep_logging('qsub %s' % job_file_name, 'qsub %s' % job_file_name, logger, 'info')
	#call("qsub %s" % job_file_name, logger)
