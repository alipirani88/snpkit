__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *

def fasttree(tree_dir, input_fasta, cluster, logger, Config):
    keep_logging(' - Running Fasttree on input: %s' % input_fasta, 'Running Fasttree on input: %s' % input_fasta, logger, 'info')
    fasttree_cmd = "%s -nt %s > %s/%s_FastTree.tree" % (ConfigSectionMap("fasttree", Config)['base_cmd'], input_fasta, tree_dir, (os.path.basename(input_fasta)).replace('.fa', ''))
    keep_logging(' - %s' % fasttree_cmd, '%s' % fasttree_cmd, logger, 'info')
    if cluster == "parallel-local" or cluster == "local":
        call("cd %s" % tree_dir, logger)
        call(fasttree_cmd, logger)
    elif cluster == "cluster":
        call("cd %s" % tree_dir, logger)
        call(fasttree_cmd, logger)
    elif cluster == "parallel-cluster":
        job_file_name = "%s/fasttree_%s.pbs" % (tree_dir, os.path.basename(input_fasta))
        job_name = os.path.basename(job_file_name)
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=76:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], tree_dir, fasttree_cmd)
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        call("qsub %s" % job_file_name, logger)