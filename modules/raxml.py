__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *


def raxml(tree_dir, input_fasta, jobrun, logger, Config):
    keep_logging('Running RAXML on input: %s' % input_fasta, 'Running RAXML on input: %s' % input_fasta, logger, 'info')
    #raxml_cmd = "%s/%s/%s %s -s %s -n %s_raxML" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("raxml", Config)['raxml_bin'], ConfigSectionMap("raxml", Config)['base_cmd'], ConfigSectionMap("raxml", Config)['parameters'], input_fasta, (os.path.basename(input_fasta)).replace('.fa', ''))

    raxml_cmd = "%s/%s/mpirun -np 2 %s/%s/%s -T 6 %s -s %s -n %s_raxML" % (
        ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("raxml", Config)['openmpi_bin'], ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("raxml", Config)['raxml_bin'],
    ConfigSectionMap("raxml", Config)['base_cmd'], ConfigSectionMap("raxml", Config)['parameters'], input_fasta,
    (os.path.basename(input_fasta)).replace('.fa', ''))

    keep_logging('%s' % raxml_cmd, '%s' % raxml_cmd, logger, 'info')
    if jobrun == "parallel-local" or jobrun == "local":
        call("cd %s" % tree_dir, logger)
        call(raxml_cmd, logger)
    elif jobrun == "cluster":
        call("cd %s" % tree_dir, logger)
        call(raxml_cmd, logger)
    elif jobrun == "parallel-cluster":
        job_file_name = "%s/raxml_%s.pbs" % (tree_dir, os.path.basename(input_fasta))
        job_name = os.path.basename(job_file_name)
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=12,mem=47000mb,walltime=250:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], tree_dir, raxml_cmd)
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        #os.system("qsub %s" % job_file_name)
        call("qsub %s" % job_file_name, logger)