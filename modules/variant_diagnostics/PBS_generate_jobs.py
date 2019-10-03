# System wide imports
from __future__ import division
import sys
import argparse
import re
import os
import csv
import subprocess
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
""" Hacky way to append. Instead Add this path to PYTHONPATH Variable """
from collections import OrderedDict
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
import thread
import glob
import readline
import errno
from datetime import datetime
import threading
import json
import ConfigParser
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *


def create_job(filter2_only_snp_vcf_dir, jobrun, vcf_filenames, unique_position_file, tmp_dir, Config):

    """
    This method takes the unique_position_file and list of final *_no_proximate_snp.vcf files and generates individual jobs/script.
    Each of these jobs/scripts will generate a *label file. These label file for each sample contains a field description for each position in unique_position_file.
    This field description denotes if the variant position made to the final variant list in a sample and if not then a reason/filter that caused it to filtered out from final list.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "parallel-cluster":
        """
        Supports only PBS clusters for now.
        """
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*vcf.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            call("qsub %s" % i, logger)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()

        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    elif jobrun == "cluster":
        #command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
        #os.system("bash %s" % command_file)
        command_array = []
        command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()

        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    elif jobrun == "local":
        """
        Generate a Command list of each job and run it on local system one at a time
        """

        command_array = []
        command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
        pbs_scripts = glob.glob(pbs_dir)


        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        call("bash %s" % command_file, logger)

def create_indel_job(filter2_only_snp_vcf_dir, jobrun, vcf_filenames, unique_position_file, tmp_dir, Config):

    """
    This method takes the unique_indel_position_file and list of final *_indel_final.vcf files and generates individual jobs/script.
    Each of these jobs/scripts will generate a *label file. These label file for each sample contains a field description of each position in unique_indel_position_file.
    This field description denotes if the variant position made to the final variant list in a sample and if not then a reason/filter that caused it to filtered out from final list.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "parallel-cluster":
        """
        Supports only PBS clusters for now.
        """
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s_indel.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = filter2_only_snp_vcf_dir + "/*vcf_indel.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            # os.system("qsub %s" % i)
            call("qsub %s" % i, logger)

    elif jobrun == "parallel-local" or jobrun == "cluster":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_indel_list.sh" % filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug_gatk.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s_indel.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = filter2_only_snp_vcf_dir + "/*vcf_indel.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()

        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    # elif jobrun == "cluster":
    #     command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
    #     os.system("bash %s" % command_file)
    elif jobrun == "local":
        """
        Generate a Command list of each job and run it on local system one at a time
        """

        command_array = []
        command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s_indel.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*vcf_indel.pbs"
        pbs_scripts = glob.glob(pbs_dir)


        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        call("bash %s" % command_file, logger)

def run_command(i):
    """Function to run each command and is run as a part of python Parallel mutiprocessing method.

    :param:
        i: command variable to run

    :return:
        done: string variable with completion status of command.
    """

    #call("%s" % i, logger)
    os.system("%s" % i)
    done = "Completed: %s" % i
    return done