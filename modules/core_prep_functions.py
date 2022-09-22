# System wide imports
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
#import pandas as pd
import errno
from datetime import datetime
import time
import threading
import json
import ConfigParser
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
from core_prep_sanity_checks import *
from core_prep_functions import *


"""core_prep methods 

    This block contains methods that are respnsible for running the first part of core_All step of the pipeline.
    This methods generates all the necessary intermediate files required for the second part of core_All step.
    Example of intermediate files: various diagnostics files/matrices where it decides why a variant was filtered out.

"""


def create_positions_filestep(vcf_filenames, temp_dir, outgroup, logger, filter2_only_snp_vcf_dir, numcores):
    """
    This method gathers SNP positions from each final *_no_proximate_snp.vcf file (these are the positions that passed variant filter parameters
    from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files. Use these *_no_proximate_snp.vcf_position files to generate a list of unique_position_file
    :param: list of final vcf filenames i.e *.vcf_no_proximate_snp.vcf . These files are the final output of variant calling step for each sample.
    :return: unique_position_file
    """

    filter2_only_snp_position_files_array = []
    for file in vcf_filenames:
        with open(file, 'rU') as csv_file:
            file_name = temp_dir + "/" + os.path.basename(file) + "_positions"
            addpositionfilenametoarray = file_name
            filter2_only_snp_position_files_array.append(addpositionfilenametoarray)
            f1 = open(file_name, 'w+')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                position = row[0]
                if not position.startswith('#'):
                    p_string = row[1] + "\n"
                    f1.write(p_string)
            f1.close()
        csv_file.close()

    """ Get Positions Specific to Outgroup Sample name """
    if outgroup:
        outgroup_vcf_filename = str(outgroup.replace('R1_001.fastq.gz', '')) + "_filter2_final.vcf_no_proximate_snp.vcf"
        outgroup_indel_vcf_filename = str(outgroup.replace('R1_001.fastq.gz', '')) + "_filter2_indel_final.vcf"
        outgroup_position_file_name = temp_dir + "/" + outgroup_vcf_filename + "_positions"
        outgroup_position_array = []
        f1 = open(outgroup_position_file_name, 'r+')
        for lines in f1:
            lines = lines.strip()
            outgroup_position_array.append(int(lines))
        f1.close()

        position_array_excluding_outgroup = []
        for filess in filter2_only_snp_position_files_array:
            if outgroup not in filess:
                f = open(filess, 'r+')
                for line in f:
                    line = line.strip()
                    position_array_excluding_outgroup.append(int(line))
                f.close()
        position_array_unique_excluding_outgroup = set(position_array_excluding_outgroup)
        position_array_sort_excluding_outgroup = sorted(position_array_unique_excluding_outgroup)
        outgroup_specific_positions = []
        f_outgroup = open("%s/outgroup_specific_positions.txt" % filter2_only_snp_vcf_dir, 'w+')
        for i in outgroup_position_array:
            if i not in position_array_sort_excluding_outgroup:
                f_outgroup.write(str(i) + '\n')
                outgroup_specific_positions.append(int(i))
        f_outgroup.close()

        # Print Checks
        print "No. of variant positions in outgroup: %s" % len(outgroup_position_array)
        print "No. of variant positions specific to outgroup: %s" % len(outgroup_specific_positions)

        position_array = []
        for filess in filter2_only_snp_position_files_array:
            f = open(filess, 'r+')
            for line in f:
                line = line.strip()
                # Changed variable to suit sorting: 25-07-2018
                position_array.append(int(line))
            f.close()

        keep_logging('- Sorting unique variant positions.', '- Sorting unique variant positions.', logger, 'info')
        position_array_unique = set(position_array)
        position_array_sort = sorted(position_array_unique)
        keep_logging('\nThe number of unique variant positions:%s' % len(position_array_sort),
                     '\nThe number of unique variant positions:%s' % len(position_array_sort), logger, 'info')
        unique_position_file = "%s/unique_positions_file" % filter2_only_snp_vcf_dir
        f = open(unique_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()

        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?',
                         'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
            exit()

        return unique_position_file

    else:

        """ Create position array containing unique positiones from positions file """
        position_array = []
        for filess in filter2_only_snp_position_files_array:
            f = open(filess, 'r+')
            for line in f:
                line = line.strip()
                position_array.append(int(line))
            f.close()

        keep_logging('- Sorting unique variant positions.', '- Sorting unique variant positions.', logger, 'info')
        position_array_unique = set(position_array)
        position_array_sort = sorted(position_array_unique)
        keep_logging('- The number of unique variant positions:%s' % len(position_array_sort),
                     '- The number of unique variant positions:%s' % len(position_array_sort), logger, 'info')
        unique_position_file = "%s/unique_positions_file" % filter2_only_snp_vcf_dir
        f = open(unique_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()

        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?',
                         'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
            exit()
        return unique_position_file

def create_indel_positions_filestep(vcf_filenames, temp_dir, outgroup, logger, filter2_only_snp_vcf_dir, numcores):
    """
    This function gathers Indel positions from each final *_indel_final.vcf (these are the positions that passed variant filter parameters
    from variant calling pipeline) and write to *_indel_final.vcf files. Use these *_indel_final.vcf_position files to generate a list of unique_position_file
    :param: list of final vcf filenames i.e *_indel_final.vcf . These files are the final output of variant calling step for each sample.
    :return: unique_indel_position_file
    """

    filter2_only_indel_position_files_array = []
    for file in vcf_filenames:
        indel_file = file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf')
        with open(indel_file, 'rU') as csv_file:
            file_name = temp_dir + "/" + os.path.basename(indel_file) + "_positions"
            addpositionfilenametoarray = file_name
            filter2_only_indel_position_files_array.append(addpositionfilenametoarray)
            f1 = open(file_name, 'w+')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                position = row[0]
                if not position.startswith('#'):
                    p_string = row[1] + "\n"
                    f1.write(p_string)
            f1.close()
        csv_file.close()

    """ Get Positions Specific to Outgroup Sample name """
    if outgroup:
        outgroup_indel_vcf_filename = str(outgroup.replace('R1_001.fastq.gz', '')) + "_filter2_indel_final.vcf"
        outgroup_position_indel_file_name = temp_dir + "/" + outgroup_indel_vcf_filename + "_positions"
        print outgroup_position_indel_file_name
        outgroup_position_indel_array = []
        f1 = open(outgroup_position_indel_file_name, 'r+')
        for lines in f1:
            lines = lines.strip()
            outgroup_position_indel_array.append(int(lines))
        f1.close()

        position_array_indel_excluding_outgroup = []
        for filess in filter2_only_indel_position_files_array:
            if outgroup not in filess:
                f = open(filess, 'r+')
                for line in f:
                    line = line.strip()
                    position_array_indel_excluding_outgroup.append(int(line))
                f.close()
        position_array_indel_unique_excluding_outgroup = set(position_array_indel_excluding_outgroup)
        position_array_sort_indel_excluding_outgroup = sorted(position_array_indel_unique_excluding_outgroup)
        outgroup_indel_specific_positions = []
        f_outgroup = open("%s/outgroup_indel_specific_positions.txt" % filter2_only_snp_vcf_dir, 'w+')
        for i in outgroup_position_indel_array:
            if i not in position_array_sort_indel_excluding_outgroup:
                f_outgroup.write(str(i) + '\n')
                outgroup_indel_specific_positions.append(int(i))
        f_outgroup.close()

        # Print Checks
        print "No. of indel variant positions in outgroup: %s" % len(outgroup_position_indel_array)
        print "No. of indel variant positions specific to outgroup: %s" % len(outgroup_indel_specific_positions)

        position_array = []
        for filess in filter2_only_indel_position_files_array:
            f = open(filess, 'r+')
            for line in f:
                line = line.strip()
                # Changed variable to suit sorting: 25-07-2018
                position_array.append(int(line))
            f.close()
        position_array_unique = set(position_array)
        position_array_sort = sorted(position_array_unique)
        keep_logging('- The number of unique indel positions:%s' % len(position_array_sort),
                     '- The number of unique indel positions:%s' % len(position_array_sort), logger, 'info')
        unique_indel_position_file = "%s/unique_indel_positions_file" % filter2_only_snp_vcf_dir
        f = open(unique_indel_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()
        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?',
                         'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
            exit()
        return unique_indel_position_file


    else:

        """ Create position array containing unique positiones from positions file """
        position_array = []
        for filess in filter2_only_indel_position_files_array:
            f = open(filess, 'r+')
            for line in f:
                line = line.strip()
                # Changed variable to suit sorting: 25-07-2018
                position_array.append(int(line))
            f.close()
        position_array_unique = set(position_array)
        position_array_sort = sorted(position_array_unique)
        keep_logging('- The number of unique indel positions:%s' % len(position_array_sort),
                     '- The number of unique indel positions:%s' % len(position_array_sort), logger, 'info')
        unique_indel_position_file = "%s/unique_indel_positions_file" % filter2_only_snp_vcf_dir
        f = open(unique_indel_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()
        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?',
                         'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
            exit()
        return unique_indel_position_file

def create_job(jobrun, vcf_filenames, unique_position_file, tmp_dir, scheduler_directives, script_Directive,
               job_name_flag, temp_dir, outgroup, logger, filter2_only_snp_vcf_dir, numcores):
    """
    This method takes the unique_position_file and list of final *_no_proximate_snp.vcf files and generates individual jobs/script.
    Each of these jobs/scripts will generate a *label file. These label file for each sample contains a field description for each position in unique_position_file.
    This field description denotes if the variant position made to the final variant list in a sample and if not then a reason/filter that caused it to filtered out from final list.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "parallel-cluster":
        """ Deprecated """
        """
        Supports only PBS clusters for now.
        """
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\npython %s/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            os.path.dirname(os.path.abspath(__file__)), filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        # os.system("mv %s/*.pbs %s/temp" % (filter2_only_snp_vcf_dir, filter2_only_snp_vcf_dir))
        pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            call("qsub %s" % i, logger)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (
            os.path.dirname(os.path.abspath(__file__)), filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s.pbs" % (i)

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

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
        if numcores:
            num_cores = int(numcores)
        else:
            # Great Lakes Integration here.
            if args.scheduler == "SLURM":
                proc = subprocess.Popen(["echo $SLURM_CPUS_PER_TASK"], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                num_cores = int(out.strip())
            elif args.scheduler == "PBS":
                num_cores = multiprocessing.cpu_count()
            else:
                num_cores = 1

        
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    elif jobrun == "cluster":
        command_array = []
        command_file = "%s/commands_list.sh" % temp_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), filter2_only_snp_vcf_dir, i, unique_position_file,
                tmp_dir)
            job_file_name = "%s.sbat" % (i)

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()
            os.system("mv %s %s" % (job_file_name, temp_dir))
        pbs_dir = "%s/*vcf.sbat" % temp_dir
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        if numcores:
            num_cores = int(numcores)
        else:
            # Slurm Changes here.
            if args.scheduler == "SLURM":
                proc = subprocess.Popen(["echo $SLURM_CPUS_PER_TASK"], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                num_cores = int(out.strip())
            elif args.scheduler == "PBS":
                num_cores = multiprocessing.cpu_count()
            else:
                num_cores = 1

        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    elif jobrun == "local":
        """
        Generate a Command list of each job and run it on local system one at a time
        """

        command_array = []
        command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), filter2_only_snp_vcf_dir, i, unique_position_file,
                tmp_dir)
            job_file_name = "%s.pbs" % (i)

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

        # os.system("mv %s/*.pbs %s/temp" % (filter2_only_snp_vcf_dir, filter2_only_snp_vcf_dir))
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

def create_indel_job(jobrun, vcf_filenames, unique_position_file, tmp_dir, scheduler_directives, script_Directive,
                     job_name_flag, temp_dir, outgroup, logger, filter2_only_snp_vcf_dir, numcores):
    """
    This method takes the unique_indel_position_file and list of final *_indel_final.vcf files and generates individual jobs/script.
    Each of these jobs/scripts will generate a *label file. These label file for each sample contains a field description of each position in unique_indel_position_file.
    This field description denotes if the variant position made to the final variant list in a sample and if not then a reason/filter that caused it to filtered out from final list.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "parallel-cluster":
        """ Deprecated """
        """
        Supports only PBS clusters for now.
        """
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\npython %s/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            os.path.dirname(os.path.abspath(__file__)), filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s_indel.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        # os.system("mv %s/*.pbs %s/temp" % (filter2_only_snp_vcf_dir, filter2_only_snp_vcf_dir))
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
        command_file = "%s/commands_indel_list.sh" % temp_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/variant_diagnostics/reason_job_indel_debug_gatk.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), filter2_only_snp_vcf_dir, i, unique_position_file,
                tmp_dir)
            job_file_name = "%s_indel.sbat" % (i)

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()
            os.system("mv %s %s" % (job_file_name, temp_dir))

        pbs_dir = "%s/*vcf_indel.sbat" % temp_dir
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        if numcores:
            num_cores = int(numcores)
        else:
            # Great Lakes Integration here.
            if args.scheduler == "SLURM":
                proc = subprocess.Popen(["echo $SLURM_CPUS_PER_TASK"], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                num_cores = int(out.strip())
            elif args.scheduler == "PBS":
                num_cores = multiprocessing.cpu_count()
            else:
                num_cores = 1

        
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    elif jobrun == "local":
        """
        Generate a Command list of each job and run it on local system one at a time
        """

        command_array = []
        command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')

        # Great Lakes Integration here
        for i in vcf_filenames:
            command = "python %s/variant_diagnostics/reason_job_indel_debug_gatk.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), filter2_only_snp_vcf_dir, i, unique_position_file,
                tmp_dir)
            job_file_name = "%s_indel.pbs" % (i)

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

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
        call("bash %s" % command_file, logger)

def run_command(i):
    """Function to run each command and is run as a part of python Parallel mutiprocessing method.

    :param:
        i: command variable to run

    :return:
        done: string variable with completion status of command.
    """

    #call("%s" % i, logger)
    os.system(i)
    done = "Completed: %s" % i
    return done