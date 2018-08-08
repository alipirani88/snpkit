from __future__ import division
import sys
import argparse
import re
import os
import csv
import subprocess
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
""" Hacky way yo append. Instead Add this path to PYTHONPATH Variable """
from collections import OrderedDict
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
import thread
import glob
import readline
import pandas as pd
import errno
from pyfasta import Fasta
from datetime import datetime
import threading
import json
from cyvcf2 import VCF
import ConfigParser
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
from tabix import *
from Bio import SeqIO
from phage_detection import *
from find_repeats import *
from mask_regions import *
from fasttree import fasttree
from gubbins import *
from raxml import raxml
from pyfasta import Fasta
from core_prep_sanity_checks import *
#from core_prep import *


# Parse Command line Arguments
parser = argparse.ArgumentParser(description='Parsing filtered VCF files and investigating Variants to determine the reason why it was filtered out from the final list')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
required.add_argument('-filter2_only_snp_vcf_filenames', action='store', dest="filter2_only_snp_vcf_filenames",
                    help='Names of filter2 only SNP vcf files with name per line.')
optional.add_argument('-jobrun', action='store', dest="jobrun",
                    help='Running a job on Cluster, Running Parallel jobs, Run jobs/commands locally (default): cluster, local, parallel-local, parallel-single-cluster')
optional.add_argument('-cluster_type', action='store', dest="cluster_type",
                    help='Type of Cluster: torque, pbs, sgd')
optional.add_argument('-cluster_resources', action='store', dest="cluster_resources",
                    help='Cluster Resources to use. for example nodes,core. Ex: 1,4')
optional.add_argument('-numcores', action='store', dest="numcores",
                    help='Number of cores to use on local system for parallel-local parameter')
optional.add_argument('-remove_temp', action='store', dest="remove_temp",
                    help='Remove Temporary files generated during the run')
optional.add_argument('-gubbins', action='store', dest="gubbins", help='yes/no for running gubbins')
required.add_argument('-reference', action='store', dest="reference",
                    help='Path to Reference Fasta file for consensus generation')
required.add_argument('-steps', action='store', dest="steps",
                    help='Analysis Steps to be performed. This should be in sequential order.'
                         'Step 1: Run pbs jobs and process all pipeline generated vcf files to generate label files'
                         'Step 2: Analyze label files and generate matrix'
                         'Step 3: DP/FQ Analysis')
required.add_argument('-results_dir', action='store', dest="results_dir",
                    help='Path to Core results directory')
required.add_argument('-config', action='store', dest="config",
                    help='Path to config file')
# optional.add_argument('-db', action='store', dest="snpeff_db",
#                     help='snpEff prebuilt reference database to use for variant annotations. The database will be downloaded in /data/ folder under snpEff install directory. Make sure if you are providing the name of pre-built snpEff reference database then the build option of snpeff section in config section is set to \"no\"')
optional.add_argument('-debug_mode', action='store', dest="debug_mode",
                    help='yes/no for debug mode')
args = parser.parse_args()

def make_sure_path_exists(out_path):
    """
    Fuction to make sure output folder exists. If not, create it.
    :param: out_path
    :return: null/exception
    """
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('\nErrors in output folder path! please change the output path or analysis name\n',
                         '\nErrors in output folder path! please change the output path or analysis name\n', logger,
                         'info')
            exit()

""" Core Prep Methods"""
def run_command(i):
    """
    Function to run each command and is run as a part of python Parallel mutiprocessing method.
    :param: command
    :return: Done status
    """
    call("%s" % i, logger)
    done = "Completed: %s" % i
    return done

""" core_prep methods """
def create_positions_filestep(vcf_filenames):

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

    """ Create position array containing unique positiones from positions file """

    position_array = []
    for filess in filter2_only_snp_position_files_array:
        f = open(filess, 'r+')
        for line in f:
            line = line.strip()
            # Changed variable to suit sorting: 25-07-2018
            position_array.append(int(line))
        f.close()
    # Check why python sorting is not working
    keep_logging('Sorting unique variant positions.\n', 'Sorting unique variant positions.\n', logger, 'info')
    position_array_unique = set(position_array)
    position_array_sort = sorted(position_array_unique)
    keep_logging('\nThe number of unique variant positions:%s' % len(position_array_sort), '\nThe number of unique variant positions:%s' % len(position_array_sort), logger, 'info')
    unique_position_file = "%s/unique_positions_file" % args.filter2_only_snp_vcf_dir
    f=open(unique_position_file, 'w+')
    for i in position_array_sort:
        # Changed variable to suit sorting: 25-07-2018
        f.write(str(i) + "\n")
    f.close()
    if len(position_array_sort) == 0:
        keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
        exit()
    return unique_position_file

def create_indel_positions_filestep(vcf_filenames):

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
    keep_logging('\nThe number of unique indel positions:%s' % len(position_array_sort), '\nThe number of unique indel positions:%s' % len(position_array_sort), logger, 'info')
    unique_indel_position_file = "%s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir
    f=open(unique_indel_position_file, 'w+')
    for i in position_array_sort:
        # Changed variable to suit sorting: 25-07-2018
        f.write(str(i) + "\n")
    f.close()
    if len(position_array_sort) == 0:
        keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
        exit()
    return unique_indel_position_file

def create_job(jobrun, vcf_filenames, unique_position_file, tmp_dir):

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
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*vcf.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        if args.numcores:
            num_cores = int(num_cores)
        else:
            num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    elif jobrun == "cluster":
        #command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
        #os.system("bash %s" % command_file)
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
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*vcf.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        if args.numcores:
            num_cores = int(num_cores)
        else:
            num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    elif jobrun == "local":
        """
        Generate a Command list of each job and run it on local system one at a time
        """

        command_array = []
        command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*vcf.pbs"
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

def create_indel_job(jobrun, vcf_filenames, unique_position_file, tmp_dir):

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
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
            job_file_name = "%s_indel.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*vcf_indel.pbs"
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
        command_file = "%s/commands_indel_list.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
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
        if args.numcores:
            num_cores = int(num_cores)
        else:
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


""" core methods """

def generate_paste_command():
    """
    This Function will take all the *label file and generate/paste it column wise to generate a matrix. These matrix will be used in downstream analysis.
    :param: null
    :return: null
    """

    """ Paste/Generate and sort SNP Filter Label Matrix """
    paste_file = args.filter2_only_snp_vcf_dir + "/paste_label_files.sh"
    f4=open(paste_file, 'w+')
    paste_command = "paste %s/unique_positions_file" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
        paste_command = paste_command + " " + label_file
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

    call("%s" % header_awk_cmd, logger)
    call("%s" % sed_header, logger)
    call("%s" % sed_header_2, logger)

    temp_paste_command = paste_command + " > %s/temp_label_final_raw.txt" % args.filter2_only_snp_vcf_dir
    paste_command = paste_command + " > %s/All_label_final_raw" % args.filter2_only_snp_vcf_dir
    f4.write(paste_command)
    f4.close()
    sort_All_label_cmd = "sort -n -k1,1 %s/All_label_final_raw > %s/All_label_final_sorted.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    paste_command_header = "cat %s/header.txt %s/All_label_final_sorted.txt > %s/All_label_final_sorted_header.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    ls = []
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
        ls.append(label_file)
    ls.insert(0, "%s/unique_positions_file" % args.filter2_only_snp_vcf_dir)

    with open('%s/All_label_final_raw.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
        outfile.write(paste_command)
    outfile.close()

    with open('%s/temp_label_final_raw.txt.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
        outfile.write(temp_paste_command)
    outfile.close()
    call("bash %s/All_label_final_raw.sh" % args.filter2_only_snp_vcf_dir, logger)
    call("bash %s/temp_label_final_raw.txt.sh" % args.filter2_only_snp_vcf_dir, logger)


    """
    remove this lines
    #subprocess.call(["%s" % paste_command], shell=True)
    #subprocess.call(["%s" % temp_paste_command], shell=True)
    #subprocess.check_call('%s' % paste_command)
    #subprocess.check_call('%s' % temp_paste_command)
    #os.system(paste_command) change
    #os.system(temp_paste_command) change
    """

    call("%s" % sort_All_label_cmd, logger)
    call("%s" % paste_command_header, logger)

    """ Assign numeric code to each variant filter reason"""
    subprocess.call(["sed -i 's/reference_unmapped_position/0/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/reference_allele/1/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/VARIANT/1TRUE/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL_DP_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP_QUAL_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL_DP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP_QUAL/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_proximate_SNP/7/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL_DP/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP_QUAL/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ/5/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ/6/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    remove_unwanted_text = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir
    call("%s" % remove_unwanted_text, logger)

def generate_indel_paste_command():
    """
    This Function will take all the *label file and generate/paste it column wise to generate a matrix. These matrix will be used in downstream analysis.
    :param: null
    :return: null
    """

    """ Paste/Generate and sort SNP Filter Label Matrix """
    paste_file = args.filter2_only_snp_vcf_dir + "/paste_indel_label_files.sh"
    f4=open(paste_file, 'w+')
    paste_command = "paste %s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf_indel_positions_label')
        paste_command = paste_command + " " + label_file
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

    #os.system(header_awk_cmd)
    #os.system(sed_header)
    #os.system(sed_header_2)

    call("%s" % header_awk_cmd, logger)
    call("%s" % sed_header, logger)
    call("%s" % sed_header_2, logger)



    temp_paste_command = paste_command + " > %s/temp_indel_label_final_raw.txt" % args.filter2_only_snp_vcf_dir
    paste_command = paste_command + " > %s/All_indel_label_final_raw" % args.filter2_only_snp_vcf_dir
    f4.write(paste_command)
    f4.close()

    call("bash %s" % paste_file, logger)

    sort_All_label_cmd = "sort -n -k1,1 %s/All_indel_label_final_raw > %s/All_indel_label_final_sorted.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    paste_command_header = "cat %s/header.txt %s/All_indel_label_final_sorted.txt > %s/All_indel_label_final_sorted_header.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    ls = []
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf_indel_positions_label')
        ls.append(label_file)
    ls.insert(0, "%s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir)

    with open('%s/All_indel_label_final_raw.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile2:
        outfile2.write(paste_command)
    outfile2.close()

    with open('%s/temp_indel_label_final_raw.txt.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile2:
        outfile2.write(temp_paste_command)
    outfile2.close()

    # Why is this not working?
    call("bash %s/All_indel_label_final_raw.sh" % args.filter2_only_snp_vcf_dir, logger)
    call("bash %s/temp_indel_label_final_raw.txt.sh" % args.filter2_only_snp_vcf_dir, logger)
    keep_logging('Finished pasting...DONE', 'Finished pasting...DONE', logger, 'info')

    """
    remove this lines
    #subprocess.call(["%s" % paste_command], shell=True)
    #subprocess.call(["%s" % temp_paste_command], shell=True)
    #subprocess.check_call('%s' % paste_command)
    #subprocess.check_call('%s' % temp_paste_command)
    #os.system(paste_command) change
    #os.system(temp_paste_command) change
    """

    call("%s" % sort_All_label_cmd, logger)
    call("%s" % paste_command_header, logger)

    """ Assign numeric code to each variant filter reason"""
    subprocess.call(["sed -i 's/reference_unmapped_position/0/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/reference_allele/1/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/VARIANT/1TRUE/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL_DP_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP_QUAL_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL_DP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP_QUAL/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_proximate_SNP/7/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL_DP/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP_QUAL/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/LowFQ/5/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    subprocess.call(["sed -i 's/HighFQ/6/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir], shell=True)
    remove_unwanted_text = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir
    call("%s" % remove_unwanted_text, logger)

def generate_position_label_data_matrix():

    """
    Generate different list of Positions using the matrix All_label_final_sorted_header.txt.

    (Defining Core Variant Position: Variant Position which was not filtered out in any of the other samples due to variant filter parameter and also this position was present in all the samples(not unmapped)).

    Filtered Position label matrix:
        List of non-core positions. These positions didn't make it to the final core list because it was filtered out in one of the samples.

    Only_ref_variant_positions_for_closely_matrix.txt :
        Those Positions where the variant was either reference allele or a variant that passed all the variant filter parameters.

    :param: null
    :return: null

    """
    def generate_position_label_data_matrix_All_label():
        position_label = OrderedDict()
        f1=open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'w+')
        f2=open("%s/Only_ref_variant_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f3=open("%s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f4=open("%s/Only_filtered_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            keep_logging('Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir, 'Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir, logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                position_label[row[0]] = row[1:]
            keep_logging('Generating different list of Positions and heatmap data matrix... \n', 'Generating different list of Positions and heatmap data matrix... \n', logger, 'info')
            print_string_header = "\t"
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"
            f2.write('\t' + print_string_header.strip() + '\n')
            f3.write('\t' + print_string_header.strip() + '\n')
            f4.write('\t' + print_string_header.strip() + '\n')
            for value in position_label:
                lll = ['0', '2', '3', '4', '5', '6', '7']
                ref_var = ['1', '1TRUE']
                if set(ref_var) & set(position_label[value]):
                    if set(lll) & set(position_label[value]):
                        print_string = ""
                        for i in position_label[value]:
                            print_string = print_string + "\t" + i
                        STRR2 = value + print_string + "\n"
                        f3.write(STRR2)
                        if position_label[value].count('1TRUE') >= 2:
                            f4.write('1\n')
                        else:
                            f4.write('0\n')
                    else:
                        strr = value + "\n"
                        f1.write(strr)
                        STRR3 =	value +	"\t" + str(position_label[value]) + "\n"
                        f2.write(STRR3)
        csv_file.close()
        f1.close()
        f2.close()
        f3.close()
        f4.close()
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/1TRUE/-1/g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)

    def temp_generate_position_label_data_matrix_All_label():

        """
        Read temp_label_final_raw.txt SNP position label data matrix for generating barplot statistics.
        """
        temp_position_label = OrderedDict()
        f33=open("%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        print_string_header = "\t"
        for i in vcf_filenames:
            print_string_header = print_string_header + os.path.basename(i) + "\t"
        f33.write('\t' + print_string_header.strip() + '\n')
        keep_logging('Reading temporary label positions file: %s/temp_label_final_raw.txt \n' % args.filter2_only_snp_vcf_dir, 'Reading temporary label positions file: %s/temp_label_final_raw.txt \n' % args.filter2_only_snp_vcf_dir, logger, 'info')
        lll = ['reference_unmapped_position', 'LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP', 'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ', 'LowFQ_proximate_SNP', 'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP', 'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP', 'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
        ref_var = ['reference_allele', 'VARIANT']
        with open("%s/temp_label_final_raw.txt" % args.filter2_only_snp_vcf_dir, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                if set(ref_var) & set(row[1:]):
                    if set(lll) & set(row[1:]):
                        print_string = ""
                        for i in row[1:]:
                            print_string = print_string + "\t" + i
                        STRR2 = row[0] + print_string + "\n"
                        f33.write(STRR2)
        csv_file.close()
        f33.close()
        """
        Read temp_Only_filtered_positions_for_closely_matrix file and generate a matrix of positions that are being filtered just because of FQ
        """
        temp_position_label_FQ = OrderedDict()
        f44=open("%s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            keep_logging('Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir, 'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir, logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)

            for row in csv_reader:
              temp_position_label_FQ[row[0]] = row[1:]
            print_string_header = "\t"
            for i in vcf_filenames:
              print_string_header = print_string_header + os.path.basename(i) + "\t"
            f44.write('\t' + print_string_header.strip() + '\n')
            for value in temp_position_label_FQ:
              lll = ['LowFQ']
              if set(lll) & set(temp_position_label_FQ[value]):
                  print_string = ""
                  for i in temp_position_label_FQ[value]:
                      print_string = print_string + "\t" + i
                  STRR2 = value + print_string + "\n"
                  f44.write(STRR2)
            f44.close()
            csv_file.close()
            f44.close()

        """
        Perform Sed on temp files. Find a faster way to do this.
        """
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ/3/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)


        """
        Read temp_Only_filtered_positions_for_closely_matrix file and generate a matrix of positions that are being filtered just because of Dp
        """
        temp_position_label_DP = OrderedDict()
        f44=open("%s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            keep_logging('Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir, 'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir, logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                temp_position_label_DP[row[0]] = row[1:]
            print_string_header = "\t"
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"
            f44.write('\t' + print_string_header.strip() + '\n')
            for value in temp_position_label_DP:
                lll = ['HighFQ_DP']
                ref_var = ['reference_allele', 'VARIANT']
                if set(lll) & set(temp_position_label_FQ[value]):
                    print_string = ""
                    for i in temp_position_label_FQ[value]:
                        print_string = print_string + "\t" + i
                    STRR2 = value + print_string + "\n"
                    f44.write(STRR2)
        f44.close()
        csv_file.close()

        """
        Perform Sed on temp files. Find a faster way to do this.
        """
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP/3/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)


    def barplot_stats():
        keep_logging('\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n', '\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n', logger, 'info')
        """
        Read each Sample columns and calculate the percentage of each label to generate barplot statistics.
        This will give a visual explanation of how many positions in each samples were filtered out because of different reason
        """

        c_reader = csv.reader(open('%s/temp_Only_filtered_positions_for_closely_matrix.txt' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
        columns = list(zip(*c_reader))
        keep_logging('Finished reading columns...', 'Finished reading columns...', logger, 'info')
        counts = 1
        end = len(vcf_filenames) + 1
        f_bar_count = open("%s/bargraph_counts.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_perc = open("%s/bargraph_percentage.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_count.write("Sample\tunmapped_positions\treference_allele\ttrue_variant\tOnly_low_FQ\tOnly_DP\tOnly_low_MQ\tother\n")
        f_bar_perc.write("Sample\tunmapped_positions_perc\ttrue_variant_perc\tOnly_low_FQ_perc\tOnly_DP_perc\tOnly_low_MQ_perc\tother_perc\n")
        for i in xrange(1, end, 1):
            """ Bar Count Statistics: Variant Position Count Statistics """
            true_variant = columns[i].count('VARIANT')
            unmapped_positions = columns[i].count('reference_unmapped_position')
            reference_allele = columns[i].count('reference_allele')
            Only_low_FQ = columns[i].count('LowFQ')
            Only_DP = columns[i].count('HighFQ_DP')
            Only_low_MQ = columns[i].count('HighFQ')
            low_FQ_other_parameters = columns[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns[i].count('LowFQ_DP_QUAL_proximate_SNP') + columns[i].count('LowFQ_QUAL_proximate_SNP') + columns[i].count('LowFQ_DP_proximate_SNP') + columns[i].count('LowFQ_proximate_SNP') + columns[i].count('LowFQ_QUAL_DP') + columns[i].count('LowFQ_DP_QUAL') + columns[i].count('LowFQ_QUAL') + columns[i].count('LowFQ_DP')
            high_FQ_other_parameters = columns[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns[i].count('HighFQ_DP_QUAL_proximate_SNP') + columns[i].count('HighFQ_QUAL_proximate_SNP') + columns[i].count('HighFQ_DP_proximate_SNP') + columns[i].count('HighFQ_proximate_SNP') + columns[i].count('HighFQ_QUAL_DP') + columns[i].count('HighFQ_DP_QUAL') + columns[i].count('HighFQ_QUAL')
            other = low_FQ_other_parameters + high_FQ_other_parameters
            total = true_variant + unmapped_positions + reference_allele + Only_low_FQ + Only_DP + low_FQ_other_parameters + high_FQ_other_parameters + Only_low_MQ
            filename_count = i - 1
            bar_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), unmapped_positions, reference_allele, true_variant, Only_low_FQ, Only_DP, Only_low_MQ, other)
            f_bar_count.write(bar_string)

            """ Bar Count Percentage Statistics: Variant Position Percentage Statistics """
            try:
                true_variant_perc = float((columns[i].count('VARIANT') * 100) / total)
            except ZeroDivisionError:
                true_variant_perc = 0
            try:
                unmapped_positions_perc = float((columns[i].count('reference_unmapped_position') * 100) / total)
            except ZeroDivisionError:
                unmapped_positions_perc = 0
            try:
                reference_allele_perc = float((columns[i].count('reference_allele') * 100) / total)
            except ZeroDivisionError:
                reference_allele_perc = 0
            try:
                Only_low_FQ_perc = float((columns[i].count('LowFQ') * 100) / total)
            except ZeroDivisionError:
                Only_low_FQ_perc = 0
            try:
                Only_DP_perc = float((columns[i].count('HighFQ_DP') * 100) / total)
            except ZeroDivisionError:
                Only_DP_perc = 0
            try:
                Only_low_MQ_perc = float((columns[i].count('HighFQ') * 100) / total)
            except ZeroDivisionError:
                Only_low_MQ_perc = 0
            try:
                low_FQ_other_parameters_perc = float(((columns[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns[i].count('LowFQ_DP_QUAL_proximate_SNP') + columns[i].count('LowFQ_QUAL_proximate_SNP') + columns[i].count('LowFQ_DP_proximate_SNP') + columns[i].count('LowFQ_proximate_SNP') + columns[i].count('LowFQ_QUAL_DP') + columns[i].count('LowFQ_DP_QUAL') + columns[i].count('LowFQ_QUAL') + columns[i].count('LowFQ_DP'))  * 100) / total)
            except ZeroDivisionError:
                low_FQ_other_parameters_perc = 0
            try:
                high_FQ_other_parameters_perc = float(((columns[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns[i].count('HighFQ_DP_QUAL_proximate_SNP') + columns[i].count('HighFQ_QUAL_proximate_SNP') + columns[i].count('HighFQ_DP_proximate_SNP') + columns[i].count('HighFQ_proximate_SNP') + columns[i].count('HighFQ_QUAL_DP') + columns[i].count('HighFQ_DP_QUAL') + columns[i].count('HighFQ_QUAL')) * 100) / total)
            except ZeroDivisionError:
                high_FQ_other_parameters_perc = 0

            other_perc = float(low_FQ_other_parameters_perc + high_FQ_other_parameters_perc)
            bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), unmapped_positions_perc, true_variant_perc, Only_low_FQ_perc, Only_DP_perc, Only_low_MQ_perc, other_perc)
            f_bar_perc.write(bar_perc_string)
        f_bar_count.close()
        f_bar_perc.close()
        bargraph_R_script = "library(ggplot2)\nlibrary(reshape)\nx1 <- read.table(\"bargraph_percentage.txt\", header=TRUE)\nx1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))\nmdf1=melt(x1,id.vars=\"Sample\")\npdf(\"%s/%s_barplot.pdf\", width = 30, height = 30)\nggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat=\"identity\") + ylab(\"Percentage of Filtered Positions\") + xlab(\"Samples\") + theme(text = element_text(size=9)) + scale_fill_manual(name=\"Reason for filtered out positions\", values=c(\"#08306b\", \"black\", \"orange\", \"darkgrey\", \"#fdd0a2\", \"#7f2704\")) + ggtitle(\"Title Here\") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = \"horizontal\")\ndev.off()" % (args.filter2_only_snp_vcf_dir, os.path.basename(os.path.normpath(args.results_dir)))
        barplot_R_file = open("%s/bargraph.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        keep_logging('Run this R script to generate bargraph plot: %s/bargraph.R' % args.filter2_only_snp_vcf_dir, 'Run this R script to generate bargraph plot: %s/bargraph.R' % args.filter2_only_snp_vcf_dir, logger, 'info')
    """ Methods Steps"""
    keep_logging('Running: Generating data matrices...', 'Running: Generating data matrices...', logger, 'info')
    generate_position_label_data_matrix_All_label()
    keep_logging('Running: Changing variables in data matrices to codes for faster processing...', 'Running: Changing variables in data matrices to codes for faster processing...', logger, 'info')
    temp_generate_position_label_data_matrix_All_label()
    keep_logging('Running: Generating Barplot statistics data matrices...', 'Running: Generating Barplot statistics data matrices...', logger, 'info')
    barplot_stats()

def generate_indel_position_label_data_matrix():

    """
    Generate different list of Positions using the matrix All_label_final_sorted_header.txt.

    (Defining Core Variant Position: Variant Position which was not filtered out in any of the other samples due to variant filter parameter and also this position was present in all the samples(not unmapped)).

    Filtered Position label matrix:
        List of non-core positions. These positions didn't make it to the final core list because it was filtered out in one of the samples.

    Only_ref_variant_positions_for_closely_matrix.txt :
        Those Positions where the variant was either reference allele or a variant that passed all the variant filter parameters.

    :param: null
    :return: null

    """
    def generate_indel_position_label_data_matrix_All_label():
        position_label = OrderedDict()
        print "Generating Only_ref_indel_positions_for_closely"
        f1=open("%s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'w+')
        f2=open("%s/Only_ref_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f3=open("%s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f4=open("%s/Only_filtered_indel_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            keep_logging('Reading All label positions file: %s/All_indel_label_final_sorted_header.txt' % args.filter2_only_snp_vcf_dir, 'Reading All label positions file: %s/All_indel_label_final_sorted_header.txt' % args.filter2_only_snp_vcf_dir, logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                position_label[row[0]] = row[1:]
            keep_logging('Generating different list of Positions and heatmap data matrix...', 'Generating different list of Positions and heatmap data matrix...', logger, 'info')
            print_string_header = "\t"
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"
            #f.write('\t' + print_string_header.strip() + '\n')
            f2.write('\t' + print_string_header.strip() + '\n')
            f3.write('\t' + print_string_header.strip() + '\n')
            f4.write('\t' + print_string_header.strip() + '\n')
            for value in position_label:
                lll = ['0', '2', '3', '4', '5', '6', '7']
                ref_var = ['1', '1TRUE']
                if set(ref_var) & set(position_label[value]):
                    if set(lll) & set(position_label[value]):
                        print_string = ""
                        for i in position_label[value]:
                            print_string = print_string + "\t" + i
                        STRR2 = value + print_string + "\n"
                        f3.write(STRR2)
                        if position_label[value].count('1TRUE') >= 2:
                            f4.write('1\n')
                        else:
                            f4.write('0\n')
                    else:
                        strr = value + "\n"
                        f1.write(strr)
                        STRR3 =	value +	"\t" + str(position_label[value]) + "\n"
                        f2.write(STRR3)
        csv_file.close()
        f1.close()
        f2.close()
        f3.close()
        f4.close()
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_indel_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/1TRUE/-1/g' %s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)

    def temp_generate_indel_position_label_data_matrix_All_label():

        """
        Read **temp_label_final_raw.txt** SNP position label data matrix for generating barplot statistics.
        """
        temp_position_label = OrderedDict()
        f33=open("%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        print_string_header = "\t"
        for i in vcf_filenames:
            print_string_header = print_string_header + os.path.basename(i) + "\t"
        f33.write('\t' + print_string_header.strip() + '\n')
        keep_logging('Reading temporary label positions file: %s/temp_label_final_raw.txt' % args.filter2_only_snp_vcf_dir, 'Reading temporary label positions file: %s/temp_label_final_raw.txt' % args.filter2_only_snp_vcf_dir, logger, 'info')
        lll = ['reference_unmapped_position', 'LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP', 'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ', 'LowFQ_proximate_SNP', 'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP', 'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP', 'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
        ref_var = ['reference_allele', 'VARIANT']
        with open("%s/temp_indel_label_final_raw.txt" % args.filter2_only_snp_vcf_dir, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                if set(ref_var) & set(row[1:]):
                    if set(lll) & set(row[1:]):
                        print_string = ""
                        for i in row[1:]:
                            print_string = print_string + "\t" + i
                        STRR2 = row[0] + print_string + "\n"
                        f33.write(STRR2)
        csv_file.close()
        f33.close()
        """
        Read temp_Only_filtered_positions_for_closely_matrix file and generate a matrix of positions that are being filtered just because of FQ
        """
        temp_position_label_FQ = OrderedDict()
        f44=open("%s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            keep_logging('Reading temporary Only_filtered_indel_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir, 'Reading temporary Only_filtered_indel_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir, logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)

            for row in csv_reader:
              temp_position_label_FQ[row[0]] = row[1:]
            print_string_header = "\t"
            for i in vcf_filenames:
              print_string_header = print_string_header + os.path.basename(i) + "\t"
            f44.write('\t' + print_string_header.strip() + '\n')
            for value in temp_position_label_FQ:
              lll = ['LowFQ']
              if set(lll) & set(temp_position_label_FQ[value]):
                  print_string = ""
                  for i in temp_position_label_FQ[value]:
                      print_string = print_string + "\t" + i
                  STRR2 = value + print_string + "\n"
                  f44.write(STRR2)
            f44.close()
            csv_file.close()
            f44.close()

        """
        Perform Sed on temp files. Find a faster way to do this.
        """
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ/3/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir], shell=True)


        """
        Read temp_Only_filtered_positions_for_closely_matrix file and generate a matrix of positions that are being filtered just because of Dp
        """
        temp_position_label_DP = OrderedDict()
        f44=open("%s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            keep_logging('Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir, 'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir, logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                temp_position_label_DP[row[0]] = row[1:]
            print_string_header = "\t"
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"
            f44.write('\t' + print_string_header.strip() + '\n')
            for value in temp_position_label_DP:
                lll = ['HighFQ_DP']
                ref_var = ['reference_allele', 'VARIANT']
                if set(lll) & set(temp_position_label_FQ[value]):
                    print_string = ""
                    for i in temp_position_label_FQ[value]:
                        print_string = print_string + "\t" + i
                    STRR2 = value + print_string + "\n"
                    f44.write(STRR2)
        f44.close()
        csv_file.close()

        """
        Perform Sed on temp files. Find a faster way to do this.
        """
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ_DP/3/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/LowFQ/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/HighFQ/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir], shell=True)


    def barplot_indel_stats():
        keep_logging('Read each Sample columns and calculate the percentage of each label to generate barplot statistics.', 'Read each Sample columns and calculate the percentage of each label to generate barplot statistics.', logger, 'info')
        """
        Read each Sample columns and calculate the percentage of each label to generate barplot statistics.
        This will give a visual explanation of how many positions in each samples were filtered out because of different reason
        """

        c_reader = csv.reader(open('%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
        columns = list(zip(*c_reader))
        keep_logging('Finished reading columns...', 'Finished reading columns...', logger, 'info')
        counts = 1
        end = len(vcf_filenames) + 1
        f_bar_count = open("%s/bargraph_indel_counts.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_perc = open("%s/bargraph_indel_percentage.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_count.write("Sample\tunmapped_positions\treference_allele\ttrue_variant\tOnly_low_FQ\tOnly_DP\tOnly_low_MQ\tother\n")
        f_bar_perc.write("Sample\tunmapped_positions_perc\ttrue_variant_perc\tOnly_low_FQ_perc\tOnly_DP_perc\tOnly_low_MQ_perc\tother_perc\n")
        for i in xrange(1, end, 1):
            """ Bar Count Statistics: Variant Position Count Statistics """
            true_variant = columns[i].count('VARIANT')
            unmapped_positions = columns[i].count('reference_unmapped_position')
            reference_allele = columns[i].count('reference_allele')
            Only_low_FQ = columns[i].count('LowFQ')
            Only_DP = columns[i].count('HighFQ_DP')
            Only_low_MQ = columns[i].count('HighFQ')
            low_FQ_other_parameters = columns[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns[i].count('LowFQ_DP_QUAL_proximate_SNP') + columns[i].count('LowFQ_QUAL_proximate_SNP') + columns[i].count('LowFQ_DP_proximate_SNP') + columns[i].count('LowFQ_proximate_SNP') + columns[i].count('LowFQ_QUAL_DP') + columns[i].count('LowFQ_DP_QUAL') + columns[i].count('LowFQ_QUAL') + columns[i].count('LowFQ_DP')
            high_FQ_other_parameters = columns[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns[i].count('HighFQ_DP_QUAL_proximate_SNP') + columns[i].count('HighFQ_QUAL_proximate_SNP') + columns[i].count('HighFQ_DP_proximate_SNP') + columns[i].count('HighFQ_proximate_SNP') + columns[i].count('HighFQ_QUAL_DP') + columns[i].count('HighFQ_DP_QUAL') + columns[i].count('HighFQ_QUAL')
            other = low_FQ_other_parameters + high_FQ_other_parameters
            total = true_variant + unmapped_positions + reference_allele + Only_low_FQ + Only_DP + low_FQ_other_parameters + high_FQ_other_parameters + Only_low_MQ
            filename_count = i - 1
            bar_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), unmapped_positions, reference_allele, true_variant, Only_low_FQ, Only_DP, Only_low_MQ, other)
            f_bar_count.write(bar_string)

            """ Bar Count Percentage Statistics: Variant Position Percentage Statistics """
            try:
                true_variant_perc = float((columns[i].count('VARIANT') * 100) / total)
            except ZeroDivisionError:
                true_variant_perc = 0
            try:
                unmapped_positions_perc = float((columns[i].count('reference_unmapped_position') * 100) / total)
            except ZeroDivisionError:
                unmapped_positions_perc = 0
            try:
                reference_allele_perc = float((columns[i].count('reference_allele') * 100) / total)
            except ZeroDivisionError:
                reference_allele_perc = 0
            try:
                Only_low_FQ_perc = float((columns[i].count('LowFQ') * 100) / total)
            except ZeroDivisionError:
                Only_low_FQ_perc = 0
            try:
                Only_DP_perc = float((columns[i].count('HighFQ_DP') * 100) / total)
            except ZeroDivisionError:
                Only_DP_perc = 0
            try:
                Only_low_MQ_perc = float((columns[i].count('HighFQ') * 100) / total)
            except ZeroDivisionError:
                Only_low_MQ_perc = 0
            try:
                low_FQ_other_parameters_perc = float(((columns[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns[i].count('LowFQ_DP_QUAL_proximate_SNP') + columns[i].count('LowFQ_QUAL_proximate_SNP') + columns[i].count('LowFQ_DP_proximate_SNP') + columns[i].count('LowFQ_proximate_SNP') + columns[i].count('LowFQ_QUAL_DP') + columns[i].count('LowFQ_DP_QUAL') + columns[i].count('LowFQ_QUAL') + columns[i].count('LowFQ_DP'))  * 100) / total)
            except ZeroDivisionError:
                low_FQ_other_parameters_perc = 0
            try:
                high_FQ_other_parameters_perc = float(((columns[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns[i].count('HighFQ_DP_QUAL_proximate_SNP') + columns[i].count('HighFQ_QUAL_proximate_SNP') + columns[i].count('HighFQ_DP_proximate_SNP') + columns[i].count('HighFQ_proximate_SNP') + columns[i].count('HighFQ_QUAL_DP') + columns[i].count('HighFQ_DP_QUAL') + columns[i].count('HighFQ_QUAL')) * 100) / total)
            except ZeroDivisionError:
                high_FQ_other_parameters_perc = 0

            other_perc = float(low_FQ_other_parameters_perc + high_FQ_other_parameters_perc)
            bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), unmapped_positions_perc, true_variant_perc, Only_low_FQ_perc, Only_DP_perc, Only_low_MQ_perc, other_perc)
            f_bar_perc.write(bar_perc_string)
        f_bar_count.close()
        f_bar_perc.close()
        bargraph_R_script = "library(ggplot2)\nlibrary(reshape)\nx1 <- read.table(\"bargraph_indel_percentage.txt\", header=TRUE)\nx1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))\nmdf1=melt(x1,id.vars=\"Sample\")\npdf(\"%s/%s_barplot_indel.pdf\", width = 30, height = 30)\nggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat=\"identity\") + ylab(\"Percentage of Filtered Positions\") + xlab(\"Samples\") + theme(text = element_text(size=9)) + scale_fill_manual(name=\"Reason for filtered out positions\", values=c(\"#08306b\", \"black\", \"orange\", \"darkgrey\", \"#fdd0a2\", \"#7f2704\")) + ggtitle(\"Title Here\") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = \"horizontal\")\ndev.off()"  % (args.filter2_only_snp_vcf_dir, os.path.basename(os.path.normpath(args.results_dir)))
        barplot_R_file = open("%s/bargraph_indel.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        keep_logging('Run this R script to generate bargraph plot: %s/bargraph_indel.R' % args.filter2_only_snp_vcf_dir, 'Run this R script to generate bargraph plot: %s/bargraph_indel.R' % args.filter2_only_snp_vcf_dir, logger, 'info')


    """ Methods Steps"""
    keep_logging('Running: Generating data matrices...', 'Running: Generating data matrices...', logger, 'info')
    generate_indel_position_label_data_matrix_All_label()
    keep_logging('Running: Changing variables in data matrices to codes for faster processing...', 'Running: Changing variables in data matrices to codes for faster processing...', logger, 'info')
    temp_generate_indel_position_label_data_matrix_All_label()
    keep_logging('Running: Generating Barplot statistics data matrices...', 'Running: Generating Barplot statistics data matrices...', logger, 'info')
    barplot_indel_stats()

def create_job_fasta(jobrun, vcf_filenames, core_vcf_fasta_dir, functional_filter):

    """ Generate jobs/scripts that creates core consensus fasta file.

    This function will generate and run scripts/jobs to create core consensus fasta file of only core variant positions.
    Input for Fasttree, Beast and pairwise variant analysis.

    :param jobrun: Based on this value all the job/scripts will run on "cluster": either on single cluster, "parallel-local": run in parallel on local system, "local": run on local system, "parallel-cluster": submit parallel jobs on cluster.
    :param vcf_filenames: list of final vcf filenames i.e *_no_proximate_snp.vcf. These files are the final output of variant calling step for each sample.
    :return:
    :raises:
    """
    if jobrun == "parallel-cluster":
        """
        Supports only PBS clusters for now.
        """
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            #os.system("qsub %s" % i)
            call("qsub %s" % i, logger)


    elif jobrun == "parallel-local" or jobrun == "cluster":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        if args.numcores:
            num_cores = int(num_cores)
        else:
            num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    # elif jobrun == "cluster":
    #     command_array = []
    #     command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
    #     f3 = open(command_file, 'w+')
    #     for i in vcf_filenames:
    #         job_name = os.path.basename(i)
    #         job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
    #         job_file_name = "%s_fasta.pbs" % (i)
    #         f1=open(job_file_name, 'w+')
    #         f1.write(job_print_string)
    #         f1.close()
    #     pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
    #     pbs_scripts = glob.glob(pbs_dir)
    #     for i in pbs_scripts:
    #         f3.write("bash %s\n" % i)
    #     f3.close()
    #     with open(command_file, 'r') as fpp:
    #         for lines in fpp:
    #             lines = lines.strip()
    #             command_array.append(lines)
    #     fpp.close()
    #     os.system("bash %s/command_file" % args.filter2_only_snp_vcf_dir)
    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_array = []
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)


        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        #os.system("bash command_file")
        call("bash %s" % command_file, logger)

def create_job_allele_variant_fasta(jobrun, vcf_filenames, core_vcf_fasta_dir, config_file):

    """ Generate jobs/scripts that creates core consensus fasta file.

    This function will generate and run scripts/jobs to create core consensus fasta file of only core variant positions.
    Input for Fasttree, Beast and pairwise variant analysis.

    :param jobrun: Based on this value all the job/scripts will run on "cluster": either on single cluster, "parallel-local": run in parallel on local system, "local": run on local system, "parallel-cluster": submit parallel jobs on cluster.
    :param vcf_filenames: list of final vcf filenames i.e *_no_proximate_snp.vcf. These files are the final output of variant calling step for each sample.
    :return:
    :raises:
    """
    if jobrun == "parallel-cluster":
        """
        Supports only PBS clusters for now.
        """
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            #os.system("qsub %s" % i)
            call("qsub %s" % i, logger)


    elif jobrun == "parallel-local" or jobrun == "cluster":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list_ref_allele_variants_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_ref_allele_variants_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        if args.numcores:
            num_cores = int(num_cores)
        else:
            num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    # elif jobrun == "cluster":
    #     command_array = []
    #     command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
    #     f3 = open(command_file, 'w+')
    #     for i in vcf_filenames:
    #         job_name = os.path.basename(i)
    #         job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
    #         job_file_name = "%s_fasta.pbs" % (i)
    #         f1=open(job_file_name, 'w+')
    #         f1.write(job_print_string)
    #         f1.close()
    #     pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
    #     pbs_scripts = glob.glob(pbs_dir)
    #     for i in pbs_scripts:
    #         f3.write("bash %s\n" % i)
    #     f3.close()
    #     with open(command_file, 'r') as fpp:
    #         for lines in fpp:
    #             lines = lines.strip()
    #             command_array.append(lines)
    #     fpp.close()
    #     os.system("bash %s/command_file" % args.filter2_only_snp_vcf_dir)
    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_array = []
        command_file = "%s/commands_list_ref_allele_variants_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_ref_allele_variants_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)


        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        #os.system("bash command_file")
        call("bash %s" % command_file, logger)

def create_job_DP(jobrun, vcf_filenames):
    """
    Based on type of jobrun; generate jobs and run accordingly.
    :param jobrun: Based on this value all the job/scripts will run on "cluster": either on single cluster, "parallel-local": run in parallel on local system, "local": run on local system, "parallel-cluster": submit parallel jobs on cluster.
    :param vcf_filenames:
    :return:
    """

    if jobrun == "parallel-cluster":
        """
        Supports only PBS clusters for now.
        """
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            #os.system("qsub %s" % i)
            call("qsub %s" % i, logger)


    elif jobrun == "parallel-local" or jobrun == "cluster" :
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
        pbs_scripts = glob.glob(pbs_dir)


        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        with open(command_file, 'r') as fpp:
            for lines in fpp:
                lines = lines.strip()
                command_array.append(lines)
        fpp.close()
        print len(command_array)
        if args.numcores:
            num_cores = int(num_cores)
        else:
            num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    # elif jobrun == "cluster":
    #     """ Test pending """
    #     command_file = "%s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir
    #     f3 = open(command_file, 'w+')
    #     for i in vcf_filenames:
    #         job_name = os.path.basename(i)
    #         job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
    #         job_file_name = "%s_DP.pbs" % (i)
    #         f1=open(job_file_name, 'w+')
    #         f1.write(job_print_string)
    #         f1.close()
    #     pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
    #     pbs_scripts = glob.glob(pbs_dir)
    #     for i in pbs_scripts:
    #         f3.write("bash %s\n" % i)
    #     f3.close()
    #     os.system("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir)

    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = "%s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        #os.system("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir)
        call("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir, logger)

def generate_vcf_files():
    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
        keep_logging('Removing Variants falling in Functional filters positions file: %s\n' % functional_class_filter_positions, 'Removing Variants falling in Functional filters positions file: %s\n' % functional_class_filter_positions, logger,
                     'info')
        # phage_positions = []
        # phage_region_positions = "%s/phage_region_positions.txt" % args.filter2_only_snp_vcf_dir
        # with open(phage_region_positions, 'rU') as fp:
        #     for line in fp:
        #         phage_positions.append(line.strip())
        # fp.close()


        functional_filter_pos_array = []
        with open(functional_class_filter_positions, 'rU') as f_functional:
            for line_func in f_functional:
                functional_filter_pos_array.append(line_func.strip())

        ref_variant_position_array = []
        ffp = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'r+')
        for line in ffp:
            line = line.strip()
            if line not in functional_filter_pos_array:
                ref_variant_position_array.append(line)
        ffp.close()

        # Adding core indel support: 2018-07-24
        ref_indel_variant_position_array = []
        ffp = open("%s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'r+')
        for line in ffp:
            line = line.strip()
            if line not in functional_filter_pos_array:
                ref_indel_variant_position_array.append(line)
        ffp.close()

    else:
        functional_filter_pos_array = []
        ref_variant_position_array = []
        ffp = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'r+')
        for line in ffp:
            line = line.strip()
            ref_variant_position_array.append(line)
        ffp.close()

        # Adding core indel support: 2018-07-24
        ref_indel_variant_position_array = []
        ffp = open("%s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'r+')
        for line in ffp:
            line = line.strip()
            if line not in functional_filter_pos_array:
                ref_indel_variant_position_array.append(line)
        ffp.close()

    print "No. of core SNPs: %s" % len(ref_variant_position_array)
    print "No. of core INDELs: %s" % len(ref_indel_variant_position_array)

    f_file = open("%s/Only_ref_variant_positions_for_closely_without_functional_filtered_positions" % args.filter2_only_snp_vcf_dir, 'w+')
    for pos in ref_variant_position_array:
        f_file.write(pos + '\n')
    f_file.close()

    # Adding core indel support: 2018-07-24
    f_file = open(
        "%s/Only_ref_indel_variant_positions_for_closely_without_functional_filtered_positions" % args.filter2_only_snp_vcf_dir,
        'w+')
    for pos in ref_indel_variant_position_array:
        f_file.write(pos + '\n')
    f_file.close()

    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    filter2_files_array = []
    for i in vcf_filenames:
        filter2_file = i.replace('_no_proximate_snp.vcf', '')
        filter2_files_array.append(filter2_file)


    filtered_out_vcf_files = []
    for i in filter2_files_array:
        print_array =[]
        with open(i) as file_open:
            for line in file_open:
                line = line.strip()
                if line.startswith("#"):
                    print_array.append(line)
                else:
                    split_array = re.split(r'\t+', line)
                    if split_array[1] in ref_variant_position_array and 'INDEL' not in split_array[7]:
                        print_array.append(line)
        file_open.close()
        file_name = i + "_core.vcf"
        keep_logging('Generating %s' % file_name, 'Generating %s' % file_name, logger, 'info')
        filtered_out_vcf_files.append(file_name)
        f1 = open(file_name, 'w+')
        for ios in print_array:
            print_string = str(ios) + "\n"
            f1.write(print_string)
        f1.close()

    filename = "%s/consensus.sh" % args.filter2_only_snp_vcf_dir
    keep_logging('Generating Consensus...', 'Generating Consensus...', logger, 'info')
    for file in filtered_out_vcf_files:
        f1 = open(filename, 'a+')
        bgzip_cmd = "%s/%s/bgzip -f %s\n" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], file)
        f1.write(bgzip_cmd)
        subprocess.call([bgzip_cmd], shell=True)
        tabix_cmd = "%s/%s/tabix -f -p vcf %s.gz\n" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], file)
        f1.write(tabix_cmd)
        subprocess.call([tabix_cmd], shell=True)
        fasta_cmd = "cat %s | %s/vcf-consensus %s.gz > %s.fa\n" % (args.reference, base_vcftools_bin, file, file.replace('_filter2_final.vcf_core.vcf', ''))
        f1.write(fasta_cmd)
        subprocess.call([fasta_cmd], shell=True)
        base = os.path.basename(file)
        header = base.replace('_filter2_final.vcf_core.vcf', '')
        sed_command = "sed -i 's/>.*/>%s/g' %s.fa\n" % (header, file.replace('_filter2_final.vcf_core.vcf', ''))
        subprocess.call([sed_command], shell=True)
        f1.write(sed_command)
    keep_logging('The consensus commands are in : %s' % filename, 'The consensus commands are in : %s' % filename, logger, 'info')
    sequence_lgth_cmd = "for i in %s/*.fa; do %s/%s/bioawk -c fastx \'{ print $name, length($seq) }\' < $i; done" % (args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bioawk", Config)['bioawk_bin'])
    #os.system(sequence_lgth_cmd)
    call("%s" % sequence_lgth_cmd, logger)

def gatk_filter2(final_raw_vcf, out_path, analysis, reference):
    gatk_filter2_parameter_expression = "MQ > 50 && QUAL > 100 && DP > 9"
    gatk_filter2_command = "java -jar %s/%s/GenomeAnalysisTK.jar -T VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s --filterExpression \"%s\" --filterName PASS_filter2" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("gatk", Config)['gatk_bin'], reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
    keep_logging('Running Command: [%s]' % gatk_filter2_command, 'Running Command: [%s]' % gatk_filter2_command, logger, 'info')
    #os.system(gatk_filter2_command)
    call("%s" % gatk_filter2_command, logger)
    filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_gatk.vcf > %s/%s_filter2_final.vcf" % (out_path, analysis, out_path, analysis)
    call("%s" % filter_flag_command, logger)
    gatk_filter2_final_vcf = "%s/%s_filter2_final.vcf" % (out_path, analysis)
    return gatk_filter2_final_vcf

def remove_proximate_snps(gatk_filter2_final_vcf_file, out_path, analysis, reference):
    all_position = []
    remove_proximate_position_array = []
    gatk_filter2_final_vcf_file_no_proximate_snp = gatk_filter2_final_vcf_file + "_no_proximate_snp.vcf"
    with open(gatk_filter2_final_vcf_file, 'rU') as csv_file:
        for line in csv_file:
            if not line.startswith('#'):
                line_array = line.split('\t')
                all_position.append(line_array[1])
    for position in all_position:
        position_index = all_position.index(position)
        next_position_index = position_index + 1

        if next_position_index < len(all_position):
            diff = int(all_position[next_position_index]) - int(position)
            if diff < 10:
                #print position + "  " + all_position[next_position_index]
                if position not in remove_proximate_position_array and all_position[next_position_index] not in remove_proximate_position_array:
                    remove_proximate_position_array.append(int(position))
                    remove_proximate_position_array.append(int(all_position[next_position_index]))
    f1=open(gatk_filter2_final_vcf_file_no_proximate_snp, 'w+')
    with open(gatk_filter2_final_vcf_file, 'rU') as csv_file2:
        for line in csv_file2:
            if line.startswith('gi') or line.startswith('MRSA_8058'): ##change this!
               line_array = line.split('\t')
               if int(line_array[1]) not in remove_proximate_position_array:
                   print_string = line
                   f1.write(print_string)
            else:
                print_string = line
                f1.write(print_string)
    gatk_filter2_final_vcf_file_no_proximate_snp_positions = gatk_filter2_final_vcf_file + "_no_proximate_snp.vcf_positions_array"
    f2=open(gatk_filter2_final_vcf_file_no_proximate_snp_positions, 'w+')
    for i in remove_proximate_position_array:
        position_print_string = str(i) + "\n"
        f2.write(position_print_string)
    return gatk_filter2_final_vcf_file_no_proximate_snp

def FQ_analysis():
    for i in vcf_filenames:
        filename_base = os.path.basename(i)
        aln_mpileup_vcf_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf_5bp_indel_removed.vcf')
        analysis = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')
        #print aln_mpileup_vcf_file
        grep_reference_file = "grep \'^##reference\' %s" % aln_mpileup_vcf_file
        proc = subprocess.Popen([grep_reference_file], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        reference_file = out.split(':')
        gatk_filter2_final_vcf_file = gatk_filter2(aln_mpileup_vcf_file, temp_dir, analysis, reference_file[1])
        #print gatk_filter2_final_vcf_file
        gatk_filter2_final_vcf_file_no_proximate_snp = remove_proximate_snps(gatk_filter2_final_vcf_file, temp_dir, analysis, reference_file[1])
        grep_fq_field = "awk -F\'\\t\' \'{print $8}\' %s | grep -o \'FQ=.*\' | sed \'s/FQ=//g\' | awk -F\';\' \'{print $1}\' > %s/%s_FQ_values" % (gatk_filter2_final_vcf_file_no_proximate_snp, os.path.dirname(i), analysis)
        #os.system(grep_fq_field)
        call("%s" % grep_fq_field, logger)
        #print grep_fq_field

def DP_analysis():
    create_job_DP(args.jobrun, vcf_filenames)
    paste_command = "paste %s/extract_DP_positions.txt" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_DP_values')
        paste_command = paste_command + " " + label_file

    paste_file = args.filter2_only_snp_vcf_dir + "/paste_DP_files.sh"
    f2=open(paste_file, 'w+')
    paste_command = paste_command + " > %s/filtered_DP_values_temp.txt" % args.filter2_only_snp_vcf_dir
    #os.system(paste_command)
    f2.write(paste_command + '\n')
    cat_header = "cat %s/header.txt %s/filtered_DP_values_temp.txt > %s/filtered_DP_values.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    #os.system(cat_header)
    f2.write(cat_header + '\n')
    sed_command = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/filtered_DP_values.txt" % (args.filter2_only_snp_vcf_dir)
    #os.system(sed_command)
    f2.write(sed_command + '\n')
    cmd = "bash %s" % paste_file
    # os.system("bash %s/paste_DP_files.sh" % args.filter2_only_snp_vcf_dir)

def DP_analysis_barplot():
    #os.system("bash %s/paste_DP_files.sh" % args.filter2_only_snp_vcf_dir)
    call("bash %s/paste_DP_files.sh" % args.filter2_only_snp_vcf_dir, logger)
    keep_logging('Generating DP barplots data...', 'Generating DP barplots data...', logger, 'info')
    c_reader = csv.reader(open('%s/filtered_DP_values.txt' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
    columns = list(zip(*c_reader))
    counts = 1
    end = len(vcf_filenames) + 1
    f_bar_count = open("%s/DP_bargraph_counts.txt" % args.filter2_only_snp_vcf_dir, 'w+')
    f_bar_perc = open("%s/DP_bargraph_percentage.txt" % args.filter2_only_snp_vcf_dir, 'w+')
    f_bar_count.write("Sample\treference_position\toneto5\tsixto10\televento14\tfifteenorabove\n")
    f_bar_perc.write("Sample\treference_position\toneto5\tsixto10\televento14\tfifteenorabove\n")
    for i in xrange(1, end, 1):
        """ Bar Count Statistics: Variant Position Count Statistics """
        reference_position = columns[i].count('NA')
        oneto5 = 0
        for k in list(columns[i][1:]):
            if k != "":
                if k != "NA":
                    if int(k) < 5:
                        oneto5 += 1
        sixto10 = 0
        for k in list(columns[i][1:]):
            if k != "":
                if k != "NA":
                    if int(k) >= 5 and int(k) <= 10:
                        sixto10 += 1
        elevento14 = 0
        for k in list(columns[i][1:]):
            if k != "":
                if k != "NA":
                    if int(k) >= 11 and int(k) <= 14:
                        elevento14 += 1
        fifteenorabove = 0
        for k in list(columns[i][1:]):
            if k != "":
                if k != "NA":
                    if int(k) >= 15:
                        fifteenorabove += 1
        total = reference_position + oneto5 + sixto10 + elevento14 + fifteenorabove
        filename_count = i - 1
        bar_string = "%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), reference_position, oneto5, sixto10, elevento14, fifteenorabove)
        f_bar_count.write(bar_string)

        """ Bar Count Percentage Statistics: Variant Position Percentage Statistics """
        try:
            reference_position_perc = float(reference_position * 100 / total)
        except ZeroDivisionError:
            reference_position_perc = 0
        try:
            oneto5_perc = float(oneto5 * 100 / total)
        except ZeroDivisionError:
            oneto5_perc = 0
        try:
            sixto10_perc = float(sixto10 * 100 / total)
        except ZeroDivisionError:
            sixto10_perc = 0
        try:
            elevento14_perc = float(elevento14 * 100 / total)
        except ZeroDivisionError:
            elevento14_perc = 0
        try:
            fifteenorabove_perc = float(fifteenorabove * 100 / total)
        except ZeroDivisionError:
            fifteenorabove_perc = 0
        bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), reference_position_perc, oneto5_perc, sixto10_perc, elevento14_perc, fifteenorabove_perc)
        f_bar_perc.write(bar_perc_string)

def extract_only_ref_variant_fasta(core_vcf_fasta_dir):
    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes" and ConfigSectionMap("functional_filters", Config)['apply_to_calls'] == "yes":
        functional_filter = "yes"
    create_job_fasta(args.jobrun, vcf_filenames, core_vcf_fasta_dir, functional_filter)

def extract_only_ref_variant_fasta_from_reference():
    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes" and \
            ConfigSectionMap("functional_filters", Config)['apply_to_calls'] == "yes":
        ffp = open("%s/Only_ref_variant_positions_for_closely_without_functional_filtered_positions" % args.filter2_only_snp_vcf_dir).readlines()
    else:
        ffp = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir).readlines()
    fasta_string = ""
    #firstLine = ffp.pop(0)
    for lines in ffp:
        lines = lines.strip()
        extract_base = "grep -v \'>\' %s | tr -d \'\\n\'| cut -b%s" % (args.reference, lines)
        proc = subprocess.Popen([extract_base], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        fasta_string = fasta_string + out
        if not out:
            print lines
            keep_logging('Error extracting reference allele', 'Error extracting reference allele', logger, 'info')
            exit()

    pattern = re.compile(r'\s+')
    fasta_string = re.sub(pattern, '', fasta_string)
    final_fasta_string = ">%s\n" % os.path.basename(args.reference.replace('.fasta', '').replace('.fa', '')) + fasta_string + "\n"
    fp = open("%s/%s_variants.fa" % (args.filter2_only_snp_vcf_dir, os.path.basename(args.reference.replace('.fasta', '').replace('.fa', ''))), 'w+')
    fp.write(final_fasta_string)
    fp.close()


# This got affected by unsorted unique_positions_file
def extract_only_ref_variant_fasta_from_reference_allele_variant():
    ffp = open("%s/unique_positions_file" % args.filter2_only_snp_vcf_dir).readlines()
    #unique_positions_array = []

    fasta_string = ""
    #firstLine = ffp.pop(0)
    for lines in ffp:
        lines = lines.strip()
        #unique_positions_array.append(lines)
        extract_base = "grep -v \'>\' %s | tr -d \'\\n\'| cut -b%s" % (args.reference, lines)
        proc = subprocess.Popen([extract_base], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        fasta_string = fasta_string + out
        if not out:
            print lines
            keep_logging('Error extracting reference allele', 'Error extracting reference allele', logger, 'info')
            exit()

    pattern = re.compile(r'\s+')
    fasta_string = re.sub(pattern, '', fasta_string)
    final_fasta_string = ">%s\n" % os.path.basename(args.reference.replace('.fasta', '').replace('.fa', '')) + fasta_string + "\n"
    fp = open("%s/%s_allele_variants.fa" % (args.filter2_only_snp_vcf_dir, os.path.basename(args.reference.replace('.fasta', '').replace('.fa', ''))), 'w+')
    fp.write(final_fasta_string)
    fp.close()

def prepare_snpEff_db(reference_basename):
    keep_logging('Preparing snpEff database requirements.', 'Preparing snpEff database requirements.', logger, 'info')
    reference_basename = (os.path.basename(args.reference)).split(".")
    if os.path.isfile("%s/%s/snpEff.config" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'])):
        #os.system("cp %s/%s/snpEff.config %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], args.filter2_only_snp_vcf_dir))
        keep_logging("cp %s/%s/snpEff.config %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], args.filter2_only_snp_vcf_dir), "cp %s/%s/snpEff.config %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], args.filter2_only_snp_vcf_dir), logger, 'debug')
        call("cp %s/%s/snpEff.config %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], args.filter2_only_snp_vcf_dir), logger)
    else:
        keep_logging("Error: %s/%s/snpEff.config doesn't exists.\nExiting..." % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']),"Error: %s/%s/snpEff.config doesn't exists.\nExiting..." % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger, 'exception')
        exit()
    make_sure_path_exists("%s/%s/data/%s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]))
    make_sure_path_exists("%s/%s/data/genomes/" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']))
    #os.system("cp %s %s/%s/data/genomes/" % (args.reference, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']))
    keep_logging("cp %s %s/%s/data/genomes/%s.fa" % (args.reference, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]), "cp %s %s/%s/data/genomes/" % (args.reference, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger, 'debug')
    call("cp %s %s/%s/data/genomes/%s.fa" % (args.reference, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]), logger)
    with open("%s/snpEff.config" % args.filter2_only_snp_vcf_dir, "a") as conf_file:
        conf_file.write("\n\n##Building Custom Database###\n%s.genome\t: %s\n\n" % (reference_basename[0], reference_basename[0]))
    conf_file.close()
    #get the gff name from config file
    if os.path.isfile("%s/%s.gff" % (os.path.dirname(args.reference), reference_basename[0])):
        keep_logging("cp %s/%s.gff %s/%s/data/%s/genes.gff" % (
        os.path.dirname(args.reference), reference_basename[0], ConfigSectionMap("bin_path", Config)['binbase'],
        ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]),
                     "cp %s/%s.gff %s/%s/data/%s/genes.gff" % (os.path.dirname(args.reference), reference_basename[0],
                                                               ConfigSectionMap("bin_path", Config)['binbase'],
                                                               ConfigSectionMap("snpeff", Config)['snpeff_bin'],
                                                               reference_basename[0]), logger, 'debug')
        keep_logging("cp %s/%s.gb* %s/%s/data/%s/genes.gbk" % (
        os.path.dirname(args.reference), reference_basename[0], ConfigSectionMap("bin_path", Config)['binbase'],
        ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]),
                     "cp %s/%s.gff %s/%s/data/%s/genes.gff" % (os.path.dirname(args.reference), reference_basename[0],
                                                               ConfigSectionMap("bin_path", Config)['binbase'],
                                                               ConfigSectionMap("snpeff", Config)['snpeff_bin'],
                                                               reference_basename[0]), logger, 'debug')
        call("cp %s/%s.gff %s/%s/data/%s/genes.gff" % (
        os.path.dirname(args.reference), reference_basename[0], ConfigSectionMap("bin_path", Config)['binbase'],
        ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]), logger)
        call("cp %s/%s.gb* %s/%s/data/%s/genes.gbk" % (
        os.path.dirname(args.reference), reference_basename[0], ConfigSectionMap("bin_path", Config)['binbase'],
        ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]), logger)
    else:
        keep_logging("Error: %s/%s.gff file doesn't exists. Make sure the GFF file has the same prefix as reference fasta file\nExiting..." % (os.path.dirname(args.reference), reference_basename[0]),
                         "Error: %s/%s.gff file doesn't exists. Make sure the GFF file has the same prefix as reference fasta file\nExiting..." % (os.path.dirname(args.reference), reference_basename[0]), logger, 'exception')
        exit()
    #keep_logging("java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), "java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger, 'debug')
    keep_logging("java -jar %s/%s/%s build -genbank -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), "java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger, 'debug')

    #call("java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger)
    call("java -jar %s/%s/%s build -genbank -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger)
    keep_logging('Finished Preparing snpEff database requirements.', 'Finished Preparing snpEff database requirements.', logger, 'info')

def variant_annotation():
    keep_logging('Annotating Variants using snpEff.', 'Annotating Variants using snpEff.', logger, 'info')

    if ConfigSectionMap("snpeff", Config)['prebuild'] == "yes":
        if ConfigSectionMap("snpeff", Config)['db']:
            print "Using pre-built snpEff database: %s" % ConfigSectionMap("snpeff", Config)['db']
            proc = subprocess.Popen(["java -jar %s/%s/%s databases | grep %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], ConfigSectionMap("snpeff", Config)['db'])],
                                    stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            if out2:
                snpeffdb = ConfigSectionMap("snpeff", Config)['db']
            else:
                print "The database name %s provided was not found. Check the name and try again" % ConfigSectionMap("snpeff", Config)['db']
                exit()
        else:
            print "snpEff db section is not set in config file"
            exit()
    else:
        reference_basename = (os.path.basename(args.reference)).split(".")
        snpeffdb = reference_basename[0]
        prepare_snpEff_db(reference_basename)

    annotate_vcf_cmd_array = []
    annotate_final_vcf_cmd_array = []
    for i in vcf_filenames:
        raw_vcf = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf')
        annotate_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], raw_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, snpeffdb, raw_vcf, raw_vcf)
        #print annotate_vcf_cmd
        annotate_vcf_cmd_array.append(annotate_vcf_cmd)
        final_vcf = i
        annotate_final_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, snpeffdb, final_vcf, final_vcf)
        annotate_final_vcf_cmd_array.append(annotate_final_vcf_cmd)
    if args.numcores:
        num_cores = int(num_cores)
    else:
        num_cores = multiprocessing.cpu_count()
    #print annotate_vcf_cmd_array
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_vcf_cmd_array)
    results_2 = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_final_vcf_cmd_array)

def indel_annotation():
    keep_logging('Annotating indels using snpEff.', 'Annotating indels using snpEff.', logger, 'info')

    if ConfigSectionMap("snpeff", Config)['prebuild'] == "yes":
        if ConfigSectionMap("snpeff", Config)['db']:
            print "Using pre-built snpEff database: %s" % ConfigSectionMap("snpeff", Config)['db']
            proc = subprocess.Popen(["java -jar %s/%s/%s databases | grep %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], ConfigSectionMap("snpeff", Config)['db'])],
                                    stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            if out2:
                snpeffdb = ConfigSectionMap("snpeff", Config)['db']
            else:
                print "The database name %s provided was not found. Check the name and try again" % ConfigSectionMap("snpeff", Config)['db']
                exit()
        else:
            print "snpEff db section is not set in config file"
            exit()
    else:
        reference_basename = (os.path.basename(args.reference)).split(".")
        snpeffdb = reference_basename[0]
        prepare_snpEff_db(reference_basename)




    annotate_vcf_cmd_array = []
    annotate_final_vcf_cmd_array = []
    for i in vcf_filenames:
        raw_vcf = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf')
        annotate_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], raw_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, snpeffdb, raw_vcf, raw_vcf)
        annotate_vcf_cmd_array.append(annotate_vcf_cmd)
        final_vcf = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf')
        annotate_final_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, snpeffdb, final_vcf, final_vcf)
        annotate_final_vcf_cmd_array.append(annotate_final_vcf_cmd)
    if args.numcores:
        num_cores = int(num_cores)
    else:
        num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_vcf_cmd_array)
    results_2 = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_final_vcf_cmd_array)

def gatk_combine_variants(files_gatk, reference, out_path, merged_file_suffix, logger, Config):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)[
        'gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    #files_gatk = "--variant " + ' --variant '.join(vcf_files_array)
    keep_logging("java -jar %s -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s" % (base_cmd, reference, files_gatk, out_path, merged_file_suffix), "java -jar %s -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s" % (base_cmd, reference, files_gatk, out_path, merged_file_suffix), logger, 'debug')
    merge_gatk_commands_file = "%s/gatk_merge.sh" % args.filter2_only_snp_vcf_dir
    with open(merge_gatk_commands_file, 'a+') as fopen:
        fopen.write("java -jar %s -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s" % (base_cmd, reference, files_gatk, out_path, merged_file_suffix) + '\n')
    fopen.close()
    os.system("bash %s" % merge_gatk_commands_file)
    # Commenting out calling gatk combine variants, problem with python subprocess, OSError: [Errno 7] Argument list too long
    #call("java -jar %s -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s" % (base_cmd, reference, files_gatk, out_path, merged_file_suffix), logger)
    #keep_logging("java -jar %s -T CombineVariants -R %s --variant %s -o %s/Final_vcf_gatk%s" % (base_cmd, reference, files_gatk, out_path, merged_file_suffix), "java -jar %s -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s" % (base_cmd, reference, files_gatk, out_path, merged_file_suffix), logger, 'debug')
    #call("java -jar %s -T CombineVariants -R %s --variant %s -o %s/Final_vcf_gatk%s" % (base_cmd, reference, files_gatk, out_path, merged_file_suffix), logger)
    return "%s/Final_vcf_gatk%s" % (out_path, merged_file_suffix)

def annotated_snp_matrix():
    """
    :return: Annotate core vcf files generated at core_prep steps.
    Read Genbank file and return a dictionary of Prokka ID mapped to Gene Name, Prokka ID mapped to Product Name.
    This dictionary will then be used to insert annotation into SNP/Indel matrix
    """

    # # """ Run snpEff annotation step """
    variant_annotation()
    #
    # # """ Run snpEff annotation step """
    indel_annotation()

    # Check if Reference genome Genbank file exists. Read the locus tag and gene annotations into a dictionary that maps locus tags to gene name/product name
    reference_basename = (os.path.basename(args.reference)).split(".")
    if os.path.isfile("%s/%s.gbf" % (os.path.dirname(args.reference), reference_basename[0])):
        handle = open("%s/%s.gbf" % (os.path.dirname(args.reference), reference_basename[0]), 'rU')
    else:
        raise IOError('%s/%s.gbf does not exist.' % (os.path.dirname(args.reference), reference_basename[0]))
        exit()
    locus_tag_to_gene_name = {}
    locus_tag_to_product = {}
    #locus_tag_to_uniprot = {}
    #locus_tag_to_ec_number = {}

    keep_logging(
        'Reading annotations from Reference genome genbank file: %s/%s.gbf' % (os.path.dirname(args.reference), reference_basename[0]),
        'Reading annotations from Reference genome genbank file: %s/%s.gbf' % (os.path.dirname(args.reference), reference_basename[0]),
        logger, 'info')
    for record in SeqIO.parse(handle, 'genbank') :
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers:
                if 'gene' in feature.qualifiers:
                    locus_tag_to_gene_name[str(feature.qualifiers['locus_tag'][0])] = str(feature.qualifiers['gene'][0])
                else:
                    locus_tag_to_gene_name[str(feature.qualifiers['locus_tag'][0])] = "null or hypothetical protein"
                if 'product' in feature.qualifiers:
                    locus_tag_to_product[str(feature.qualifiers['locus_tag'][0])] = str(feature.qualifiers['product'][0])
                else:
                    locus_tag_to_product[str(feature.qualifiers['locus_tag'][0])] = "null or hypothetical protein"
                # elif 'uniprot' in feature.qualifiers:
                #     locus_tag_to_product[str(feature.qualifiers['locus_tag'][0])] = str(feature.qualifiers['product'][0])
            else:
                keep_logging(
                    'Error: locus_tag specifications for the below feature doesnt exists. Please check the format of genbank file\n%s' % str(feature),
                    'Error: locus_tag specifications for the below feature doesnt exists. Please check the format of genbank file\n%s' % str(feature),
                    logger, 'exception')


    """ Merge Annotated final vcf file """
    keep_logging('Merging Final Annotated VCF files into %s/Final_vcf_no_proximate_snp.vcf using bcftools' % args.filter2_only_snp_vcf_dir, 'Merging Final Annotated VCF files into %s/Final_vcf_no_proximate_snp.vcf using bcftools' % args.filter2_only_snp_vcf_dir, logger, 'info')

    files_for_tabix = glob.glob("%s/*.vcf_no_proximate_snp.vcf_ANN.vcf" % args.filter2_only_snp_vcf_dir)
    tabix(files_for_tabix, "vcf", logger, Config)
    files_for_tabix = glob.glob("%s/*_filter2_indel_final.vcf_ANN.vcf" % args.filter2_only_snp_vcf_dir)
    tabix(files_for_tabix, "vcf", logger, Config)

    files = ' '.join(vcf_filenames)


    #Merge with bcftools, ** Deprecated **
    merge_commands_file = "%s/bcftools_merge.sh" % args.filter2_only_snp_vcf_dir
    #"_filter2_indel_final.vcf_ANN.vcf.gz")) + '\n')
    with open(merge_commands_file, 'w+') as fopen:
        fopen.write("%s/%s/bcftools merge -i ANN:join -m both -o %s/Final_vcf_no_proximate_snp.vcf -O v %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bcftools", Config)['bcftools_bin'], args.filter2_only_snp_vcf_dir, files.replace("_filter2_final.vcf_no_proximate_snp.vcf", "_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz")) + '\n')
        fopen.write("%s/%s/bcftools merge -i ANN:join -m both -o %s/Final_vcf_indel.vcf -O v %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bcftools", Config)['bcftools_bin'], args.filter2_only_snp_vcf_dir,files.replace("_filter2_final.vcf_no_proximate_snp.vcf","_filter2_indel_final.vcf_ANN.vcf.gz")) + '\n')


    fopen.close()

    os.system("bash %s" % merge_commands_file)




    # Merge with Gatk combine variants method
    merged_file_suffix = "_no_proximate_snp.vcf"

    annotated_no_proximate_snp_file = "%s/annotated_no_proximate_snp_list.txt" % args.filter2_only_snp_vcf_dir
    annotated_no_proximate_snp_indel_file = "%s/annotated_no_proximate_snp_indel_list.txt" % args.filter2_only_snp_vcf_dir

    with open(annotated_no_proximate_snp_file, 'w+') as fopen:
        for i in vcf_filenames:
            fopen.write(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz') + '\n')
    fopen.close()

    with open(annotated_no_proximate_snp_indel_file, 'w+') as fopen:
        for i in vcf_filenames:
            fopen.write(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf_ANN.vcf.gz') + '\n')
    fopen.close()

    #files_gatk = "--variant " + ' --variant '.join(vcf_filenames)
    files_gatk = ""
    for i in vcf_filenames:
        files_gatk = files_gatk + " --variant " + i
    final_gatk_snp_merged_vcf = gatk_combine_variants(files_gatk.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz'), args.reference, args.filter2_only_snp_vcf_dir, merged_file_suffix, logger, Config)
    merged_file_suffix = "_indel.vcf"
    final_gatk_indel_merged_vcf = gatk_combine_variants(files_gatk.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                                                         '_filter2_indel_final.vcf_ANN.vcf.gz'),
                                                      args.reference, args.filter2_only_snp_vcf_dir, merged_file_suffix,
                                                      logger, Config)

    # Tabix index GATK combined Final vcf file
    files_for_tabix = glob.glob("%s/Final_vcf_*.vcf" % args.filter2_only_snp_vcf_dir)
    tabix(files_for_tabix, "vcf", logger, Config)


    # Extract ANN information from bcftools Final vcf file
    snp_var_ann_dict = {}
    indel_var_ann_dict = {}

    for variants in VCF("%s/Final_vcf_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir):
        snp_var_ann_dict[variants.POS] = variants.INFO.get('ANN')

    for variants in VCF("%s/Final_vcf_indel.vcf.gz" % args.filter2_only_snp_vcf_dir):
        indel_var_ann_dict[variants.POS] = variants.INFO.get('ANN')

    # ######Make changes in generating a new All_label_final_sorted.txt file that should have similar columns to Final_gatk vcf file
    # # Get all the unique variant positions from the matrix All_label_final_sorted/All_indel_label_final_sorted
    # position_label = OrderedDict()
    # with open("%s/All_label_final_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
    #     keep_logging('Reading All label positions file: %s/All_label_final_sorted.txt' % args.filter2_only_snp_vcf_dir, 'Reading All label positions file: %s/All_label_final_sorted.txt' % args.filter2_only_snp_vcf_dir, logger, 'info')
    #     csv_reader = csv.reader(csv_file, delimiter='\t')
    #     for row in csv_reader:
    #         position_label[row[0]] = ','.join(row[1:])
    # csv_file.close()
    #
    # position_indel_label = OrderedDict()
    # with open("%s/All_indel_label_final_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
    #     keep_logging('Reading All label positions file: %s/All_indel_label_final_sorted.txt' % args.filter2_only_snp_vcf_dir, 'Reading All label positions file: %s/All_indel_label_final_sorted.txt' % args.filter2_only_snp_vcf_dir, logger, 'info')
    #     csv_reader = csv.reader(csv_file, delimiter='\t')
    #     for row in csv_reader:
    #         if row[0] not in position_label.keys():
    #             position_indel_label[row[0]] = ','.join(row[1:])
    #         else:
    #             position_indel_label[row[0]] = ','.join(row[1:])
    #             keep_logging('Warning: position %s already present as a SNP' % row[0], 'Warning: position %s already present as a SNP' % row[0], logger, 'info')
    # csv_file.close()

    # Generate a header of file names
    print_string_header = "\t"
    for i in vcf_filenames:
        print_string_header = print_string_header + os.path.basename(i) + "\t"

    # Get the final core variant positions
    core_positions = []
    if ConfigSectionMap("functional_filters", Config)['apply_to_calls'] == "yes":
        core_positions_file = "%s/Only_ref_variant_positions_for_closely_without_functional_filtered_positions" % args.filter2_only_snp_vcf_dir
    else:
        core_positions_file = "%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir
    with open(core_positions_file) as fp:
        for line in fp:
            line = line.strip()
            core_positions.append(line)
        fp.close()

    indel_core_positions = []
    if ConfigSectionMap("functional_filters", Config)['apply_to_calls'] == "yes":
        core_positions_file = "%s/Only_ref_indel_variant_positions_for_closely_without_functional_filtered_positions" % args.filter2_only_snp_vcf_dir
    else:
        core_positions_file = "%s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir
    with open(core_positions_file) as fp:
        for line in fp:
            line = line.strip()
            indel_core_positions.append(line)
        fp.close()

    functional_filter_pos_array = []
    with open(functional_class_filter_positions, 'rU') as f_functional:
        for line_func in f_functional:
            functional_filter_pos_array.append(line_func.strip())

    # GET PHAGE/Repetitive region/mask region positions
    phage_positions = []
    repetitive_positions = []
    mask_positions = []
    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
        if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
            phage_region_positions = "%s/phage_region_positions.txt" % args.filter2_only_snp_vcf_dir
            if os.path.isfile(phage_region_positions):
                with open(phage_region_positions, 'rU') as fphage:
                    for line in fphage:
                        phage_positions.append(line.strip())
                fphage.close()

        # GET REPETITIVE REGIONS
        if ConfigSectionMap("functional_filters", Config)['find_repetitive_region'] == "yes":
            repetitive_positions_file = "%s/repeat_region_positions.txt" % args.filter2_only_snp_vcf_dir
            if os.path.isfile(repetitive_positions_file):
                with open(repetitive_positions_file, 'rU') as frep:
                    for line in frep:
                        repetitive_positions.append(line.strip())
                frep.close()

        # GET MASK REGIONS
        if ConfigSectionMap("functional_filters", Config)['mask_region'] == "yes":
            mask_positions_file = "%s/mask_positions.txt" % args.filter2_only_snp_vcf_dir
            if os.path.isfile(mask_positions_file):
                with open(mask_positions_file, 'rU') as fmask:
                    for line in fmask:
                        mask_positions.append(line.strip())
                fmask.close()


    # Prepare SNP/Indel Matrix strings to print

    # Read merged vcf files using cyvcf library
    final_merge_anno_file = VCF("%s/Final_vcf_gatk_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir)

    header_print_string = "Type of SNP at POS > ALT; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos"
    for sample in final_merge_anno_file.samples:
        # header_print_string = header_print_string + "," + sample
        header_print_string = header_print_string + "\t" + sample
    header_print_string = header_print_string + "\n"
    # header_print_string = header_print_string.replace(':::,', ':::')
    # changing to tab-delimited
    #header_print_string = header_print_string.replace(':::,', '\t')

    ######Make changes in generating a new All_label_final_sorted.txt file that should have similar columns to Final_gatk vcf file
    # Get all the unique variant positions from the matrix All_label_final_sorted/All_indel_label_final_sorted
    paste_label_command = "paste %s/unique_positions_file " % args.filter2_only_snp_vcf_dir
    paste_indel_label_command = "paste %s/unique_indel_positions_file " % args.filter2_only_snp_vcf_dir
    for filename_base in final_merge_anno_file.samples:
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        elif ".1.fastq.gz" in filename_base:
            second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
            first_part_split = filename_base.split('.1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
        sample_label_file = "%s/%s_filter2_final.vcf_no_proximate_snp.vcf_positions_label" % (args.filter2_only_snp_vcf_dir, first_part)
        sample_indel_label_file = "%s/%s_filter2_indel_final.vcf_indel_positions_label" % (args.filter2_only_snp_vcf_dir, first_part)
        paste_label_command = paste_label_command + sample_label_file + " "
        paste_indel_label_command = paste_indel_label_command + sample_indel_label_file + " "

    paste_label_command = paste_label_command + " > %s/All_label_final_ordered.txt" % args.filter2_only_snp_vcf_dir
    paste_indel_label_command = paste_indel_label_command + " > %s/All_indel_label_final_ordered.txt" % args.filter2_only_snp_vcf_dir
    sort_ordered_label_cmd = "sort -n -k1,1 %s/All_label_final_ordered.txt > %s/All_label_final_ordered_sorted.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    sort_ordered_indel_label_cmd = "sort -n -k1,1 %s/All_indel_label_final_ordered.txt > %s/All_indel_label_final_ordered_sorted.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    # call("%s" % paste_label_command, logger)
    # call("%s" % sort_ordered_label_cmd, logger)
    # call("%s" % paste_indel_label_command, logger)
    # call("%s" % sort_ordered_indel_label_cmd, logger)

    os.system(paste_label_command)
    os.system(sort_ordered_label_cmd)
    os.system(paste_indel_label_command)
    os.system(sort_ordered_indel_label_cmd)

    with open('%s/All_label_final_ordered.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
        outfile.write(paste_label_command + '\n')
        outfile.write(sort_ordered_label_cmd + '\n')
        outfile.write(paste_indel_label_command + '\n')
        outfile.write(sort_ordered_indel_label_cmd + '\n')
    outfile.close()

    os.system("bash %s/All_label_final_ordered.sh" % args.filter2_only_snp_vcf_dir)

    position_label = OrderedDict()
    with open("%s/All_label_final_ordered_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
        keep_logging('Reading All label positions file: %s/All_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
                     'Reading All label positions file: %s/All_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
                     logger, 'info')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            position_label[row[0]] = ','.join(row[1:])
    csv_file.close()

    position_indel_label = OrderedDict()
    with open("%s/All_indel_label_final_ordered_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
        keep_logging(
            'Reading All label positions file: %s/All_indel_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
            'Reading All label positions file: %s/All_indel_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
            logger, 'info')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if row[0] not in position_label.keys():
                position_indel_label[row[0]] = ','.join(row[1:])
            else:
                position_indel_label[row[0]] = ','.join(row[1:])
                keep_logging('Warning: position %s already present as a SNP' % row[0],
                             'Warning: position %s already present as a SNP' % row[0], logger, 'info')
    csv_file.close()

    # Open files to write the strings
    fp_code = open("%s/SNP_matrix_code.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele = open("%s/SNP_matrix_allele.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele_new = open("%s/SNP_matrix_allele_new.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_code.write(header_print_string)
    fp_allele.write(header_print_string)
    fp_allele_new.write(header_print_string)

    # Parse variant positions from the loaded cyvcf VCF object and generate the print strings
    for variants in VCF("%s/Final_vcf_gatk_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir):
        print_string = ""

        # Initiate and assign Functional Field filter string => PHAGE/REPEAT/MASK/NULL
        functional_field = ""
        if str(variants.POS) in phage_positions:
            functional_field = functional_field + "PHAGE_"
        else:
            functional_field = functional_field + "NULL_"
        if str(variants.POS) in repetitive_positions:
            functional_field = functional_field + "REPEATS_"
        else:
            functional_field = functional_field + "NULL_"
        if str(variants.POS) in mask_positions:
            functional_field = functional_field + "MASK"
        else:
            functional_field = functional_field + "NULL"

        # Initiate variant code string => REF allele = 0, core = 1, Filtered = 2, unmapped = -1, True but non-core = 3
        code_string = position_label[str(variants.POS)]
        code_string = code_string.replace('reference_allele', '0')
        code_string = code_string.replace('reference_unmapped_position', '-1')
        code_string = code_string.replace('LowFQ_QUAL_DP_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_DP_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_DP_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_QUAL_DP', '2')
        code_string = code_string.replace('LowFQ_DP_QUAL', '2')
        code_string = code_string.replace('LowFQ_QUAL', '2')
        code_string = code_string.replace('LowFQ_DP', '2')
        code_string = code_string.replace('HighFQ_QUAL_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_DP_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_QUAL_DP', '2')
        code_string = code_string.replace('HighFQ_DP_QUAL', '2')
        code_string = code_string.replace('HighFQ_QUAL', '2')
        code_string = code_string.replace('HighFQ_DP', '2')
        code_string = code_string.replace('LowFQ', '2')
        code_string = code_string.replace('HighFQ', '2')


        if str(variants.POS) in core_positions:
            code_string = code_string.replace('VARIANT', '1')
        # Adding functional class status code to SNP matrix: 2018-07-24
        elif str(variants.POS) in functional_filter_pos_array:
            code_string = code_string.replace('VARIANT', '2')
        else:
            code_string = code_string.replace('VARIANT', '3')

        # Assign type of snp: coding / non-coding
        if variants.INFO.get('ANN'):
            if "protein_coding" in variants.INFO.get('ANN'):
                snp_type = "Coding SNP"
            else:
                snp_type = "Non-coding SNP"
        else:
            if len(variants.ALT) > 1:
                #print variants.ALT
                #print ';'.join(set(snp_var_ann_dict[variants.POS].split(',')))
                if "protein_coding" in set(snp_var_ann_dict[variants.POS].split(',')):
                    snp_type = "Coding SNP"
                else:
                    snp_type = "Non-coding SNP"
            else:
                snp_type = "Non-coding SNP"

        print_string = print_string + snp_type + " at %s > " % str(variants.POS) + str(",".join(variants.ALT)) + " functional=%s" % functional_field

        # Get ANN field from variant INFO column and save it as an array. Split and Go through each elements, add bells and whistles
        if variants.INFO.get('ANN'):
            ann_array = (variants.INFO.get('ANN')).split(',')
            ann_string = ";"
            for i in list(set(ann_array)):
                i_split = i.split('|')
                #ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"
                tag = str(i_split[4]).replace('CHR_START-', '')
                tag = str(tag).replace('-CHR_END', '')

                if "-" in tag:
                    #print tag
                    extra_tags = ""
                    tag_split = tag.split('-')
                    for i in tag_split:
                        if i in locus_tag_to_gene_name.keys():
                            extra_tags = extra_tags + locus_tag_to_gene_name[i] + ","
                        else:
                            extra_tags = extra_tags + "None" + ","
                    extra_tags_prot = ""
                    for i in tag_split:
                        if i in locus_tag_to_product.keys():
                            extra_tags_prot = extra_tags_prot + locus_tag_to_product[i] + ","
                        else:
                            extra_tags_prot = extra_tags_prot + "None" + ","
                    ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags, extra_tags_prot]) + ";"
                elif tag == "":
                    print list(set(ann_array))
                else:
                    if tag in locus_tag_to_gene_name.keys() and tag in locus_tag_to_product.keys():
                        extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                    else:
                        print "tag key not found: %s" % tag
                        extra_tags = "NULL" + "|" + "NULL"
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                    ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
        else:
            if len(variants.ALT) > 1:
                #print variants.ALT
                #print ';'.join(set(snp_var_ann_dict[variants.POS].split(',')))
                ann_string = ";%s" % ';'.join(set(snp_var_ann_dict[variants.POS].split(',')))
            else:
                ann_string = ";None"


        print_string = print_string + ann_string

        gt_string = ""
        for gt in variants.gt_bases:
            gt = gt.replace('./.', '.')
            gt_string = gt_string + "," + gt
        gt_string = gt_string.replace('A/A', 'A')
        gt_string = gt_string.replace('G/G', 'G')
        gt_string = gt_string.replace('C/C', 'C')
        gt_string = gt_string.replace('T/T', 'T')
        gt_string = gt_string.replace('.', variants.REF)

        # #print print_string + gt_string + '\n'
        # fp_allele.write(print_string + gt_string + '\n')
        # #print print_string + "," + code_string + '\n'
        # fp_code.write(print_string + "," + code_string + '\n')

        final_allele_string = print_string + gt_string.replace(',', '\t') + '\n'
        final_code_string = print_string + "\t" + code_string.replace(',', '\t') + '\n'
        final_allele_string = final_allele_string.replace(',|', '|')
        # final_allele_string = final_allele_string.replace(',;,', ':::')
        # final_allele_string = final_allele_string.replace(';,', ':::')
        final_allele_string = final_allele_string.replace(',;,', ':::')
        final_allele_string = final_allele_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(',|', '|')
        # final_code_string = final_code_string.replace(',;,', ':::')
        # final_code_string = final_code_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(',;,', ':::')
        final_code_string = final_code_string.replace(';,', ':::')
        fp_allele.write(final_allele_string)
        fp_code.write(final_code_string)

        # print code_string
        # print gt_string[1:]
        ntd_string = ""
        count = 0
        code_string_array = code_string.split(',')
        gt_string_array = gt_string[1:].split(',')
        # print code_string_array
        # print gt_string_array
        for i in gt_string_array:
            # Fixing this: 2018-07-21
            if str(code_string_array[count]) == "0" or str(code_string_array[count]) == "1" or str(code_string_array[count]) == "3":
                ntd_string = ntd_string + "\t" + str(i)
            if code_string_array[count] == "-1":
                ntd_string = ntd_string + "\t" + "-"
            # Fixing this: 2018-07-21
            if str(code_string_array[count]) == "2":
                ntd_string = ntd_string + "\t" + "N"
            count += 1
        print_string = print_string + ntd_string + "\n"
        print_string.replace(',;,', '\t')
        print_string.replace(';,', '\t')
        fp_allele_new.write(print_string)

    fp_code.close()
    fp_allele.close()
    fp_allele_new.close()


    ## Indel SNP matrix generation steps
    header_print_string = "Type of SNP at POS > ALT; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos"
    final_merge_anno_file = VCF("%s/Final_vcf_gatk_indel.vcf.gz" % args.filter2_only_snp_vcf_dir)
    for sample in final_merge_anno_file.samples:
        # header_print_string = header_print_string + "," + sample
        header_print_string = header_print_string + "\t" + sample
    header_print_string = header_print_string + "\n"
    #header_print_string = header_print_string.replace(':::,', ':::')
    #header_print_string = header_print_string.replace(':::,', '\t')
    fp_code = open("%s/Indel_matrix_code.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele = open("%s/Indel_matrix_allele.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_code.write(header_print_string)
    fp_allele.write(header_print_string)

    for variants in VCF("%s/Final_vcf_gatk_indel.vcf.gz" % args.filter2_only_snp_vcf_dir):
        print_string = ""

        functional_field = ""
        if str(variants.POS) in phage_positions:
            functional_field = functional_field + "PHAGE_"
        else:
            functional_field = functional_field + "NULL_"
        if str(variants.POS) in repetitive_positions:
            functional_field = functional_field + "REPEATS_"
        else:
            functional_field = functional_field + "NULL_"
        if str(variants.POS) in mask_positions:
            functional_field = functional_field + "MASK"
        else:
            functional_field = functional_field + "NULL"

        code_string = position_indel_label[str(variants.POS)]
        code_string = code_string.replace('reference_allele', '0')
        code_string = code_string.replace('reference_unmapped_position', '-1')
        code_string = code_string.replace('LowFQ_QUAL_DP_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_DP_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_DP_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_proximate_SNP', '2')
        code_string = code_string.replace('LowFQ_QUAL_DP', '2')
        code_string = code_string.replace('LowFQ_DP_QUAL', '2')
        code_string = code_string.replace('LowFQ_QUAL', '2')
        code_string = code_string.replace('LowFQ_DP', '2')
        code_string = code_string.replace('HighFQ_QUAL_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_DP_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_QUAL_DP', '2')
        code_string = code_string.replace('HighFQ_DP_QUAL', '2')
        code_string = code_string.replace('HighFQ_QUAL', '2')
        code_string = code_string.replace('HighFQ_DP', '2')
        code_string = code_string.replace('LowFQ', '2')
        code_string = code_string.replace('HighFQ', '2')

        if str(variants.POS) in indel_core_positions:
            code_string = code_string.replace('VARIANT', '1')
        # Adding functional class status code to SNP matrix: 2018-07-24
        elif str(variants.POS) in functional_filter_pos_array:
            code_string = code_string.replace('VARIANT', '2')
        else:
            code_string = code_string.replace('VARIANT', '3')





        if variants.INFO.get('ANN'):
            if "protein_coding" in variants.INFO.get('ANN'):
                snp_type = "Coding INDEL"
            else:
                snp_type = "Non-coding INDEL"
        else:
            if len(variants.ALT) > 1:
                #print variants.ALT
                #print ';'.join(set(snp_var_ann_dict[variants.POS].split(',')))
                if "protein_coding" in set(indel_var_ann_dict[variants.POS].split(',')):
                    snp_type = "Coding INDEL"
                else:
                    snp_type = "Non-coding INDEL"
            else:
                snp_type = "Non-coding INDEL"
        print_string = print_string + snp_type + " at %s > " % str(variants.POS) + str(",".join(variants.ALT)) + " functional=%s" % functional_field

        if variants.INFO.get('ANN'):
            ann_array = (variants.INFO.get('ANN')).split(',')
            ann_string = ";"
            # for i in list(set(ann_array)):
            #     i_split = i.split('|')
            #     ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"
            # print_string = print_string + ann_string
            for i in list(set(ann_array)):
                i_split = i.split('|')
                #ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"
                tag = str(i_split[4]).replace('CHR_START-', '')
                tag = str(tag).replace('-CHR_END', '')
                tag = str(tag).replace('&', '-')
                #print tag
                if "-" in tag:
                    #print tag
                    extra_tags = ""
                    tag_split = tag.split('-')
                    for i in tag_split:
                        if i in locus_tag_to_gene_name.keys():
                            extra_tags = extra_tags + locus_tag_to_gene_name[i] + ","
                        else:
                            extra_tags = extra_tags + "None" + ","
                    extra_tags_prot = ""
                    for i in tag_split:
                        if i in locus_tag_to_product.keys():
                            extra_tags_prot = extra_tags_prot + locus_tag_to_product[i] + ","
                        else:
                            extra_tags_prot = extra_tags_prot + "None" + ","
                    ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags, extra_tags_prot]) + ";"
                else:
                    if tag in locus_tag_to_gene_name.keys() and tag in locus_tag_to_product.keys():
                        extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                    else:
                        print "tag key not found: %s" % tag
                        extra_tags = "NULL" + "|" + "NULL"
                    #extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                    ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
        else:

            if len(variants.ALT) > 1:
                #print variants.ALT
                #print ';'.join(set(indel_var_ann_dict[variants.POS].split(',')))
                ann_string = ";%s" % ';'.join(set(indel_var_ann_dict[variants.POS].split(',')))
            else:
                ann_string = ";None"

        print_string = print_string + ann_string

        gt_string = ""
        for gt in variants.gt_bases:
            gt = gt.replace('./.', '.')
            if "/" in gt:
                gt_split = gt.split('/')
                gt = gt_split[0]
            gt_string = gt_string + "," + gt
        gt_string = gt_string.replace('.', variants.REF)



        final_allele_string = print_string + gt_string.replace(',', '\t') + '\n'
        final_code_string = print_string + "\t" + code_string.replace(',', '\t') + '\n'
        final_allele_string = final_allele_string.replace(',|', '|')
        # final_allele_string = final_allele_string.replace(',;,', ':::')
        # final_allele_string = final_allele_string.replace(';,', ':::')
        final_allele_string = final_allele_string.replace(',;,', ':::')
        final_allele_string = final_allele_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(',|', '|')
        # final_code_string = final_code_string.replace(',;,', ':::')
        # final_code_string = final_code_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(',;,', ':::')
        final_code_string = final_code_string.replace(';,', ':::')
        fp_allele.write(final_allele_string)
        fp_code.write(final_code_string)
    fp_code.close()
    fp_allele.close()

def core_prep_snp(core_vcf_fasta_dir):
    """ Generate SNP Filter Label Matrix """
    generate_paste_command()

    """ Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix. """
    generate_position_label_data_matrix()

    """ Generate VCF files from final list of variants in Only_ref_variant_positions_for_closely; generate commands for consensus generation """
    generate_vcf_files()

    """ Generate consensus fasta file from core vcf files """
    extract_only_ref_variant_fasta_from_reference()

    """ Generate consensus fasta file with only reference and variant position bases """
    extract_only_ref_variant_fasta(core_vcf_fasta_dir)

    """ Analyze the positions that were filtered out only due to insufficient depth"""
    DP_analysis()

def core_prep_indel(core_vcf_fasta_dir):
    """ Generate SNP Filter Label Matrix """
    generate_indel_paste_command()

    """ Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix. """
    generate_indel_position_label_data_matrix()

""" report methods """
def alignment_report(data_matrix_dir):
    keep_logging('Generating Alignment report...', 'Generating Alignment report...', logger, 'info')
    varcall_dir = os.path.dirname(args.results_dir)
    print varcall_dir
    report_string = ""
    header = "Sample,QC-passed reads,Mapped reads,% mapped reads,mean depth,%_bases_above_5,%_bases_above_10,%_bases_above_15,unmapped_positions,READ_PAIR_DUPLICATES,READ_PAIR_OPTICAL_DUPLICATES,unmapped reads,% unmapped reads"
    fp = open("%s/Report_alignment.txt" % (data_matrix_dir), 'w+')
    fp.write(header + '\n')
    for vcf in vcf_filenames:
        sample = os.path.basename(vcf.replace('_filter2_final.vcf_no_proximate_snp.vcf', ''))
        #print sample
        report_string = sample + ","
        qc = (subprocess.check_output("grep \'QC-passed\' %s/%s/%s_alignment_stats | sed \'s/ + 0 in total (QC-passed reads + QC-failed reads)//g\'" % (varcall_dir, sample, sample), shell=True)).strip()
        mapped = (subprocess.check_output("grep \'mapped (\' %s/%s/%s_alignment_stats | awk -F\' \' \'{print $1}\'" % (varcall_dir, sample, sample), shell=True)).strip()
        replace = "%:-nan%)"
        perc_mapped = (subprocess.check_output("grep \'mapped (\' %s/%s/%s_alignment_stats | awk -F\' \' \'{print $5}\' | sed \'s/%s//g\' | sed \'s/(//g\'" % (varcall_dir, sample, sample, replace), shell=True)).strip()
        depth_of_coverage = (subprocess.check_output("awk -F\'\\t\' \'{OFS=\",\"};FNR==2{print $3,$7,$8,$9}\' %s/%s/%s_depth_of_coverage.sample_summary" % (varcall_dir, sample, sample), shell=True)).strip()
        unmapped_positions = (subprocess.check_output("wc -l %s/%s/%s_unmapped.bed_positions | cut -d\' \' -f1" % (varcall_dir, sample, sample), shell=True)).strip()
        opt_dup = (subprocess.check_output("awk -F\'\\t\' \'{OFS=\",\"};FNR==8{print $7,$8,$5}\' %s/%s/%s_markduplicates_metrics" % (varcall_dir, sample, sample), shell=True)).strip()
        perc_unmapped = str(100 - float(perc_mapped))
        myList = ','.join(map(str, (sample, qc, mapped, perc_mapped, depth_of_coverage, unmapped_positions, opt_dup, perc_unmapped)))
        #print myList
        fp.write(myList + '\n')
    fp.close()
    keep_logging('Alignment report can be found in %s/Report_alignment.txt' % data_matrix_dir, 'Alignment report can be found in %s/Report_alignment.txt' % data_matrix_dir, logger, 'info')

def variant_report(data_matrix_dir):
    keep_logging('Generating Variants report...', 'Generating Variants report...', logger, 'info')
    varcall_dir = os.path.dirname(os.path.abspath(args.results_dir))
    report_string = ""
    header = "Sample,Total Unique Variants,core SNPs,unmapped_positions,reference_allele,true_variant,Only_low_FQ,Only_DP,Only_low_MQ,other,unmapped_positions_perc,true_variant_perc,Only_low_FQ_perc,Only_DP_perc,Only_low_MQ_perc,other_perc"
    fp = open("%s/Report_variants.txt" % (data_matrix_dir), 'w+')
    fp.write(header + '\n')

    for vcf in vcf_filenames:
        sample = os.path.basename(vcf.replace('_filter2_final.vcf_no_proximate_snp.vcf', ''))
        report_string = sample + ","
        unmapped_positions = (subprocess.check_output("wc -l %s/core_temp_dir/unique_positions_file | cut -d\' \' -f1" % (varcall_dir), shell=True)).strip()
        core_snps = (subprocess.check_output("wc -l %s/core_temp_dir/Only_ref_variant_positions_for_closely | cut -d\' \' -f1" % (varcall_dir), shell=True)).strip()
        filtered_snp_count = (subprocess.check_output("grep -w \'^%s\' %s/core_temp_dir/bargraph_counts.txt | awk -F\'\\t\' \'{OFS=\",\"};{print $2,$3,$4,$5,$6,$7}\'" % (sample, varcall_dir), shell=True)).strip()
        filtered_snp_perc = (subprocess.check_output("grep -w \'^%s\' %s/core_temp_dir/bargraph_percentage.txt | awk -F\'\\t\' \'{OFS=\",\"};{print $2,$3,$4,$5,$6,$7}\'" % (sample, varcall_dir), shell=True)).strip()
        myList = ','.join(map(str, (sample, unmapped_positions, core_snps, filtered_snp_count, filtered_snp_perc)))
        fp.write(myList + '\n')
    fp.close()
    keep_logging('Variant call report can be found in %s/Report_variants.txt' % data_matrix_dir, 'Variant call report can be found in %s/Report_variants.txt' % data_matrix_dir, logger, 'info')

def gubbins(gubbins_dir, input_fasta, jobrun, logger, Config):
    keep_logging('\nRunning Gubbins on input: %s\n' % input_fasta, '\nRunning Gubbins on input: %s\n' % input_fasta,
                 logger,
                 'info')



    #gubbins_cmd = "%s/%s --prefix %s/%s %s" % (
    # ConfigSectionMap("gubbins", Config)['gubbins_bin'], ConfigSectionMap("gubbins", Config)['base_cmd'], gubbins_dir,
    # (os.path.basename(input_fasta)).replace('.fa', ''), input_fasta)


    gubbins_cmd = "%s --prefix %s/%s %s" % (
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
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=76:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], gubbins_dir, gubbins_cmd)
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        #os.system("qsub %s" % job_file_name)
        call("qsub %s" % job_file_name, logger)


"""
Pending inclusion

class FuncThread(threading.Thread):
    def __init__(self, target, *args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)
    def run(self):
        self._target(*self._args)

def someOtherFunc(data, key):
    print "someOtherFunc was called : data=%s; key=%s" % (str(data), str(key))

Pending inclusion
"""



if __name__ == '__main__':
    """
    Main Function for Variant Calling Core Pipeline 
    :param:
    :return:
    
    This function runs "core_prep" step to generate intermediate files required for extracting core variants at "core" step. 
    Using these core variants, a "report" step will generate the final reports and output results of the pipeline as well as runs "tree" step to generate fasttree and raxml results 
    using the core variants consensus in Date_Time_core_results folder.
    Steps:
    1. core_prep
    2. core
    3. report
    4. tree
    """


    # Start Timer to use it for generating folder names and Log prefixes.
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    global logger
    analysis_name_log = "step_" + str(args.steps)
    logger = generate_logger(args.filter2_only_snp_vcf_dir, analysis_name_log, log_unique_time)
    keep_logging('test log', 'test log', logger, 'info')
    keep_logging('\nThe Script started at: %s' % start_time, '\nThe Script started at: %s' % start_time, logger, 'info')
    print_details = "This step will parse final vcf files(*_no_proximate_snp.vcf) generated at the end of Variant Calling Pipeline. At the end of this step, the following results will be generated and placed in output directory:\n\n" \
          "1. Final Core SNP Positions list(Variant positions that were not filtered out in any of the samples and passed all the filters)\n" \
          "2. SNP Positions that were filtered out with labels indicating the reason (Depth, FQ, MQ, Unmapped in one or other samples, Proximate SNPS, Quality of Variant) why they were filtered out.\n" \
          "3. Barplot Statistics about the filtered variants and their reason for getting filtered.\n" \
          "4. Final Consensus fasta file using only Core SNP Positions\n"
    keep_logging('%s' % print_details, '%s' % print_details, logger, 'info')

    # Create temporary Directory core_temp_dir/temp for storing temporary intermediate files. Check if core_temp_dir contains all the required files to run these pipeline.
    global temp_dir
    temp_dir = args.filter2_only_snp_vcf_dir + "/temp"

    # Read Config file into Config object that will be used to extract configuration settings set up in config file.
    global config_file
    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
    global Config
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    keep_logging('Path to config file: %s' % config_file, 'Path to config file: %s' % config_file, logger, 'info')

    make_sure_path_exists(temp_dir)

    # Read filenames. Core variants and final results will be extracted considering only these files.
    filter2_only_snp_vcf_filenames = args.filter2_only_snp_vcf_filenames
    vcf_filenames_temp = []
    with open(filter2_only_snp_vcf_filenames) as fp:
        for line in fp:
            line = line.strip()
            line = args.filter2_only_snp_vcf_dir + line
            vcf_filenames_temp.append(line)
        fp.close()
    vcf_filenames = sorted(vcf_filenames_temp)

    make_sure_files_exists(vcf_filenames, Config, logger)





    log_file_handle = "%s/%s_%s.log.txt" % (args.filter2_only_snp_vcf_dir, log_unique_time, analysis_name_log)


    # Start Variant Calling Core Pipeline steps based on steps argument supplied.
    if "1" in args.steps:
        """ 
        core_prep step
        """

        # Gather SNP positions from each final *_no_proximate_snp.vcf file (that passed the variant filter parameters from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files for use in downstream methods
        keep_logging('Gathering SNP position information from each final *_no_proximate_snp.vcf file...', 'Gathering SNP position information from each final *_no_proximate_snp.vcf file...', logger, 'info')

        # Extract All the unique SNO and Indel position list from final filtered *_no_proximate_snp.vcf files.
        unique_position_file = create_positions_filestep(vcf_filenames)
        unique_indel_position_file = create_indel_positions_filestep(vcf_filenames)

        # bgzip and tabix all the vcf files in core_temp_dir.
        files_for_tabix = glob.glob("%s/*.vcf" % args.filter2_only_snp_vcf_dir)
        tabix(files_for_tabix, "vcf", logger, Config)

        if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
            keep_logging('Functional class filter is set to yes. Preparing Functional class filters\n', 'Functional class filter is set to yes. Preparing Functional class filters\n', logger,
                         'info')
            if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
                # Submit Phaster jobs to find ProPhage region in reference genome.
                run_phaster(args.reference, args.filter2_only_snp_vcf_dir, logger, Config)

        # Get the cluster option; create and run jobs based on given parameter. The jobs will parse all the intermediate vcf file to extract information such as if any unique variant position was unmapped in a sample, if it was filtered out dur to DP,MQ, FQ, proximity to indel, proximity to other SNPs and other variant filter parameters set in config file.
        tmp_dir = "/tmp/temp_%s/" % log_unique_time
        create_job(args.jobrun, vcf_filenames, unique_position_file, tmp_dir)
        create_indel_job(args.jobrun, vcf_filenames, unique_indel_position_file, tmp_dir)

        call("cp %s %s/Logs/core_prep/" % (log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)
    if "2" in args.steps:
        """ 
        core step 
        """

        # Set variables; check if the output from core_prep steps (*label files) exists and was completed without any errors.
        snp_unique_positions_file = args.filter2_only_snp_vcf_dir + "/unique_positions_file"
        indel_unique_positions_file = args.filter2_only_snp_vcf_dir + "/unique_indel_positions_file"
        uniq_snp_positions = sum(1 for line in open('%s' % snp_unique_positions_file))
        uniq_indel_positions = sum(1 for line in open('%s' % indel_unique_positions_file))
        if not os.path.isfile(snp_unique_positions_file) and not os.path.isfile(indel_unique_positions_file):
            keep_logging('Error finding unique_positions_file/unique_indel_positions_file. Please rerun core_prep step.','Error finding unique_positions_file/unique_indel_positions_file. Please rerun core_prep step.', logger,'exception')
            exit()
        make_sure_label_files_exists(vcf_filenames, uniq_snp_positions, uniq_indel_positions, Config, logger)

        # Set up Report and results directories to transfer the final results.
        data_matrix_dir = args.results_dir + '/data_matrix'
        core_vcf_fasta_dir = args.results_dir + '/core_snp_consensus'
        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(core_vcf_fasta_dir)

        functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % args.filter2_only_snp_vcf_dir

        # 06/05 Moving this to variant_call.py
        # #Parse Phaster results file to extract phage region.
        # if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
        #     keep_logging('Preparing Functional class filters\n', 'Preparing Functional class filters\n', logger,
        #                  'info')
        #     functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % args.filter2_only_snp_vcf_dir
        #     f1 = open(functional_class_filter_positions, 'w+')
        #     if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
        #         phage_region_positions = parse_phaster(args.reference, args.filter2_only_snp_vcf_dir, logger, Config)
        #         with open(phage_region_positions, 'rU') as fp:
        #             for line in fp:
        #                 f1.write(line)
        #         fp.close()
        #     if ConfigSectionMap("functional_filters", Config)['find_repetitive_region'] == "yes":
        #         # Find repeat regions in reference genome
        #         repeat_region_positions = nucmer_repeat(args.reference, args.filter2_only_snp_vcf_dir, logger, Config)
        #         with open(repeat_region_positions, 'rU') as fp:
        #             for line in fp:
        #                 f1.write(line)
        #         fp.close()
        #     if ConfigSectionMap("functional_filters", Config)['mask_region'] == "yes":
        #         # Mask custom region/Positions
        #         if ConfigSectionMap("functional_filters", Config)['mask_file']:
        #             mask_file = ConfigSectionMap("functional_filters", Config)['mask_file']
        #             mask_extension = os.path.splitext(mask_file)[1]
        #             if mask_extension == ".bed":
        #                 mask_positions_file = mask_regions(mask_file, args.filter2_only_snp_vcf_dir, logger, Config)
        #                 keep_logging(
        # 'Mask positions in this file %s will be filtered out' % mask_positions_file,
        # 'Mask positions in this file %s will be filtered out' % mask_positions_file,
        # logger, 'info')
        #             else:
        #                 #mask_positions_file = mask_file
        #                 os.system("cp %s %s/mask_positions.txt" % (mask_file, args.filter2_only_snp_vcf_dir))
        #                 mask_positions_file = "%s/mask_positions.txt" % args.filter2_only_snp_vcf_dir
        #                 keep_logging(
        # 'Mask positions in this file %s will be filtered out' % mask_positions_file,
        # 'Mask positions in this file %s will be filtered out' % mask_positions_file,
        # logger, 'info')
        #             with open(mask_positions_file   , 'rU') as fp:
        #                 for line in fp:
        #                     f1.write(line)
        #             fp.close()
        #     f1.close()

        # Run core steps. Generate SNP and data Matrix results. Extract core SNPS and consensus files.
        core_prep_indel(core_vcf_fasta_dir)


        core_prep_snp(core_vcf_fasta_dir)

        # Moving this up before core_prep_snp; for some weird reason, it is failing to generate Only_ref_indel
        #core_prep_indel(core_vcf_fasta_dir)

        # Annotate core variants. Generate SNP and Indel matrix.
        annotated_snp_matrix()

        # Read new allele matrix and generate fasta; generate a seperate function
        keep_logging('Generating Fasta from Variant Alleles...\n', 'Generating Fasta from Variant Alleles...\n', logger, 'info')

        create_job_allele_variant_fasta(args.jobrun, vcf_filenames, args.filter2_only_snp_vcf_dir, config_file)

        extract_only_ref_variant_fasta_from_reference_allele_variant()

        call("cp %s %s/Logs/core/" % (
            log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)

    if "3" in args.steps:
        """ 
        report step 
        """

        keep_logging('Step 3: Generate Reports and Results folder.', 'Step 3: Generate Reports and Results folder.', logger, 'info')

        # Generate DP barplots data and Analyze the FQ values of all the unique variant
        #DP_analysis_barplot()
        #FQ_analysis()

        # Set up Report and results directories to transfer the final results.
        data_matrix_dir = args.results_dir + '/data_matrix'
        data_matrix_snpeff_dir = data_matrix_dir + '/snpEff_results'
        core_vcf_fasta_dir = args.results_dir + '/core_snp_consensus'
        consensus_var_dir = core_vcf_fasta_dir + '/consensus_variant_positions'
        core_vcf_dir = core_vcf_fasta_dir + '/core_vcf'
        consensus_allele_var_dir = core_vcf_fasta_dir + '/consensus_allele_variant_positions'
        consensus_ref_allele_var_dir = core_vcf_fasta_dir + '/consensus_ref_allele_variant_positions'
        consensus_ref_var_dir = core_vcf_fasta_dir + '/consensus_ref_variant_positions'
        consensus_ref_allele_unmapped_variant_dir = core_vcf_fasta_dir + '/consensus_ref_allele_unmapped_variant'
        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(data_matrix_snpeff_dir)
        make_sure_path_exists(core_vcf_fasta_dir)
        make_sure_path_exists(consensus_var_dir)
        make_sure_path_exists(core_vcf_dir)
        make_sure_path_exists(consensus_allele_var_dir)
        make_sure_path_exists(consensus_ref_allele_var_dir)
        make_sure_path_exists(consensus_ref_var_dir)
        make_sure_path_exists(consensus_ref_allele_unmapped_variant_dir)
        reference_base = os.path.basename(args.reference).split('.')[0]
        # Move results to the results directory
        move_data_matrix_results = "cp -r %s/*.csv %s/*.txt %s/temp* %s/All* %s/Only* %s/*.R %s/R_scripts/generate_diagnostics_plots.R %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, os.path.dirname(os.path.abspath(__file__)), data_matrix_dir)
        #move_core_vcf_fasta_results = "cp %s/*_core.vcf.gz %s/*.fa %s/*_variants.fa %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, core_vcf_fasta_dir)
        move_core_vcf_fasta_results = "cp %s/*_core.vcf.gz %s/*.fa %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, core_vcf_fasta_dir)
        move_consensus_var_fasta_results = "mv %s/*_variants.fa %s/" % (core_vcf_fasta_dir, consensus_var_dir)
        move_consensus_ref_var_fasta_results = "mv %s/*.fa %s/" % (core_vcf_fasta_dir, consensus_ref_var_dir)
        move_core_vcf = "mv %s/*_core.vcf.gz %s/" % (core_vcf_fasta_dir, core_vcf_dir)
        move_consensus_allele_var_fasta_results = "mv %s/*allele_variants.fa %s/" % (consensus_var_dir, consensus_allele_var_dir)
        move_consensus_ref_allele_var_fasta_results = "mv %s/*_ref_allele_variants.fa %s/" % (
            consensus_allele_var_dir, consensus_ref_allele_var_dir)
        move_consensus_ref_allele_unmapped_var_fasta_results = "mv %s/*_ref_allele_unmapped_variants.fa %s/" % (
            consensus_var_dir, consensus_ref_allele_unmapped_variant_dir)
        move_snpeff_results = "mv %s/*ANN* %s/" % (data_matrix_dir, data_matrix_snpeff_dir)
        copy_reference = "cp %s %s/%s.fa" % (args.reference, consensus_ref_var_dir, reference_base)
        copy_reference_2 = "cp %s %s/%s.fa" % (args.reference, consensus_ref_allele_var_dir, reference_base)

        call("%s" % move_data_matrix_results, logger)
        call("%s" % move_core_vcf_fasta_results, logger)
        call("%s" % move_consensus_var_fasta_results, logger)
        call("%s" % move_consensus_ref_var_fasta_results, logger)
        call("%s" % move_core_vcf, logger)
        call("%s" % move_consensus_allele_var_fasta_results, logger)
        call("%s" % move_consensus_ref_allele_var_fasta_results, logger)
        call("%s" % move_consensus_ref_allele_unmapped_var_fasta_results, logger)
        call("%s" % copy_reference, logger)
        call("%s" % copy_reference_2, logger)
        call("%s" % move_snpeff_results, logger)
        subprocess.call(["sed -i 's/title_here/%s/g' %s/generate_diagnostics_plots.R" % (os.path.basename(args.results_dir), data_matrix_dir)], shell=True)

        # Sanity Check if the variant consensus files generated are of same length
        count = 0
        for line in open("%s/Only_ref_variant_positions_for_closely_matrix.txt" % data_matrix_dir).xreadlines():
            count += 1
            ref_variants = count - 1
        variant_consensus_files = glob.glob("%s/*_variants.fa" % core_vcf_fasta_dir)
        for f in variant_consensus_files:
            cmd2 = "%s/%s/bioawk -c fastx '{ print length($seq) }' < %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bioawk", Config)['bioawk_bin'], f)
            proc = subprocess.Popen([cmd2], stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()

            try:
                int(out2) != int(ref_variants)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    keep_logging('Error generating variant consensus position file: %s' % f,
                                 'Error generating variant consensus position file: %s' % f, logger, 'info')
                    keep_logging('Error generating variant consensus position file: %s' % f, 'Error generating variant consensus position file: %s' % f, logger, 'exception')
                    exit()

        # """ Generate alignment report """
        # alignment_report(data_matrix_dir)
        #
        # """ Generate core snps report """
        # variant_report(data_matrix_dir)

        """ Generating Gubbins MFA files"""
        #parse_phaster(args.reference)
        reference_base = os.path.basename(args.reference).split('.')[0]
        gubbins_dir = args.results_dir + '/gubbins'
        tree_dir = args.results_dir + '/trees'

        make_sure_path_exists(gubbins_dir)
        make_sure_path_exists(tree_dir)


        prepare_ref_var_consensus_input = "%s/gubbins/%s_%s_ref_var_consensus.fa" % (args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''), reference_base)
        prepare_var_consensus_input = "%s/gubbins/%s_%s_var_consensus.fa" % (args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''), reference_base)
        prepare_allele_var_consensus_input = "%s/gubbins/%s_%s_allele_var_consensus.fa" % (
        args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''),
        reference_base)
        prepare_ref_allele_var_consensus_input = "%s/gubbins/%s_%s_ref_allele_var_consensus.fa" % (
            args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''),
            reference_base)
        prepare_ref_allele_unmapped_consensus_input = "%s/gubbins/%s_%s_ref_allele_unmapped_consensus.fa" % (args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''), reference_base)
        prepare_ref_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_ref_variant_positions/*.fa > %s" % (args.results_dir, prepare_ref_var_consensus_input)
        prepare_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_variant_positions/*_variants.fa > %s" % (args.results_dir, prepare_var_consensus_input)
        prepare_allele_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_allele_variant_positions/*_allele_variants.fa > %s" % (
        args.results_dir, prepare_allele_var_consensus_input)
        prepare_ref_allele_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_ref_allele_variant_positions/*.fa > %s" % (
            args.results_dir, prepare_ref_allele_var_consensus_input)
        prepare_ref_allele_unmapped_consensus_input_cmd = "cat %s %s/core_snp_consensus/consensus_ref_allele_unmapped_variant/*.fa > %s" % (args.reference, args.results_dir, prepare_ref_allele_unmapped_consensus_input)
        call("%s" % prepare_ref_var_consensus_input_cmd, logger)
        call("%s" % prepare_var_consensus_input_cmd, logger)
        call("%s" % prepare_allele_var_consensus_input_cmd, logger)
        call("%s" % prepare_ref_allele_var_consensus_input_cmd, logger)
        call("%s" % prepare_ref_allele_unmapped_consensus_input_cmd, logger)
        # os.system(prepare_ref_var_consensus_input_cmd)
        # os.system(prepare_var_consensus_input_cmd)

        print_details = "Results for core pipeline can be found in: %s\n" \
              "Description of Results:\n" \
              "1. data_matrix folder contains all the data matrices and other temporary files generated during the core pipeline. bargraph_counts.txt and bargraph_percentage.txt: contains counts/percentage of unique positions filtered out due to different filter parameters for each sample. Run bargraph.R to plot bargraph statistics." \
              "2. core_snp_consensus contains all the core vcf and fasta files. *_core.vcf.gz: core vcf files, *.fa and *_variants.fa: core consensus fasta file and core consensus fasta with only variant positions." % (args.results_dir)
        keep_logging(print_details, print_details, logger, 'info')

        call("cp %s %s/Logs/report/" % (
            log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)

    if "4" in args.steps:
        """ 
        tree step 
        """

        keep_logging('Step 4: Generate FastTree and RAxML trees for core consensus fasta files.', 'Step 4: Generate FastTree and RAxML trees for core consensus fasta files.', logger, 'info')

        #parse_phaster(args.reference)
        reference_base = os.path.basename(args.reference).split('.')[0]
        gubbins_dir = args.results_dir + '/gubbins'
        tree_dir = args.results_dir + '/trees'

        make_sure_path_exists(gubbins_dir)
        make_sure_path_exists(tree_dir)


        prepare_ref_var_consensus_input = "%s/gubbins/%s_%s_ref_var_consensus.fa" % (args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''), reference_base)
        prepare_var_consensus_input = "%s/gubbins/%s_%s_var_consensus.fa" % (args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''), reference_base)
        prepare_allele_var_consensus_input = "%s/gubbins/%s_%s_allele_var_consensus.fa" % (
        args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''),
        reference_base)
        prepare_ref_allele_var_consensus_input = "%s/gubbins/%s_%s_ref_allele_var_consensus.fa" % (
            args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''),
            reference_base)
        prepare_ref_allele_unmapped_consensus_input = "%s/gubbins/%s_%s_ref_allele_unmapped_consensus.fa" % (
        args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''),
        reference_base)
        prepare_ref_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_ref_variant_positions/*.fa > %s" % (args.results_dir, prepare_ref_var_consensus_input)
        prepare_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_variant_positions/*_variants.fa > %s" % (args.results_dir, prepare_var_consensus_input)
        prepare_allele_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_allele_variant_positions/*_allele_variants.fa > %s" % (
        args.results_dir, prepare_allele_var_consensus_input)
        prepare_ref_allele_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_ref_allele_variant_positions/*.fa > %s" % (
            args.results_dir, prepare_ref_allele_var_consensus_input)
        prepare_ref_allele_unmapped_consensus_input_cmd = "cat %s %s/core_snp_consensus/consensus_ref_allele_unmapped_variant/*.fa > %s" % (
        args.reference, args.results_dir, prepare_ref_allele_unmapped_consensus_input)
        call("%s" % prepare_ref_var_consensus_input_cmd, logger)
        call("%s" % prepare_var_consensus_input_cmd, logger)
        call("%s" % prepare_allele_var_consensus_input_cmd, logger)
        call("%s" % prepare_ref_allele_var_consensus_input_cmd, logger)
        call("%s" % prepare_ref_allele_unmapped_consensus_input_cmd, logger)
        # os.system(prepare_ref_var_consensus_input_cmd)
        # os.system(prepare_var_consensus_input_cmd)

        fasttree(tree_dir, prepare_ref_var_consensus_input, args.jobrun, logger, Config)
        fasttree(tree_dir, prepare_var_consensus_input, args.jobrun, logger, Config)
        fasttree(tree_dir, prepare_allele_var_consensus_input, args.jobrun, logger, Config)
        fasttree(tree_dir, prepare_ref_allele_var_consensus_input, args.jobrun, logger, Config)

        raxml(tree_dir, prepare_ref_var_consensus_input, args.jobrun, logger, Config)
        raxml(tree_dir, prepare_var_consensus_input, args.jobrun, logger, Config)
        raxml(tree_dir, prepare_allele_var_consensus_input, args.jobrun, logger, Config)
        raxml(tree_dir, prepare_ref_allele_var_consensus_input, args.jobrun, logger, Config)

        if args.gubbins and args.gubbins == "yes":
            gubbins(gubbins_dir, prepare_ref_var_consensus_input, args.jobrun, logger, Config)
            #gubbins(gubbins_dir, prepare_ref_allele_var_consensus_input, logger, Config)
            gubbins(gubbins_dir, prepare_ref_allele_unmapped_consensus_input,args.jobrun, logger, Config)
        call("cp %s %s/Logs/tree/" % (
            log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)

    if "5" in args.steps:
        """ 
        Debugging Purposes only: Run only SNP matrix annotation step 
        """

        keep_logging('Step 5: Running SNP matrix annotation step.', 'Step 5: Running SNP matrix annotation step.', logger, 'info')

        functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % args.filter2_only_snp_vcf_dir

        # Annotate core variants. Generate SNP and Indel matrix.
        annotated_snp_matrix()

        # # Read new allele matrix and generate fasta; generate a seperate function
        keep_logging('Generating Fasta from Variant Alleles...\n', 'Generating Fasta from Variant Alleles...\n', logger, 'info')

        create_job_allele_variant_fasta(args.jobrun, vcf_filenames, args.filter2_only_snp_vcf_dir, config_file)

        extract_only_ref_variant_fasta_from_reference_allele_variant()

        call("cp %s %s/Logs/core/" % (
            log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)

    time_taken = datetime.now() - start_time_2
    if args.remove_temp:
        del_command = "rm -r %s" % temp_dir
        os.system(del_command)















