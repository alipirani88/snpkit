from __future__ import division
import argparse
import re
import os
import csv
import subprocess
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
from Bio import SeqIO

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
optional.add_argument('-debug_mode', action='store', dest="debug_mode",
                    help='yes/no for debug mode')
args = parser.parse_args()


def create_positions_filestep(vcf_filenames):

    """
    Gather SNP positions from each final *_no_proximate_snp.vcf file (that passed the variant filter parameters
    from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files for use in downstream methods

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
            position_array.append(line)
        f.close()
    position_array_unique = set(position_array)
    position_array_sort = sorted(position_array_unique)
    print "\nThe number of unique variant positions: " + str(len(position_array_sort)) + "\n"
    unique_position_file = "%s/unique_positions_file" % args.filter2_only_snp_vcf_dir
    f=open(unique_position_file, 'w+')
    for i in position_array_sort:
        f.write(i + "\n")
    f.close()
    if len(position_array_sort) == 0:
        print "ERROR: No unique positions found. Check if vcf files are empty?"
        exit()
    return unique_position_file


def create_indel_positions_filestep(vcf_filenames):

    """
    Gather SNP positions from each final *_no_proximate_snp.vcf file (that passed the variant filter parameters
    from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files for use in downstream methods

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
            position_array.append(line)
        f.close()
    position_array_unique = set(position_array)
    position_array_sort = sorted(position_array_unique)
    print "\nThe number of unique indel positions: " + str(len(position_array_sort)) + "\n"
    unique_indel_position_file = "%s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir
    f=open(unique_indel_position_file, 'w+')
    for i in position_array_sort:
        f.write(i + "\n")
    f.close()
    if len(position_array_sort) == 0:
        print "ERROR: No unique positions found. Check if vcf files are empty?"
        exit()
    return unique_indel_position_file



def make_sure_path_exists(out_path):
    """
    Make sure the output folder exists or create at given path
    :param out_path:
    :return:
    """
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()

def run_command(i):
    print "Running: %s" % i
    os.system(i)
    done = "done: %s" % i
    return done


def create_job(jobrun, vcf_filenames, unique_position_file, tmp_dir):

    """
    Based on type of jobrun; generate jobs and run accordingly.
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
            print "Running: qsub %s" % i
            os.system("qsub %s" % i)

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
        command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
        os.system("bash %s" % command_file)
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
        #print "Running local mode: bash %s" % command_file
        os.system("bash %s" % command_file)

def create_indel_job(jobrun, vcf_filenames, unique_position_file, tmp_dir):

    """
    Based on type of jobrun; generate jobs and run accordingly.
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
            print "Running: qsub %s" % i
            os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
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

    elif jobrun == "cluster":
        command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
        os.system("bash %s" % command_file)
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
        #print "Running local mode: bash %s" % command_file
        os.system("bash %s" % command_file)

def create_job_fasta(jobrun, vcf_filenames, core_vcf_fasta_dir):

    """
    Based on type of jobrun; generate jobs and run accordingly.
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
            job_file_name = "%s_fasta.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            print "Running: qsub %s" % i
            os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
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

    elif jobrun == "cluster":
        command_array = []
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
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
        os.system("bash %s/command_file" % args.filter2_only_snp_vcf_dir)
    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_array = []
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
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
        os.system("bash command_file")

def create_job_DP(jobrun, vcf_filenames):
    """
    Based on type of jobrun; generate jobs and run accordingly.
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
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            print "Running: qsub %s" % i
            os.system("qsub %s" % i)


    elif jobrun == "parallel-local":
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

    elif jobrun == "cluster":
        """ Test pending """
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
        os.system("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir)

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
        os.system("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir)

def generate_paste_command():

    """ Generate SNP Filter Label Matrix """
    paste_file = args.filter2_only_snp_vcf_dir + "/paste_label_files.sh"
    f4=open(paste_file, 'w+')
    paste_command = "paste %s/unique_positions_file" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
        paste_command = paste_command + " " + label_file
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

    os.system(header_awk_cmd)
    os.system(sed_header)
    os.system(sed_header_2)

    temp_paste_command = paste_command + " > %s/temp_label_final_raw.txt" % args.filter2_only_snp_vcf_dir
    paste_command = paste_command + " > %s/All_label_final_raw" % args.filter2_only_snp_vcf_dir
    f4.write(paste_command)
    f4.close()
    sort_All_label_cmd = "sort -n -k1,1 %s/All_label_final_raw > %s/All_label_final_sorted.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    paste_command_header = "cat %s/header.txt %s/All_label_final_sorted.txt > %s/All_label_final_sorted_header.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    #print temp_paste_command

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
    os.system("bash %s/All_label_final_raw.sh" % args.filter2_only_snp_vcf_dir)
    os.system("bash %s/temp_label_final_raw.txt.sh" % args.filter2_only_snp_vcf_dir)
    print "Finished pasting...DONE"

    """
    remove this lines
    #subprocess.call(["%s" % paste_command], shell=True)
    #subprocess.call(["%s" % temp_paste_command], shell=True)
    #subprocess.check_call('%s' % paste_command)
    #subprocess.check_call('%s' % temp_paste_command)
    #os.system(paste_command) change
    #os.system(temp_paste_command) change
    """
    os.system(sort_All_label_cmd)
    os.system(paste_command_header)
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
    os.system(remove_unwanted_text)


def generate_indel_paste_command():

    """ Generate SNP Filter Label Matrix """
    paste_file = args.filter2_only_snp_vcf_dir + "/paste_indel_label_files.sh"
    f4=open(paste_file, 'w+')
    paste_command = "paste %s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf_indel_positions_label')
        paste_command = paste_command + " " + label_file
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

    os.system(header_awk_cmd)
    os.system(sed_header)
    os.system(sed_header_2)

    temp_paste_command = paste_command + " > %s/temp_indel_label_final_raw.txt" % args.filter2_only_snp_vcf_dir
    paste_command = paste_command + " > %s/All_indel_label_final_raw" % args.filter2_only_snp_vcf_dir
    f4.write(paste_command)
    f4.close()
    sort_All_label_cmd = "sort -n -k1,1 %s/All_indel_label_final_raw > %s/All_indel_label_final_sorted.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    paste_command_header = "cat %s/header.txt %s/All_indel_label_final_sorted.txt > %s/All_indel_label_final_sorted_header.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    #print temp_paste_command

    ls = []
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf_indel_positions_label')
        ls.append(label_file)
    ls.insert(0, "%s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir)

    with open('%s/All_indel_label_final_raw.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
        outfile.write(paste_command)
    outfile.close()

    with open('%s/temp_indel_label_final_raw.txt.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
        outfile.write(temp_paste_command)
    outfile.close()
    os.system("bash %s/All_indel_label_final_raw.sh" % args.filter2_only_snp_vcf_dir)
    os.system("bash %s/temp_indel_label_final_raw.txt.sh" % args.filter2_only_snp_vcf_dir)
    print "Finished pasting...DONE"

    """
    remove this lines
    #subprocess.call(["%s" % paste_command], shell=True)
    #subprocess.call(["%s" % temp_paste_command], shell=True)
    #subprocess.check_call('%s' % paste_command)
    #subprocess.check_call('%s' % temp_paste_command)
    #os.system(paste_command) change
    #os.system(temp_paste_command) change
    """
    os.system(sort_All_label_cmd)
    os.system(paste_command_header)
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
    os.system(remove_unwanted_text)


def generate_position_label_data_matrix():

    """
    Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix.

    Filtered Position label matrix:
        Too bad! This positions where atleast one variant was observed in atleast one sample
        This position didn't made it to the final Only_ref_variant_positions_for_closely_matrix list,
        because it was either unmapped(non-core) in one or more of the samples or was filtered out one or more of the sample due to Variant Filtered Parameter

    Only_ref_variant_positions_for_closely_matrix.txt :
        Those Positions where the variant was either reference allele or a variant that passed all the variant filter parameters.
        Yeah! This ones made it to final vcf file and are core variants
        (Core variant Position: Variant Position which was not filtered out in any of the other samples due to variant filter parameter and also this position was present in all the samples(not unmapped)).

    """
    def generate_position_label_data_matrix_All_label():
        position_label = OrderedDict()
        f1=open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'w+')
        f2=open("%s/Only_ref_variant_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f3=open("%s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f4=open("%s/Only_filtered_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            print "Reading All label positions file: %s/All_label_final_sorted_header.txt \n" % args.filter2_only_snp_vcf_dir
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                position_label[row[0]] = row[1:]
            print "Generating different list of Positions and heatmap data matrix... \n"
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
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir], shell=True)
        subprocess.call(["sed -i 's/1TRUE/-1/g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir], shell=True)

    def temp_generate_position_label_data_matrix_All_label():

        """
        Read **temp_label_final_raw.txt** SNP position label data matrix for generating barplot statistics.
        """
        temp_position_label = OrderedDict()
        f33=open("%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        print_string_header = "\t"
        for i in vcf_filenames:
            print_string_header = print_string_header + os.path.basename(i) + "\t"
        f33.write('\t' + print_string_header.strip() + '\n')
        print "Reading temporary label positions file: %s/temp_label_final_raw.txt \n" % args.filter2_only_snp_vcf_dir
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
            print "Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n" % args.filter2_only_snp_vcf_dir
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
            print "Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n" % args.filter2_only_snp_vcf_dir
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
        print "\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n"
        """
        Read each Sample columns and calculate the percentage of each label to generate barplot statistics.
        This will give a visual explanation of how many positions in each samples were filtered out because of different reason
        """

        c_reader = csv.reader(open('%s/temp_Only_filtered_positions_for_closely_matrix.txt' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
        columns = list(zip(*c_reader))
        print "Finished reading columns..."
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
        bargraph_R_script = "library(ggplot2)\nlibrary(reshape)\nx1 <- read.table(\"bargraph_percentage.txt\", header=TRUE)\nx1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))\nmdf1=melt(x1,id.vars=\"Sample\")\npdf(\"barplot.pdf\", width = 30, height = 30)\nggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat=\"identity\") + ylab(\"Percentage of Filtered Positions\") + xlab(\"Samples\") + theme(text = element_text(size=9)) + scale_fill_manual(name=\"Reason for filtered out positions\", values=c(\"#08306b\", \"black\", \"orange\", \"darkgrey\", \"#fdd0a2\", \"#7f2704\")) + ggtitle(\"Title Here\") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = \"horizontal\")\ndev.off()"
        barplot_R_file = open("%s/bargraph.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        print "Run this R script to generate bargraph plot: %s" % barplot_R_file
    """ Methods Steps"""
    print "Running: Generating data matrices..."
    generate_position_label_data_matrix_All_label()
    print "Running: Changing variables in data matrices to codes for faster processing..."
    temp_generate_position_label_data_matrix_All_label()
    print "Running: Generating Barplot statistics data matrices..."
    barplot_stats()

def generate_indel_position_label_data_matrix():

    """
    Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix.

    Filtered Position label matrix:
        Too bad! This positions where atleast one variant was observed in atleast one sample
        This position didn't made it to the final Only_ref_variant_positions_for_closely_matrix list,
        because it was either unmapped(non-core) in one or more of the samples or was filtered out one or more of the sample due to Variant Filtered Parameter

    Only_ref_variant_positions_for_closely_matrix.txt :
        Those Positions where the variant was either reference allele or a variant that passed all the variant filter parameters.
        Yeah! This ones made it to final vcf file and are core variants
        (Core variant Position: Variant Position which was not filtered out in any of the other samples due to variant filter parameter and also this position was present in all the samples(not unmapped)).

    """
    def generate_indel_position_label_data_matrix_All_label():
        position_label = OrderedDict()
        f1=open("%s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'w+')
        f2=open("%s/Only_ref_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f3=open("%s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f4=open("%s/Only_filtered_indel_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
            print "Reading All label positions file: %s/All_indel_label_final_sorted_header.txt \n" % args.filter2_only_snp_vcf_dir
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                position_label[row[0]] = row[1:]
            print "Generating different list of Positions and heatmap data matrix... \n"
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
                    print "here"
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
        print "Reading temporary label positions file: %s/temp_label_final_raw.txt \n" % args.filter2_only_snp_vcf_dir
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
            print "Reading temporary Only_filtered_indel_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt \n" % args.filter2_only_snp_vcf_dir
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
            print "Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt \n" % args.filter2_only_snp_vcf_dir
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
        print "\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n"
        """
        Read each Sample columns and calculate the percentage of each label to generate barplot statistics.
        This will give a visual explanation of how many positions in each samples were filtered out because of different reason
        """

        c_reader = csv.reader(open('%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
        columns = list(zip(*c_reader))
        print "Finished reading columns..."
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
        bargraph_R_script = "library(ggplot2)\nlibrary(reshape)\nx1 <- read.table(\"bargraph_indel_percentage.txt\", header=TRUE)\nx1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))\nmdf1=melt(x1,id.vars=\"Sample\")\npdf(\"barplot.pdf\", width = 30, height = 30)\nggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat=\"identity\") + ylab(\"Percentage of Filtered Positions\") + xlab(\"Samples\") + theme(text = element_text(size=9)) + scale_fill_manual(name=\"Reason for filtered out positions\", values=c(\"#08306b\", \"black\", \"orange\", \"darkgrey\", \"#fdd0a2\", \"#7f2704\")) + ggtitle(\"Title Here\") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = \"horizontal\")\ndev.off()"
        barplot_R_file = open("%s/bargraph_indel.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        print "Run this R script to generate bargraph plot: %s/bargraph_indel.R" % args.filter2_only_snp_vcf_dir


    """ Methods Steps"""
    print "Running: Generating data matrices..."
    generate_indel_position_label_data_matrix_All_label()
    print "Running: Changing variables in data matrices to codes for faster processing..."
    temp_generate_indel_position_label_data_matrix_All_label()
    print "Running: Generating Barplot statistics data matrices..."
    barplot_indel_stats()


def generate_vcf_files():
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    filter2_files_array = []
    for i in vcf_filenames:
        filter2_file = i.replace('_no_proximate_snp.vcf', '')
        filter2_files_array.append(filter2_file)
    ref_variant_position_array = []
    ffp = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'r+')
    for line in ffp:
        line = line.strip()
        ref_variant_position_array.append(line)
    ffp.close()

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
        print "Generating %s" % file_name
        filtered_out_vcf_files.append(file_name)
        f1 = open(file_name, 'w+')
        for ios in print_array:
            print_string = str(ios) + "\n"
            f1.write(print_string)
        f1.close()

    filename = "%s/consensus.sh" % args.filter2_only_snp_vcf_dir
    print "\nGenerating Consensus...\n"
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
    print "The consensus commands are in : %s" % filename
    sequence_lgth_cmd = "for i in %s/*.fa; do %s/%s/bioawk -c fastx \'{ print $name, length($seq) }\' < $i; done" % (args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bioawk", Config)['bioawk_bin'])
    os.system(sequence_lgth_cmd)

def gatk_filter2(final_raw_vcf, out_path, analysis, reference):
    gatk_filter2_parameter_expression = "MQ > 50 && QUAL > 100 && DP > 9"
    gatk_filter2_command = "java -jar %s/%s/GenomeAnalysisTK.jar -T VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s --filterExpression \"%s\" --filterName PASS_filter2" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("gatk", Config)['gatk_bin'], reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
    print "\n\nRunning Command: [%s]\n\n" % gatk_filter2_command
    os.system(gatk_filter2_command)
    filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_gatk.vcf > %s/%s_filter2_final.vcf" % (out_path, analysis, out_path, analysis)
    os.system(filter_flag_command)
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
        os.system(grep_fq_field)
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
    os.system("bash %s/paste_DP_files.sh" % args.filter2_only_snp_vcf_dir)
    print "Generating DP barplots data..."
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
    create_job_fasta(args.jobrun, vcf_filenames, core_vcf_fasta_dir)

def extract_only_ref_variant_fasta_from_reference():
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
            print "Error extracting reference allele"

    pattern = re.compile(r'\s+')
    fasta_string = re.sub(pattern, '', fasta_string)
    final_fasta_string = ">%s\n" % os.path.basename(args.reference.replace('.fasta', '')) + fasta_string + "\n"
    fp = open("%s/%s_variants.fa" % (args.filter2_only_snp_vcf_dir, os.path.basename(args.reference.replace('.fasta', ''))), 'w+')
    fp.write(final_fasta_string)
    fp.close()

def prepare_snpEff_db(reference_basename):
    keep_logging('Preparing snpEff database requirements.', 'Preparing snpEff database requirements.', logger, 'info')
    os.system("cp %s/%s/snpEff.config %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], args.filter2_only_snp_vcf_dir))
    make_sure_path_exists("%s/%s/data/%s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]))
    make_sure_path_exists("%s/%s/data/genomes/" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']))
    os.system("cp %s %s/%s/data/genomes/" % (args.reference, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']))
    with open("%s/snpEff.config" % args.filter2_only_snp_vcf_dir, "a") as conf_file:
        conf_file.write("\n\n##Building Custom Database###\n%s.genome\t: %s\n\n" % (reference_basename[0], reference_basename[0]))
    conf_file.close()
    #get the gff name from config file
    os.system("cp %s/%s.gff %s/%s/data/%s/genes.gff" % (os.path.dirname(args.reference), reference_basename[0], ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], reference_basename[0]))
    os.system("java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']))

    keep_logging('Finished Preparing snpEff database requirements.', 'Finished Preparing snpEff database requirements.', logger, 'info')

def variant_annotation():
    keep_logging('Annotating Variants using snpEff.', 'Annotating Variants using snpEff.', logger, 'info')
    reference_basename = (os.path.basename(args.reference)).split(".")
    print reference_basename[0]
    prepare_snpEff_db(reference_basename)
    annotate_vcf_cmd_array = []
    annotate_final_vcf_cmd_array = []
    for i in vcf_filenames:
        raw_vcf = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf')
        annotate_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], raw_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, reference_basename[0], raw_vcf, raw_vcf)
        annotate_vcf_cmd_array.append(annotate_vcf_cmd)
        final_vcf = i
        annotate_final_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, reference_basename[0], final_vcf, final_vcf)
        annotate_final_vcf_cmd_array.append(annotate_final_vcf_cmd)
    if args.numcores:
        num_cores = int(num_cores)
    else:
        num_cores = multiprocessing.cpu_count()
    print "\n\nhere\n\n"
    print annotate_final_vcf_cmd_array
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_vcf_cmd_array)
    results_2 = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_final_vcf_cmd_array)

def indel_annotation():
    keep_logging('Annotating indels using snpEff.', 'Annotating indels using snpEff.', logger, 'info')
    reference_basename = (os.path.basename(args.reference)).split(".")
    print reference_basename[0]
    prepare_snpEff_db(reference_basename)
    annotate_vcf_cmd_array = []
    annotate_final_vcf_cmd_array = []
    for i in vcf_filenames:
        raw_vcf = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf')
        annotate_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], raw_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, reference_basename[0], raw_vcf, raw_vcf)
        annotate_vcf_cmd_array.append(annotate_vcf_cmd)
        final_vcf = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf')
        annotate_final_vcf_cmd = "java -Xmx4g -jar %s/%s/%s -csvStats %s_ANN.csv -dataDir %s/%s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir, reference_basename[0], final_vcf, final_vcf)
        annotate_final_vcf_cmd_array.append(annotate_final_vcf_cmd)
    if args.numcores:
        num_cores = int(num_cores)
    else:
        num_cores = multiprocessing.cpu_count()
    print annotate_final_vcf_cmd_array
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_vcf_cmd_array)
    results_2 = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_final_vcf_cmd_array)

def annotated_snp_matrix():
    """
    :return: Read Genbank file and return a dictionary of Prokka ID mapped to Gene Name, Prokka ID mapped to Product Name
    """
    reference_basename = (os.path.basename(args.reference)).split(".")
    handle = open("%s/%s.gbf" % (os.path.dirname(args.reference), reference_basename[0]), 'rU')
    locus_tag_to_gene_name = {}
    locus_tag_to_product = {}
    #locus_tag_to_uniprot = {}
    #locus_tag_to_ec_number = {}

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



    """ Merge Annotated final vcf file """
    print "Merging Final Annotated VCF files into %s/Final_vcf_no_proximate_snp.vcf" % args.filter2_only_snp_vcf_dir
    os.system("for i in %s/*.vcf_no_proximate_snp.vcf_ANN.vcf; do bgzip -c $i > $i.gz; done" % args.filter2_only_snp_vcf_dir)
    os.system("for i in %s/*.vcf_no_proximate_snp.vcf_ANN.vcf.gz; do tabix $i; done" % args.filter2_only_snp_vcf_dir)
    os.system("for i in %s/*_filter2_indel_final.vcf_ANN.vcf; do bgzip -c $i > $i.gz; done" % args.filter2_only_snp_vcf_dir)
    os.system("for i in %s/*_filter2_indel_final.vcf_ANN.vcf.gz; do tabix $i; done" % args.filter2_only_snp_vcf_dir)


    os.system("bcftools merge -i ANN:join -m both -o %s/Final_vcf_no_proximate_snp.vcf -O v %s/*.vcf_no_proximate_snp.vcf_ANN.vcf.gz" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
    os.system("bcftools merge -i ANN:join -m both -o %s/Final_vcf_indel.vcf -O v %s/*_filter2_indel_final.vcf_ANN.vcf.gz" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))

    os.system("bgzip -c %s/Final_vcf_no_proximate_snp.vcf > %s/Final_vcf_no_proximate_snp.vcf.gz" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
    os.system("tabix %s/Final_vcf_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir)
    os.system("bgzip -c %s/Final_vcf_indel.vcf > %s/Final_vcf_indel.vcf.gz" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
    os.system("tabix %s/Final_vcf_indel.vcf.gz" % args.filter2_only_snp_vcf_dir)



    position_label = OrderedDict()
    with open("%s/All_label_final_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
        print "Reading All label positions file: %s/All_label_final_sorted.txt \n" % args.filter2_only_snp_vcf_dir
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            position_label[row[0]] = ','.join(row[1:])
    csv_file.close()

    with open("%s/All_indel_label_final_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
        print "Reading All label positions file: %s/All_indel_label_final_sorted.txt \n" % args.filter2_only_snp_vcf_dir
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if row[0] not in position_label.keys():
                position_label[row[0]] = ','.join(row[1:])
            else:
                print "Warning: position %s already present as a SNP" % row[0]
    csv_file.close()

    print_string_header = "\t"
    for i in vcf_filenames:
        print_string_header = print_string_header + os.path.basename(i) + "\t"

    core_positions = []
    with open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir) as fp:
        for line in fp:
            line = line.strip()
            core_positions.append(line)
        fp.close()
    with open("%s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir) as fp:
        for line in fp:
            line = line.strip()
            core_positions.append(line)
        fp.close()


    header_print_string = "Type of SNP at POS > ALT; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos"
    final_merge_anno_file = VCF("%s/Final_vcf_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir)
    for sample in final_merge_anno_file.samples:
        header_print_string = header_print_string + "," + sample
    header_print_string = header_print_string + "\n"

    fp_code = open("%s/SNP_matrix_code.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele = open("%s/SNP_matrix_allele.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_code.write(header_print_string)
    fp_allele.write(header_print_string)

    for variants in VCF("%s/Final_vcf_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir):
        print_string = ""

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
        else:
            code_string = code_string.replace('VARIANT', '3')


        if "protein_coding" in variants.INFO.get('ANN'):
            snp_type = "Coding SNP"
        else:
            snp_type = "Non-coding SNP"
        print_string = print_string + snp_type + " at %s > " % str(variants.POS) + str(",".join(variants.ALT))

        ann_array = (variants.INFO.get('ANN')).split(',')
        ann_string = ";"
        for i in list(set(ann_array)):
            i_split = i.split('|')
            ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"
            tag = str(i_split[3]).replace('CHR_START-', '')
            if "-" in tag:
                extra_tags = ""
                tag_split = tag.split('-')
                for i in tag_split:
                    extra_tags = extra_tags + locus_tag_to_gene_name[i] + ","
                extra_tags_prot = ""
                for i in tag_split:
                    extra_tags_prot = extra_tags_prot + locus_tag_to_product[i] + ","
                ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags, extra_tags_prot]) + ";"
            else:
                extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"

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

        final_allele_string = print_string + gt_string + '\n'
        final_code_string = print_string + "," + code_string + '\n'
        final_allele_string = final_allele_string.replace(',|', '|')
        final_allele_string = final_allele_string.replace(',;,', ';')
        final_code_string = final_code_string.replace(',|', '|')
        final_code_string = final_code_string.replace(',;,', ';')

        fp_allele.write(final_allele_string)
        fp_code.write(final_code_string)

    fp_code.close()
    fp_allele.close()



    ##Indel
    header_print_string = "Type of SNP at POS > ALT; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos"
    final_merge_anno_file = VCF("%s/Final_vcf_indel.vcf.gz" % args.filter2_only_snp_vcf_dir)
    for sample in final_merge_anno_file.samples:
        header_print_string = header_print_string + "," + sample
    header_print_string = header_print_string + "\n"

    fp_code = open("%s/Indel_matrix_code.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele = open("%s/Indel_matrix_allele.csv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_code.write(header_print_string)
    fp_allele.write(header_print_string)

    for variants in VCF("%s/Final_vcf_indel.vcf.gz" % args.filter2_only_snp_vcf_dir):
        print_string = ""

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
        else:
            code_string = code_string.replace('VARIANT', '3')


        if "protein_coding" in variants.INFO.get('ANN'):
            snp_type = "Coding SNP"
        else:
            snp_type = "Non-coding SNP"
        print_string = print_string + snp_type + " at %s > " % str(variants.POS) + str(",".join(variants.ALT))

        ann_array = (variants.INFO.get('ANN')).split(',')
        ann_string = ";"
        # for i in list(set(ann_array)):
        #     i_split = i.split('|')
        #     ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"
        # print_string = print_string + ann_string
        for i in list(set(ann_array)):
            i_split = i.split('|')
            ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"
            tag = str(i_split[3]).replace('CHR_START-', '')
            if "-" in tag:
                extra_tags = ""
                tag_split = tag.split('-')
                for i in tag_split:
                    extra_tags = extra_tags + locus_tag_to_gene_name[i] + ","
                extra_tags_prot = ""
                for i in tag_split:
                    extra_tags_prot = extra_tags_prot + locus_tag_to_product[i] + ","
                ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags, extra_tags_prot]) + ";"
            else:
                extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"

        print_string = print_string + ann_string

        gt_string = ""
        for gt in variants.gt_bases:
            gt = gt.replace('./.', '.')
            if "/" in gt:
                gt_split = gt.split('/')
                gt = gt_split[0]
            gt_string = gt_string + "," + gt
        gt_string = gt_string.replace('.', variants.REF)
        final_allele_string = print_string + gt_string + '\n'
        final_code_string = print_string + "," + code_string + '\n'
        final_allele_string = final_allele_string.replace(',|', '|')
        final_allele_string = final_allele_string.replace(',;,', ';')
        final_code_string = final_code_string.replace(',|', '|')
        final_code_string = final_code_string.replace(',;,', ';')

        fp_allele.write(final_allele_string)
        fp_code.write(final_code_string)
    fp_code.close()
    fp_allele.close()


def alignment_report(data_matrix_dir):
    print "\nGenerating Alignment report...\n"
    varcall_dir = os.path.dirname(os.path.abspath(args.results_dir))
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
    print "Alignment report can be found in %s/Report_alignment.txt" % data_matrix_dir



def variant_report(data_matrix_dir):
    print "\nGenerating Variants report...\n"
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
    print "Variant call report can be found in %s/Report_variants.txt" % data_matrix_dir


def fasttree(tree_dir, input_fasta, cluster):
    keep_logging('Running Fasttree on input: %s' % input_fasta, 'Running Fasttree on input: %s' % input_fasta, logger, 'info')
    fasttree_cmd = "%s/%s/%s -nt %s > %s/%s_FastTree.tree" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("fasttree", Config)['fasttree_bin'], ConfigSectionMap("fasttree", Config)['base_cmd'], input_fasta, tree_dir, (os.path.basename(input_fasta)).replace('.fa', ''))
    keep_logging('%s' % fasttree_cmd, '%s' % fasttree_cmd, logger, 'info')
    if cluster == "parallel-local" or cluster == "local":
        os.system("cd %s" % tree_dir)
        os.system(fasttree_cmd)
    elif cluster == "parallel-cluster":
        job_name = os.path.basename(tree_dir)
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=76:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], tree_dir, fasttree_cmd)
        job_file_name = "%s/fasttree_%s.pbs" % (tree_dir, os.path.basename(input_fasta))
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        os.system("qsub %s" % job_file_name)

def raxml(tree_dir, input_fasta):
    keep_logging('Running RAXML on input: %s' % input_fasta, 'Running RAXML on input: %s' % input_fasta, logger, 'info')
    raxml_cmd = "%s/%s/%s %s -s %s -n %s_raxML" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("raxml", Config)['raxml_bin'], ConfigSectionMap("raxml", Config)['base_cmd'], ConfigSectionMap("raxml", Config)['parameters'], input_fasta, (os.path.basename(input_fasta)).replace('.fa', ''))
    keep_logging('%s' % raxml_cmd, '%s' % raxml_cmd, logger, 'info')
    if args.jobrun == "parallel-local" or args.jobrun == "local":
        os.system("cd %s" % tree_dir)
        os.system(raxml_cmd)
    elif args.jobrun == "parallel-cluster":
        job_name = os.path.basename(tree_dir)
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=76:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], tree_dir, raxml_cmd)
        job_file_name = "%s/raxml_%s.pbs" % (tree_dir, os.path.basename(input_fasta))
        f1=open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        os.system("qsub %s" % job_file_name)


def gubbins(gubbins_dir, input_fasta):
    print "\nRunning Gubbins on input: %s\n" % input_fasta
    os.system("cd %s" % ConfigSectionMap("gubbins", Config)['gubbins_bin'])
    gubbins_cmd = "%s/%s --prefix %s/%s %s" % (ConfigSectionMap("gubbins", Config)['gubbins_bin'], ConfigSectionMap("gubbins", Config)['base_cmd'], gubbins_dir, (os.path.basename(input_fasta)).replace('.fa', ''), input_fasta)
    print gubbins_cmd
    os.system(gubbins_cmd)

def core_prep_snp(core_vcf_fasta_dir):

    """ Run snpEff annotation step """
    variant_annotation()

    # """ Generate SNP Filter Label Matrix """
    generate_paste_command()
    #
    # """ Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix. """
    generate_position_label_data_matrix()
    #
    # """ Generate VCF files from final list of variants in Only_ref_variant_positions_for_closely; generate commands for consensus generation """
    generate_vcf_files()
    #
    # """ Generate consensus fasta file from core vcf files """
    extract_only_ref_variant_fasta_from_reference()
    #
    # """ Generate consensus fasta file with only reference and variant position bases """
    extract_only_ref_variant_fasta(core_vcf_fasta_dir)
    #
    # """ Analyze the positions that were filtered out only due to insufficient depth"""
    DP_analysis()

def core_prep_indel(core_vcf_fasta_dir):

    """ Run snpEff annotation step """
    indel_annotation()

    # """ Generate SNP Filter Label Matrix """
    generate_indel_paste_command()

    # """ Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix. """
    generate_indel_position_label_data_matrix()



"""
Pending inclusion
"""

class FuncThread(threading.Thread):
    def __init__(self, target, *args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)
    def run(self):
        self._target(*self._args)

def someOtherFunc(data, key):
    print "someOtherFunc was called : data=%s; key=%s" % (str(data), str(key))

def run_phaster(reference_genome):
    print "\nRunning Phaster on input reference genome: %s\n" % reference_genome
    out_name = (os.path.basename(reference_genome)).split('.')
    phaster_post_cmd = "wget --post-file=\"%s\" \"http://phaster.ca/phaster_api\" -O %s/%s" % (reference_genome, args.filter2_only_snp_vcf_dir, str(out_name[0]) + "_phaster_post.json")

    print "Running: %s\n" % phaster_post_cmd
    #os.system(phaster_post_cmd)
    with open('%s/%s' % (args.filter2_only_snp_vcf_dir, str(out_name[0]) + "_phaster_post.json")) as json_data:
        data = json.load(json_data)
        print "Status: %s\njob_id: %s\n" % (data["status"], data["job_id"])

def parse_phaster(reference_genome):
    out_name = (os.path.basename(reference_genome)).split('.')
    with open('%s/%s' % (args.filter2_only_snp_vcf_dir, str(out_name[0]) + "_phaster_post.json")) as json_data:
        data = json.load(json_data)
        phaster_get_cmd = "wget \"http://phaster.ca/phaster_api?acc=%s\" -O %s/%s" % (data["job_id"], args.filter2_only_snp_vcf_dir, str(out_name[0]) + "_phaster_get.json")
        print phaster_get_cmd

    with open('%s/%s' % (args.filter2_only_snp_vcf_dir, str(out_name[0]) + "_phaster_get.json")) as json_get_data:
        get_data = json.load(json_get_data)
        print get_data["zip"]
        phaster_zip_cmd = "wget \"http://%s\" -O %s/%s_phaster_get.zip" % (str(get_data["zip"]), args.filter2_only_snp_vcf_dir, str(out_name[0]))
        phaster_unzip_cmd = "unzip %s/%s_phaster_get.zip" % (args.filter2_only_snp_vcf_dir, str(out_name[0]))
        print phaster_zip_cmd
        print phaster_unzip_cmd
        # for key, value in get_data.items():
        #     print get_data["zip"][0]
"""
Pending inclusion
"""



#Main Steps
if __name__ == '__main__':

    """Start Timer"""
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    global logger
    analysis_name_log = "step_" + str(args.steps)
    logger = generate_logger(args.filter2_only_snp_vcf_dir, analysis_name_log, log_unique_time)
    keep_logging('The Script started at: %s' % start_time, 'The Script started at: %s' % start_time, logger, 'info')
    print_details = "This step will parse final vcf files(*_no_proximate_snp.vcf) generated at the end of Variant Calling Pipeline. At the end of this step, the following results will be generated and placed in output directory:\n\n" \
          "1. Final Core SNP Positions list(Variant positions that were not filtered out in any of the samples and passed all the filters)\n" \
          "2. SNP Positions that were filtered out with labels indicating the reason (Depth, FQ, MQ, Unmapped in one or other samples, Proximate SNPS, Quality of Variant) why they were filtered out.\n" \
          "3. Barplot Statistics about the filtered variants and their reason for getting filtered.\n" \
          "4. Final Consensus fasta file using only Core SNP Positions\n"
    keep_logging('%s' % print_details, '%s' % print_details, logger, 'info')

    """ Create Temp Directory for storing unwanted temp files generated while running script """
    temp_dir = args.filter2_only_snp_vcf_dir + "/temp"
    make_sure_path_exists(temp_dir)

    filter2_only_snp_vcf_filenames = args.filter2_only_snp_vcf_filenames
    vcf_filenames = []
    with open(filter2_only_snp_vcf_filenames) as fp:
        for line in fp:
            line = line.strip()
            line = args.filter2_only_snp_vcf_dir + line
            vcf_filenames.append(line)
        fp.close()
    global config_file
    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
    global Config
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    keep_logging('Path to config file: %s' % config_file, 'Path to config file: %s' % config_file, logger, 'info')



    ### Start the core SNP pipeline steps
    """ core_prep step """
    if "1" in args.steps:
        keep_logging('Gathering SNP position information from each final *_no_proximate_snp.vcf file...', 'Gathering SNP position information from each final *_no_proximate_snp.vcf file...', logger, 'info')

        """
        Gather SNP positions from each final *_no_proximate_snp.vcf file (that passed the variant filter parameters
        from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files for use in downstream methods
        """
        unique_position_file = create_positions_filestep(vcf_filenames)

        unique_indel_position_file = create_indel_positions_filestep(vcf_filenames)

        tmp_dir = "/tmp/temp_%s/" %log_unique_time
        
        bgzip_cmd = "for i in %s/*.vcf; do bgzip -c $i > $i%s; done" % (args.filter2_only_snp_vcf_dir, ".gz")
        tabix_cmd = "for i in %s/*.vcf.gz; do tabix -f $i; done" % (args.filter2_only_snp_vcf_dir)
        os.system(bgzip_cmd)
        os.system(tabix_cmd)

        """ Get the cluster option; create and run jobs based on given parameter """
        create_job(args.jobrun, vcf_filenames, unique_position_file, tmp_dir)

        create_indel_job(args.jobrun, vcf_filenames, unique_indel_position_file, tmp_dir)
        """ Find ProPhage region in reference genome """
        #run_phaster(args.reference)


    """ core step """
    if "2" in args.steps:
        # #Adhoc
        data_matrix_dir = args.results_dir + '/data_matrix'
        core_vcf_fasta_dir = args.results_dir + '/core_snp_consensus'
        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(core_vcf_fasta_dir)

        core_prep_snp(core_vcf_fasta_dir)

        core_prep_indel(core_vcf_fasta_dir)

        annotated_snp_matrix()

        keep_logging('Wait for individual cluster jobs to finish before running the third step', 'Wait for individual cluster jobs to finish before running the third step', logger, 'info')

    """ report step """
    if "3" in args.steps:
        keep_logging('Step 3: Generate Reports and Results folder.', 'Step 3: Generate Reports and Results folder.', logger, 'info')
        """ Generate DP barplots data """
        DP_analysis_barplot()

        """ Analyze the FQ values of all the unique variant """
        FQ_analysis()

        data_matrix_dir = args.results_dir + '/data_matrix'
        core_vcf_fasta_dir = args.results_dir + '/core_snp_consensus'
        consensus_var_dir = core_vcf_fasta_dir + '/consensus_variant_positions'
        consensus_ref_var_dir = core_vcf_fasta_dir + '/consensus_ref_variant_positions'

        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(core_vcf_fasta_dir)
        make_sure_path_exists(consensus_var_dir)
        make_sure_path_exists(consensus_ref_var_dir)

        move_data_matrix_results = "cp -r %s/*.txt %s/temp* %s/All* %s/Only* %s/*.R %s/R_scripts/generate_diagnostics_plots.R %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, os.path.dirname(os.path.abspath(__file__)), data_matrix_dir)
        #move_core_vcf_fasta_results = "cp %s/*_core.vcf.gz %s/*.fa %s/*_variants.fa %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, core_vcf_fasta_dir)
        move_core_vcf_fasta_results = "cp %s/*_core.vcf.gz %s/*.fa %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, core_vcf_fasta_dir)
        move_consensus_var_fasta_results = "mv %s/*_variants.fa %s/" % (core_vcf_fasta_dir, consensus_var_dir)
        move_consensus_ref_var_fasta_results = "mv %s/*.fa %s/" % (core_vcf_fasta_dir, consensus_ref_var_dir)


        os.system(move_data_matrix_results)
        os.system(move_core_vcf_fasta_results)
        os.system(move_consensus_var_fasta_results)
        os.system(move_consensus_ref_var_fasta_results)

        subprocess.call(["sed -i 's/title_here/%s/g' %s/generate_diagnostics_plots.R" % (os.path.basename(args.results_dir), data_matrix_dir)], shell=True)

        # Check if the variant consensus files generated are of same length
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
                    print "Error generating variant consensus position file: %s\n" % f
                    keep_logging('Error generating variant consensus position file: %s' % f, 'Error generating variant consensus position file: %s' % f, logger, 'exception')

        """ Generate alignment report """
        alignment_report(data_matrix_dir)

        """ Generate core snps report """
        variant_report(data_matrix_dir)

        print_details = "Results for core pipeline can be found in: %s\n" \
              "Description of Results:\n" \
              "1. data_matrix folder contains all the data matrices and other temporary files generated during the core pipeline. bargraph_counts.txt and bargraph_percentage.txt: contains counts/percentage of unique positions filtered out due to different filter parameters for each sample. Run bargraph.R to plot bargraph statistics." \
              "2. core_snp_consensus contains all the core vcf and fasta files. *_core.vcf.gz: core vcf files, *.fa and *_variants.fa: core consensus fasta file and core consensus fasta with only variant positions." % (args.results_dir)
        keep_logging(print_details, print_details, logger, 'info')

    """ tree step """
    if "4" in args.steps:
        keep_logging('Step 4: Ongoing Testing.', 'Step 4: Ongoing Testing.', logger, 'info')
        #parse_phaster(args.reference)

        gubbins_dir = args.results_dir + '/gubbins'
        tree_dir = args.results_dir + '/trees'

        make_sure_path_exists(gubbins_dir)
        make_sure_path_exists(tree_dir)

        prepare_ref_var_consensus_input = "%s/gubbins/%s_ref_var_consensus.fa" % (args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''))
        prepare_var_consensus_input = "%s/gubbins/%s_var_consensus.fa" % (args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''))

        prepare_ref_var_consensus_input_cmd = "cat %s %s/core_snp_consensus/consensus_ref_variant_positions/*.fa > %s" % (args.reference, args.results_dir, prepare_ref_var_consensus_input)
        prepare_var_consensus_input_cmd = "cat %s/core_snp_consensus/consensus_variant_positions/*.fa > %s" % (args.results_dir, prepare_var_consensus_input)


        os.system(prepare_ref_var_consensus_input_cmd)
        os.system(prepare_var_consensus_input_cmd)

        fasttree(tree_dir, prepare_ref_var_consensus_input, args.jobrun)
        fasttree(tree_dir, prepare_var_consensus_input, args.jobrun)

        raxml(tree_dir, prepare_ref_var_consensus_input)
        raxml(tree_dir, prepare_var_consensus_input)

        # Disabling Gubbins function due to installation issues
        #gubbins(gubbins_dir, prepare_ref_var_consensus_input)

    time_taken = datetime.now() - start_time_2
    if args.remove_temp:
        del_command = "rm -r %s" % temp_dir
        os.system(del_command)















