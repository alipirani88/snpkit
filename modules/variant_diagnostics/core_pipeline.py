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
import ConfigParser
from config_settings import ConfigSectionMap
#from logging_subprocess import *
#from log_modules import *

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


def create_job(jobrun, vcf_filenames):

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
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/reason_job.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i)
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
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/reason_job.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i)
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
            job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/reason_job.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, i)
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir)
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
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            print "Running: qsub %s" % i
            #os.system("qsub %s" % i)


    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')


        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
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
        command_file = "%s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
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
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
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
        aln_mpileup_vcf_file = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf_5bp_indel_removed.vcf')
        analysis = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')
        grep_reference_file = "grep \'^##reference\' %s" % aln_mpileup_vcf_file
        proc = subprocess.Popen([grep_reference_file], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        reference_file = out.split(':')
        gatk_filter2_final_vcf_file = gatk_filter2(aln_mpileup_vcf_file, temp_dir, analysis, reference_file[1])
        gatk_filter2_final_vcf_file_no_proximate_snp = remove_proximate_snps(gatk_filter2_final_vcf_file, temp_dir, analysis, reference_file[1])
        grep_fq_field = "awk -F\'\t\' \'{print $8}\' %s | grep -o \'FQ=.*\' | sed \'s/FQ=//g\' | awk -F\';\' \'{print $1}\' > %s_FQ_values" % (gatk_filter2_final_vcf_file_no_proximate_snp, analysis)
        os.system(grep_fq_field)

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

    """
    Test this again
    """

    # print "Analyzing positions that were filtered out due to Depth..."
    # extract_DP_positions = "awk -F\'\\t\' \'{print $1}\' temp_Only_filtered_positions_for_closely_matrix_DP.txt | sed \'/^$/d\' > extract_DP_positions.txt"
    # os.system(extract_DP_positions)
    #
    # for i in vcf_filenames:
    #     filename_base = os.path.basename(i)
    #     aln_mpileup_vcf_file = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf_5bp_indel_removed.vcf')
    #     analysis = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')
    #     grep_reference_file = "grep \'^##reference\' %s" % aln_mpileup_vcf_file
    #     proc = subprocess.Popen([grep_reference_file], stdout=subprocess.PIPE, shell=True)
    #     (out, err) = proc.communicate()
    #     out = out.strip()
    #     reference_file = out.split(':')
    #     #gatk_filter2_final_vcf_file = gatk_filter2(aln_mpileup_vcf_file, temp_dir, analysis, reference_file[1])
    #     #gatk_filter2_final_vcf_file_no_proximate_snp = remove_proximate_snps(gatk_filter2_final_vcf_file, temp_dir, analysis, reference_file[1])
    #     DP_values_file = "%s/%s_DP_values" % (args.filter2_only_snp_vcf_dir, analysis)
    #     f2=open(DP_values_file, 'w+')
    #
    #
    #     with open("%s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_filess:
    #         csv_readerr = csv.reader(csv_filess, delimiter='\t')
    #         next(csv_readerr, None)
    #         for rows in csv_readerr:
    #             #print rows
    #             #grep_dp_field = "grep -wP \'^\S+\s+%s\s+\b\' %s | awk -F\'\\t\' \'{print $8}\' | grep -o \'DP=.*\' | sed \'s/DP=//g\' | awk -F\';\' \'{print $1}\'" % (rows[0], aln_mpileup_vcf_file)
    #             grep_dp_field = "grep -w \'%s\' %s" % (rows[0], aln_mpileup_vcf_file)
    #             awk_dp_field = "awk -F\'\t\' \'$2 == %s\' %s | awk -F\'\t\' \'{print $8}\' | awk -F\';\' \'{print $1}\' | sed \'s/DP=//g\'" % (rows[0], aln_mpileup_vcf_file)
    #             #print grep_dp_field
    #             #proc = subprocess.Popen([grep_dp_field], stdout=subprocess.PIPE, shell=True)
    #             #(out2, err2) = proc.communicate()
    #             #out_split = out.split('\n')
    #             #out = out.strip()
    #             proc = subprocess.Popen([awk_dp_field], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #             (out2, err2) = proc.communicate()
    #             #print out2.strip()
    #             if out2:
    #                 #print out2.strip()
    #                 if "INDEL" in out2:
    #                     #print awk_dp_field
    #                     out2 == "NA"
    #                 f2.write(out2.strip() + '\n')
    #                 # if len(out_split) > 1:
    #                 #     print out_split[0]
    #                 # # for i in out:
    #                 # #     print i
    #                 # line_split = out.split('\t')
    #                 # #print line_split
    #                 # if line_split[1] == rows[0]:
    #                 #     DP_field = line_split[7].split(';')
    #                 #     DP_value = DP_field[0].replace('DP=', '')
    #                 #print out
    #             else:
    #                 f2.write("NA\n")
    #                 #print "NA"
    #
    # paste_command = "paste %s/extract_DP_positions.txt" % args.filter2_only_snp_vcf_dir
    # for i in vcf_filenames:
    #     label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_DP_values')
    #     paste_command = paste_command + " " + label_file
    #
    # paste_file = args.filter2_only_snp_vcf_dir + "/paste_DP_files.sh"
    # f2=open(paste_file, 'w+')
    # paste_command = paste_command + " > %s/filtered_DP_values_temp.txt" % args.filter2_only_snp_vcf_dir
    # #os.system(paste_command)
    # f2.write(paste_command + '\n')
    # cat_header = "cat %s/header.txt %s/filtered_DP_values_temp.txt > %s/filtered_DP_values.txt" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    # #os.system(cat_header)
    # f2.write(cat_header + '\n')
    # sed_command = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/filtered_DP_values.txt" % (args.filter2_only_snp_vcf_dir)
    # #os.system(sed_command)
    # f2.write(sed_command + '\n')
    # cmd = "bash %s" % paste_file
    # os.system("bash %s/paste_DP_files.sh" % args.filter2_only_snp_vcf_dir)
    #
    # #os.system(cmd) change
    # #subprocess.call(["%s" % cmd], shell=True)
    # #subprocess.check_call('%s' % cmd)

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
    final_fasta_string = ">%s\n" % os.path.basename(args.reference.replace('.fasta', '')) + fasta_string
    fp = open("%s/%s_variants.fa" % (args.filter2_only_snp_vcf_dir, os.path.basename(args.reference.replace('.fasta', ''))), 'w+')
    fp.write(final_fasta_string)
    fp.close()

# def make_sure_path_exists(out_path):
#     try:
#         os.makedirs(out_path)
#     except OSError as exception:
#         if exception.errno != errno.EEXIST:
#             print "Errors in output folder path! please change the output path or analysis name."
#             #keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
#             exit()

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

    print "\nThe Script started at: %s\n" % start_time

    print "\nThis step will parse final vcf files(*_no_proximate_snp.vcf) generated at the end of Variant Calling Pipeline. At the end of this step, the following results will be generated and placed in output directory:\n\n" \
          "1. Final Core SNP Positions list(Variant positions that were not filtered out in any of the samples and passed all the filters)\n" \
          "2. SNP Positions that were filtered out with labels indicating the reason (Depth, FQ, MQ, Unmapped in one or other samples, Proximate SNPS, Quality of Variant) why they were filtered out.\n" \
          "3. Barplot Statistics about the filtered variants and their reason for getting filtered.\n" \
          "4. Final Consensus fasta file using only Core SNP Positions\n"

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
    print config_file
    if "1" in args.steps:
        print "Gathering SNP position information from each final *_no_proximate_snp.vcf file..."
        """
        Gather SNP positions from each final *_no_proximate_snp.vcf file (that passed the variant filter parameters
        from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files for use in downstream methods
        """
        create_positions_filestep(vcf_filenames)

        """ Get the cluster option; create and run jobs based on given parameter """
        create_job(args.jobrun, vcf_filenames)

        """ Find ProPhage region in reference genome """
        #run_phaster(args.reference)

    if "2" in args.steps:
        #Adhoc
        data_matrix_dir = args.results_dir + '/data_matrix'
        core_vcf_fasta_dir = args.results_dir + '/core_snp_consensus'
        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(core_vcf_fasta_dir)

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


        print "Wait for individual cluster jobs to finish before running the third step"

    if "3" in args.steps:
        print "\nStep 3: Generate Reports and Results folder.\n"

        data_matrix_dir = args.results_dir + '/data_matrix'
        core_vcf_fasta_dir = args.results_dir + '/core_snp_consensus'

        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(core_vcf_fasta_dir)

        move_data_matrix_results = "cp -r %s/*.txt %s/temp* %s/All* %s/Only* %s/*.R %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, data_matrix_dir)
        move_core_vcf_fasta_results = "cp %s/*_core.vcf.gz %s/*.fa %s/*_variants.fa %s/" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, core_vcf_fasta_dir)

        os.system(move_data_matrix_results)
        os.system(move_core_vcf_fasta_results)

        # Check if the variant consensus files generated are of same length
        count = 0
        for line in open("%s/Only_ref_variant_positions_for_closely_matrix.txt" % data_matrix_dir).xreadlines(  ): count += 1
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


        """ Generate DP barplots data """
        DP_analysis_barplot()

        """ Analyze the FQ values of all the unique variant """
        #FQ_analysis()

        print "\nResults for core pipeline can be found in: %s\n" \
              "\nDescription of Results:\n" \
              "\n1. data_matrix folder contains all the data matrices and other temporary files generated during the core pipeline. bargraph_counts.txt and bargraph_percentage.txt: contains counts/percentage of unique positions filtered out due to different filter parameters for each sample. Run bargraph.R to plot bargraph statistics." \
              "\n2. core_snp_consensus contains all the core vcf and fasta files. *_core.vcf.gz: core vcf files, *.fa and *_variants.fa: core consensus fasta file and core consensus fasta with only variant positions." % (args.results_dir)


    if "4" in args.steps:
        print "Step 4: Require Testing."
        #parse_phaster(args.reference)

    time_taken = datetime.now() - start_time_2
    if args.remove_temp:
        del_command = "rm -r %s" % temp_dir
        os.system(del_command)















