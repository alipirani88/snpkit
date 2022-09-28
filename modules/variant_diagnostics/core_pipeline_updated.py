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
import pandas as pd
import numpy as np
import errno
from pyfasta import Fasta
from datetime import datetime
import time
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
from core_prep_functions import *
from iqtree import iqtree
from memory_profiler import profile

# Parse Command line Arguments
parser = argparse.ArgumentParser(
    description='Parsing filtered VCF files and investigating Variants to determine the reason why it was filtered out from the final list')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                      help='Directory where all the filter2 only SNP vcf files are saved.')
required.add_argument('-filter2_only_snp_vcf_filenames', action='store', dest="filter2_only_snp_vcf_filenames",
                      help='Names of filter2 only SNP vcf files with name per line.')
optional.add_argument('-jobrun', action='store', dest="jobrun",
                      help='Running a job on Cluster, Running Parallel jobs, Run jobs/commands locally (default): cluster, local, parallel-local, parallel-single-cluster')
optional.add_argument('-scheduler', action='store', dest="scheduler",
                      help='Type of Cluster: PBS, SLURM')
optional.add_argument('-cluster_resources', action='store', dest="cluster_resources",
                      help='Cluster Resources to use. for example nodes,core. Ex: 1,4')
optional.add_argument('-numcores', action='store', dest="numcores",
                      help='Number of cores to use on local system for parallel-local parameter')
optional.add_argument('-remove_temp', action='store', dest="remove_temp",
                      help='Remove Temporary files generated during the run')
optional.add_argument('-gubbins', action='store', dest="gubbins", help='yes/no for running gubbins')
optional.add_argument('-outgroup', action='store', dest="outgroup", help='outgroup sample name')
required.add_argument('-reference', action='store', dest="reference",
                      help='Path to Reference Fasta file for consensus generation')
optional.add_argument('-mask', action='store_true', dest="mask",
                      help='Mask Gubbins detected recombinant region in WGA and run Iqtree on masked alignment')
optional.add_argument('-readme', action='store', dest="readme",
                      help='Generate a README file describing this run. Fill up the form /variant_calling_pipeline/readme_metadata_form.txt residing under the code directory or provide it with this argument. DEFAULT - /variant_calling_pipeline/readme_metadata_form.txt')
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

""" Generic Methods """


def make_sure_path_exists(out_path):
    """This function checks if the args out_path exists and generates an empty directory if it doesn't.

    :param:
        out_path: Directory path to check or create a new directory.

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

def get_outgroup():
    """
    Prepare Outgroup Sample name from the argument.
    """
    if args.outgroup:
        if "R1_001_final.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        elif "_R1.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        elif "R1.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
            outgroup = re.sub("_S.*", "", first_part)

        elif "1_combine.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        elif "1_sequence.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        elif "_forward.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        elif "R1_001.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        elif "_1.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        elif ".1.fastq.gz" in args.outgroup:
            first_part_split = args.outgroup.split('.1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            outgroup = re.sub("_S.*_", "", first_part)

        keep_logging(
            'Using %s as Outgroup Sample Name' % outgroup,
            'Using %s as Outgroup Sample Name' % outgroup,
            logger, 'info')

        return outgroup
    else:
        keep_logging('- Note: Outgroup Sample Name not provided. Ignoring Outgroup steps.', '- Note: Outgroup Sample Name not provided. Ignoring Outgroup steps.', logger, 'info')
        outgroup = ""


## Great Lakes Integration Changes
def get_scheduler_directive(scheduler, Config):
    """ Generate Cluster Directive lines for a scheduler provided with args.scheduler"""
    # Scheduler Changes here; current changes
    if scheduler and scheduler == "SLURM":
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                               % (ConfigSectionMap("slurm", Config)['email'],
                                  ConfigSectionMap("slurm", Config)['notification'],
                                  ConfigSectionMap("slurm", Config)['partition'],
                                  ConfigSectionMap("slurm", Config)['flux_account'],
                                  ConfigSectionMap("slurm", Config)['resources'])
    elif scheduler and scheduler == "PBS":
        script_Directive = "#PBS"
        job_name_flag = "-N"
        scheduler_directives = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                               % (ConfigSectionMap("scheduler", Config)['email'],
                                  ConfigSectionMap("scheduler", Config)['notification'],
                                  ConfigSectionMap("scheduler", Config)['resources'],
                                  ConfigSectionMap("scheduler", Config)['queue'],
                                  ConfigSectionMap("scheduler", Config)['flux_account'])
    else:
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                               % (ConfigSectionMap("slurm", Config)['email'],
                                  ConfigSectionMap("slurm", Config)['notification'],
                                  ConfigSectionMap("slurm", Config)['partition'],
                                  ConfigSectionMap("slurm", Config)['flux_account'],
                                  ConfigSectionMap("slurm", Config)['resources'])
    return scheduler_directives, script_Directive, job_name_flag

def run_command(i):
    """Function to run each command and is run as a part of python Parallel mutiprocessing method.

    :param:
        i: command variable to run

    :return:
        done: string variable with completion status of command.
    """
    print "Running: %s" % i
    call("%s" % i, logger)
    # A subprocess exception is raised if the command finish abnormally.
    # An exception is raised in call method.
    # If none of the exceptions are raised, return done status.
    done = "Completed: %s" % i
    return done


"""core_prep methods

    This block contains methods that are respnsible for running the first part of core_All step of the pipeline.
    This methods generates all the necessary intermediate files required for the second part of core_All step.
    Example of intermediate files: various diagnostics files/matrices where it decides why a variant was filtered out.

    Updates - 2020-02-05 Moving core_prep methods to core_prep_functions python library
"""

"""core methods 

    This block contains methods that are respnsible for running the second part of core_All step of the pipeline.
    It uses intermediate files generated during the first step, finds core SNPs and annotates variants using snpEff.
    It will generate all types of SNP matrices that is required for downstream pathways / Association analysis.
    Output:
        - 

"""


def generate_paste_command():
    """
    This Function will take all the *label file and generate/paste it column wise to generate a matrix. These matrix will be used in downstream analysis.
    :param: null
    :return: null
    """

    """ Paste/Generate and sort SNP Filter Label Matrix """
    # Initialize Paste command. First column should be unique positions. 
    paste_command = "paste %s/unique_positions_file" % args.filter2_only_snp_vcf_dir
    
    # Add VCF label filenames to the Paste commands.
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                               '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
        paste_command = paste_command + " " + label_file
    
    # Generate a header line. First column should be a tab.
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s | sed \'s/[[:space:]]$//g\' > %s/header.txt" % (
    args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    call("%s" % header_awk_cmd, logger)
    call("%s" % sed_header, logger)
    call("%s" % sed_header_2, logger)

    
    paste_command = paste_command + " > %s/temp_label_final_raw.txt" % args.filter2_only_snp_vcf_dir
    call(paste_command, logger)
    
    paste_command_header = "cat %s/header.txt %s/temp_label_final_raw.txt > %s/All_label_final_sorted_header.txt" % (
    args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    
    call("%s" % paste_command_header, logger)

    """ Assign numeric code to each variant filter reason"""
    subprocess.call([
                        "sed -i 's/reference_unmapped_position/0/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
    subprocess.call(
        [
            "sed -i 's/reference_allele_indel_proximate/1/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/reference_allele/1/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(["sed -i 's/VARIANT/1TRUE/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/LowFQ_QUAL_DP_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/LowFQ_DP_QUAL_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/LowFQ_QUAL_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(
        ["sed -i 's/LowFQ_DP_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/LowFQ_proximate_SNP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/LowFQ_QUAL_DP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/LowFQ_DP_QUAL/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(["sed -i 's/LowFQ_QUAL/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(["sed -i 's/LowFQ_DP/2/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(
        ["sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/HighFQ_proximate_SNP/7/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/HighFQ_QUAL_DP/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/HighFQ_DP_QUAL/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(["sed -i 's/HighFQ_QUAL/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(["sed -i 's/HighFQ_DP/3/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(["sed -i 's/LowFQ/5/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(["sed -i 's/HighFQ/6/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    remove_unwanted_text = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir
    call("%s" % remove_unwanted_text, logger)

def generate_paste_command_outgroup():
    """
    This Function will take all the *label file and generate/paste it column wise to generate a matrix. These matrix will be used in downstream analysis.
    :param: null
    :return: null
    """

    if args.outgroup:
        """ Paste/Generate and sort SNP Filter Label Matrix """
        paste_file = args.filter2_only_snp_vcf_dir + "/temp/paste_label_files_outgroup.sh"
        f4 = open(paste_file, 'w+')
        paste_command = "paste %s/unique_positions_file" % args.filter2_only_snp_vcf_dir
        for i in vcf_filenames:
            if "%s_filter2_final.vcf_no_proximate_snp.vcf" % outgroup not in i:
                label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                       '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
                paste_command = paste_command + " " + label_file

        """Exclude outgroup sample name in header

        header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
        sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
        sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

        """

        header_awk_cmd = "grep -v \'%s\' %s | awk \'{ORS=\"\t\";}{print $1}\' > %s/header_outgroup.txt" % (
        outgroup, args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
        sed_header = "sed -i \'s/^/\t/\' %s/header_outgroup.txt" % args.filter2_only_snp_vcf_dir
        sed_header_2 = "sed -i -e \'$a\\' %s/header_outgroup.txt" % args.filter2_only_snp_vcf_dir

        call("%s" % header_awk_cmd, logger)
        call("%s" % sed_header, logger)
        call("%s" % sed_header_2, logger)

        paste_command = paste_command + " > %s/temp_label_final_raw_outgroup.txt" % args.filter2_only_snp_vcf_dir
        f4.write(paste_command)
        f4.close()
        paste_command_header = "cat %s/header_outgroup.txt %s/temp_label_final_raw_outgroup.txt > %s/All_label_final_sorted_header_outgroup.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
        call("%s" % paste_command_header, logger)

        """ Assign numeric code to each variant filter reason"""
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
        subprocess.call([
            "sed -i 's/reference_allele_indel_proximate/1/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
            shell=True)
        subprocess.call([
                            "sed -i 's/reference_allele/1/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/VARIANT/1TRUE/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_DP_proximate_SNP/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_QUAL_proximate_SNP/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_proximate_SNP/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_proximate_SNP/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_proximate_SNP/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_DP/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_QUAL/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call(
            ["sed -i 's/LowFQ_QUAL/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
            shell=True)
        subprocess.call(
            ["sed -i 's/LowFQ_DP/2/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
            shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_proximate_SNP/7/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_DP/3/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_QUAL/3/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL/3/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call(
            ["sed -i 's/HighFQ_DP/3/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
            shell=True)
        subprocess.call(
            ["sed -i 's/LowFQ/5/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
            shell=True)
        subprocess.call(
            ["sed -i 's/HighFQ/6/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
            shell=True)
        remove_unwanted_text = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir
        call("%s" % remove_unwanted_text, logger)

def generate_indel_paste_command():
    """
    This Function will take all the *label file and generate/paste it column wise to generate a matrix. These matrix will be used in downstream analysis.
    :param: null
    :return: null
    """

    """ Paste/Generate and sort SNP Filter Label Matrix """
    paste_command = "paste %s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                               '_filter2_indel_final.vcf_indel_positions_label')
        paste_command = paste_command + " " + label_file
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s | sed \'s/[[:space:]]$//g\' > %s/header.txt" % (
    args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

    call("%s" % header_awk_cmd, logger)
    call("%s" % sed_header, logger)
    call("%s" % sed_header_2, logger)

    
    paste_command = paste_command + " > %s/temp_indel_label_final_raw.txt" % args.filter2_only_snp_vcf_dir

    call(paste_command, logger)

    paste_command_header = "cat %s/header.txt %s/temp_indel_label_final_raw.txt > %s/All_indel_label_final_sorted_header.txt" % (
    args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    call("%s" % paste_command_header, logger)

    """ Assign numeric code to each variant filter reason"""
    subprocess.call([
                        "sed -i 's/reference_unmapped_position/0/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
    subprocess.call(
        [
            "sed -i 's/reference_allele_indel_proximate/1/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/reference_allele/1/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/VARIANT/1TRUE/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call([
                        "sed -i 's/LowAF_QUAL_DP_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/LowAF_DP_QUAL_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/LowAF_QUAL_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/LowAF_DP_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/LowAF_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(
        ["sed -i 's/LowAF_QUAL_DP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/LowAF_DP_QUAL/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/LowAF_QUAL/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/LowAF_DP/2/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call([
                        "sed -i 's/HighAF_QUAL_DP_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/HighAF_DP_QUAL_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/HighAF_QUAL_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/HighAF_DP_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call([
                        "sed -i 's/HighAF_proximate_SNP/7/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(
        ["sed -i 's/HighAF_QUAL_DP/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/HighAF_DP_QUAL/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/HighAF_QUAL/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(
        ["sed -i 's/HighAF_DP/3/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    subprocess.call(["sed -i 's/LowAF/5/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
                    shell=True)
    subprocess.call(
        ["sed -i 's/HighAF/6/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
        shell=True)
    remove_unwanted_text = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir
    call("%s" % remove_unwanted_text, logger)

def generate_indel_paste_command_outgroup():
    """
    This Function will take all the *label file and generate/paste it column wise to generate a matrix. These matrix will be used in downstream analysis.
    :param: null
    :return: null
    """

    if args.outgroup:
        """ Paste/Generate and sort SNP Filter Label Matrix """
        # define a file name where the paste commands will be saved.
        paste_file = args.filter2_only_snp_vcf_dir + "/temp/paste_indel_label_files_outgroup.sh"
        f4 = open(paste_file, 'w+')

        # initiate paste command string
        paste_command = "paste %s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir

        # Generate paste command
        for i in vcf_filenames:
            if "%s_filter2_final.vcf_no_proximate_snp.vcf" % outgroup not in i:
                label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                       '_filter2_indel_final.vcf_indel_positions_label')
                paste_command = paste_command + " " + label_file
        # Change header awk command to exclude outgroup
        # header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
        header_awk_cmd = "grep -v \'%s\' %s | awk \'{ORS=\"\t\";}{print $1}\' > %s/header_outgroup.txt" % (
        outgroup, args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
        sed_header = "sed -i \'s/^/\t/\' %s/header_outgroup.txt" % args.filter2_only_snp_vcf_dir
        sed_header_2 = "sed -i -e \'$a\\' %s/header_outgroup.txt" % args.filter2_only_snp_vcf_dir

        call("%s" % header_awk_cmd, logger)
        call("%s" % sed_header, logger)
        call("%s" % sed_header_2, logger)

        
        paste_command = paste_command + " > %s/temp_indel_label_final_raw_outgroup.txt" % args.filter2_only_snp_vcf_dir
        f4.write(paste_command)
        f4.close()

        call(paste_command, logger)
    
        paste_command_header = "cat %s/header_outgroup.txt %s/temp_indel_label_final_raw_outgroup.txt > %s/All_indel_label_final_sorted_header_outgroup.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

        call("%s" % paste_command_header, logger)
        
        """ Assign numeric code to each variant filter reason"""
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
        subprocess.call([
            "sed -i 's/reference_allele_indel_proximate/1/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
            shell=True)
        subprocess.call([
                            "sed -i 's/reference_allele/1/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/VARIANT/1TRUE/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_DP_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_QUAL_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_proximate_SNP/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_DP/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_QUAL/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP/2/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_DP_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_QUAL_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_proximate_SNP/4/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_proximate_SNP/7/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_DP/3/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_QUAL/3/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL/3/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP/3/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF/5/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF/6/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        remove_unwanted_text = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir
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
    method_start_time = datetime.now()

    def generate_position_label_data_matrix_All_label():
        position_label = OrderedDict()
        f1 = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'w+')
        other_codes = ['0', '2', '3', '4', '5', '6', '7']
        ref_var = ['1', '1TRUE']
        
        if args.outgroup:
            Only_ref_variant_positions_for_closely_outgroup = []
            Only_filtered_variant_positions_for_closely_outgroup = []

            print_string_header = "\t"
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"

            with open("%s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir,
                      'rU') as csv_file:
                keep_logging(
                    '- Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir,
                    '- Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir,
                    logger, 'info')
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    position_label[row[0]] = row[1:]
                csv_file.close()
                
                for value in position_label:
                    if set(ref_var) & set(position_label[value]):
                        if set(other_codes) & set(position_label[value]):
                            if int(value) not in outgroup_specific_positions:
                                print_string = ""
                                for i in position_label[value]:
                                    print_string = print_string + "\t" + i
                                STRR2 = value + print_string + "\n"
                                Only_filtered_variant_positions_for_closely_outgroup.append(value)
                        else:
                            if int(value) not in outgroup_specific_positions:
                                strr = value + "\n"
                                f1.write(strr)
                                if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
                                    if line not in functional_filter_pos_array:
                                        Only_ref_variant_positions_for_closely_outgroup.append(value)
                                else:
                                    Only_ref_variant_positions_for_closely_outgroup.append(value)
            f1.close()
            

        else:
            All_label_final_sorted_header = pd.read_csv("%s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, sep='\t', header=0)
            All_label_final_sorted_header.replace('1TRUE','1', inplace=True, regex=True)
            Code_count = All_label_final_sorted_header.iloc[:, 1:].apply(pd.Series.value_counts, axis=1)
            numberofsamples = len(All_label_final_sorted_header.columns) - 1
            frames = [All_label_final_sorted_header['Unnamed: 0'], Code_count]
            results = pd.concat(frames, axis=1, join='inner')
            Only_ref_variant_positions_for_closely = results.loc[results['1'] == numberofsamples, 'Unnamed: 0']
            Only_filtered_variant_positions_for_closely = results.loc[results['1'] != numberofsamples, 'Unnamed: 0']
            Only_ref_variant_positions_for_closely.to_csv('Only_ref_variant_positions_for_closely', index=False, sep='\n')
            Only_filtered_variant_positions_for_closely.to_csv('Only_filtered_positions_for_closely', index=False, sep='\n')

        return Only_ref_variant_positions_for_closely

    def temp_generate_position_label_data_matrix_All_label():

        """
        Read temp_label_final_raw.txt SNP position label data matrix for generating barplot statistics.
        """
        temp_position_label = OrderedDict()
        f33 = open("%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        print_string_header = "\t"

        if args.outgroup:
            for i in vcf_filenames:
                if "%s_filter2_final.vcf_no_proximate_snp.vcf" % outgroup not in i:
                    print_string_header = print_string_header + os.path.basename(i) + "\t"
        else:
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"

        f33.write('\t' + print_string_header.strip() + '\n')

        """ GET individual PHAGE/Repetitive/masked region positions to assign functional class group string """

        phage_positions = []

        phage_region_positions = "%s/phage_region_positions.txt" % args.filter2_only_snp_vcf_dir
        if os.path.isfile(phage_region_positions):
            with open(phage_region_positions, 'rU') as fphage:
                for line in fphage:
                    phage_positions.append(line.strip())
            fphage.close()
        else:
            raise IOError('%s/phage_region_positions.txt does not exist.' % args.filter2_only_snp_vcf_dir)
            exit()

        """ End: Generate a list of functional class positions from Phaster, Mummer and Custom Masking results/files"""

        f_open_temp_Only_filtered_positions_for_closely_matrix = open(
            "%s/temp_Only_filtered_positions_for_closely_matrix_exclude_phage.txt" % args.filter2_only_snp_vcf_dir,
            'w+')

        f_open_temp_Only_filtered_positions_for_closely_matrix.write('\t' + print_string_header.strip() + '\n')

        keep_logging(
            'Reading temporary label positions file: %s/temp_label_final_raw.txt \n' % args.filter2_only_snp_vcf_dir,
            'Reading temporary label positions file: %s/temp_label_final_raw.txt \n' % args.filter2_only_snp_vcf_dir,
            logger, 'info')
        other_codes = ['reference_unmapped_position', 'LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP',
               'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ', 'LowFQ_proximate_SNP',
               'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP',
               'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP',
               'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
        # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
        ref_var = ['reference_allele', 'reference_allele_indel_proximate', 'VARIANT']

        if args.outgroup:
            with open("%s/temp_label_final_raw_outgroup.txt" % args.filter2_only_snp_vcf_dir, 'r') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    if set(ref_var) & set(row[1:]):
                        if set(other_codes) & set(row[1:]):
                            if int(row[0]) not in outgroup_specific_positions:

                                print_string = ""
                                for i in row[1:]:
                                    print_string = print_string + "\t" + i
                                STRR2 = row[0] + print_string + "\n"
                                f33.write(STRR2)

                                if str(row[0]) not in phage_positions:
                                    print_string_2 = ""
                                    for i in row[1:]:
                                        print_string_2 = print_string_2 + "\t" + i
                                    STRR3 = row[0] + print_string_2 + "\n"
                                    f_open_temp_Only_filtered_positions_for_closely_matrix.write(STRR3)
            csv_file.close()
            f33.close()
            f_open_temp_Only_filtered_positions_for_closely_matrix.close()

        else:
            program_starts = time.time()
            # print "Debugging this method."
            with open("%s/temp_label_final_raw.txt" % args.filter2_only_snp_vcf_dir, 'r') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    if set(ref_var) & set(row[1:]):
                        if set(other_codes) & set(row[1:]):
                            print_string = ""
                            for i in row[1:]:
                                print_string = print_string + "\t" + i
                            STRR2 = row[0] + print_string + "\n"
                            f33.write(STRR2)
                            if str(row[0]) not in phage_positions:
                                print_string_2 = ""
                                for i in row[1:]:
                                    print_string_2 = print_string_2 + "\t" + i
                                STRR3 = row[0] + print_string_2 + "\n"
                                f_open_temp_Only_filtered_positions_for_closely_matrix.write(STRR3)

            csv_file.close()
            f33.close()
            f_open_temp_Only_filtered_positions_for_closely_matrix.close()
            now = time.time()
            print "Time taken to iterate the loop once - {0} seconds".format(now - program_starts)
        

    def barplot_stats():
        keep_logging(
            '\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n',
            '\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n',
            logger, 'info')
        """
        Read each Sample columns and calculate the percentage of each label to generate barplot statistics.
        This will give a visual explanation of how many positions in each samples were filtered out because of different reason
        """

        print "Exluding Phage regions from temp_Only_filtered_positions_for_closely_matrix.txt file. The results will be outputed to temp_Only_filtered_positions_for_closely_matrix_exclude_phage.txt"

        # temp_Only_filtered_positions_for_closely_matrix_exclude_phage = "%s/temp_Only_filtered_positions_for_closely_matrix_exclude_phage.txt" % args.filter2_only_snp_vcf_dir
        temp_Only_filtered_positions_for_closely_matrix_exclude_phage = "%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir
        print temp_Only_filtered_positions_for_closely_matrix_exclude_phage
        # c_reader = csv.reader(open('%s/temp_Only_filtered_positions_for_closely_matrix.txt' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
        c_reader_2 = csv.reader(
            open(temp_Only_filtered_positions_for_closely_matrix_exclude_phage, 'r'), delimiter='\t')
        columns_2 = list(zip(*c_reader_2))
        # print len(columns_2)
        keep_logging('Finished reading columns...', 'Finished reading columns...', logger, 'info')
        counts = 1

        if args.outgroup:
            end = len(vcf_filenames) + 1
            end = end - 1
        else:
            end = len(vcf_filenames) + 1

        f_bar_count = open("%s/bargraph_counts.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_perc = open("%s/bargraph_percentage.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_count.write(
            "Sample\tunmapped_positions\treference_allele\ttrue_variant\tOnly_low_FQ\tOnly_DP\tOnly_low_MQ\tother\n")
        f_bar_perc.write(
            "Sample\tunmapped_positions_perc\ttrue_variant_perc\tOnly_low_FQ_perc\tOnly_DP_perc\tOnly_low_MQ_perc\tother_perc\n")

        for i in xrange(1, end, 1):
            """ Bar Count Statistics: Variant Position Count Statistics """
            true_variant = columns_2[i].count('VARIANT')
            unmapped_positions = columns_2[i].count('reference_unmapped_position')
            # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
            reference_allele = columns_2[i].count('reference_allele') + columns_2[i].count('reference_allele_indel_proximate')
            Only_low_FQ = columns_2[i].count('LowFQ')
            Only_DP = columns_2[i].count('HighFQ_DP')
            Only_low_MQ = columns_2[i].count('HighFQ')
            low_FQ_other_parameters = columns_2[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns_2[i].count(
                'LowFQ_DP_QUAL_proximate_SNP') + columns_2[i].count('LowFQ_QUAL_proximate_SNP') + columns_2[i].count(
                'LowFQ_DP_proximate_SNP') + columns_2[i].count('LowFQ_proximate_SNP') + columns_2[i].count(
                'LowFQ_QUAL_DP') + columns_2[i].count('LowFQ_DP_QUAL') + columns_2[i].count('LowFQ_QUAL') + columns_2[
                                          i].count('LowFQ_DP')
            high_FQ_other_parameters = columns_2[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns_2[i].count(
                'HighFQ_DP_QUAL_proximate_SNP') + columns_2[i].count('HighFQ_QUAL_proximate_SNP') + columns_2[i].count(
                'HighFQ_DP_proximate_SNP') + columns_2[i].count('HighFQ_proximate_SNP') + columns_2[i].count(
                'HighFQ_QUAL_DP') + columns_2[i].count('HighFQ_DP_QUAL') + columns_2[i].count('HighFQ_QUAL')
            other = low_FQ_other_parameters + high_FQ_other_parameters

            total = true_variant + unmapped_positions + reference_allele + Only_low_FQ + Only_DP + low_FQ_other_parameters + high_FQ_other_parameters + Only_low_MQ

            filename_count = i - 1

            if args.outgroup:
                bar_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(
                    vcf_filenames_outgroup[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                                                                   unmapped_positions, reference_allele, true_variant,
                                                                   Only_low_FQ, Only_DP, Only_low_MQ, other)
                f_bar_count.write(bar_string)
            else:
                bar_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(
                    vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                                                                   unmapped_positions, reference_allele, true_variant,
                                                                   Only_low_FQ, Only_DP, Only_low_MQ, other)
            # f_bar_count.write(bar_string)
            """ Bar Count Percentage Statistics: Variant Position Percentage Statistics """
            try:
                true_variant_perc = float((columns_2[i].count('VARIANT') * 100) / total)
            except ZeroDivisionError:
                true_variant_perc = 0
            try:
                unmapped_positions_perc = float((columns_2[i].count('reference_unmapped_position') * 100) / total)
            except ZeroDivisionError:
                unmapped_positions_perc = 0
            try:
                # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
                reference_allele_perc = float(((columns_2[i].count('reference_allele') + columns_2[i].count('reference_allele_indel_proximate')) * 100) / total)
            except ZeroDivisionError:
                reference_allele_perc = 0
            try:
                Only_low_FQ_perc = float((columns_2[i].count('LowFQ') * 100) / total)
            except ZeroDivisionError:
                Only_low_FQ_perc = 0
            try:
                Only_DP_perc = float((columns_2[i].count('HighFQ_DP') * 100) / total)
            except ZeroDivisionError:
                Only_DP_perc = 0
            try:
                Only_low_MQ_perc = float((columns_2[i].count('HighFQ') * 100) / total)
            except ZeroDivisionError:
                Only_low_MQ_perc = 0
            try:
                low_FQ_other_parameters_perc = float(((columns_2[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns_2[
                    i].count('LowFQ_DP_QUAL_proximate_SNP') + columns_2[i].count('LowFQ_QUAL_proximate_SNP') +
                                                       columns_2[i].count('LowFQ_DP_proximate_SNP') + columns_2[
                                                           i].count('LowFQ_proximate_SNP') + columns_2[i].count(
                            'LowFQ_QUAL_DP') + columns_2[i].count('LowFQ_DP_QUAL') + columns_2[i].count('LowFQ_QUAL') +
                                                       columns_2[i].count('LowFQ_DP')) * 100) / total)
            except ZeroDivisionError:
                low_FQ_other_parameters_perc = 0
            try:
                high_FQ_other_parameters_perc = float(((columns_2[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns_2[
                    i].count('HighFQ_DP_QUAL_proximate_SNP') + columns_2[i].count('HighFQ_QUAL_proximate_SNP') +
                                                        columns_2[i].count('HighFQ_DP_proximate_SNP') + columns_2[
                                                            i].count('HighFQ_proximate_SNP') + columns_2[i].count(
                            'HighFQ_QUAL_DP') + columns_2[i].count('HighFQ_DP_QUAL') + columns_2[i].count(
                            'HighFQ_QUAL')) * 100) / total)
            except ZeroDivisionError:
                high_FQ_other_parameters_perc = 0

            other_perc = float(low_FQ_other_parameters_perc + high_FQ_other_parameters_perc)
            if args.outgroup:
                bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(
                    vcf_filenames_outgroup[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                                                                    unmapped_positions_perc, true_variant_perc,
                                                                    Only_low_FQ_perc, Only_DP_perc, Only_low_MQ_perc,
                                                                    other_perc)
            else:
                bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(
                    vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                                                                    unmapped_positions_perc, true_variant_perc,
                                                                    Only_low_FQ_perc, Only_DP_perc, Only_low_MQ_perc,
                                                                    other_perc)
            f_bar_count.write(bar_string)
            f_bar_perc.write(bar_perc_string)
        f_bar_count.close()
        f_bar_perc.close()
        bargraph_R_script = "library(ggplot2)\nlibrary(reshape)\nx1 <- read.table(\"bargraph_percentage.txt\", header=TRUE)\nx1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))\nmdf1=melt(x1,id.vars=\"Sample\")\npdf(\"%s/%s_barplot.pdf\", width = 30, height = 30)\nggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat=\"identity\") + ylab(\"Percentage of Filtered Positions\") + xlab(\"Samples\") + theme(text = element_text(size=9)) + scale_fill_manual(name=\"Reason for filtered out positions\", values=c(\"#08306b\", \"black\", \"orange\", \"darkgrey\", \"#fdd0a2\", \"#7f2704\")) + ggtitle(\"Title Here\") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = \"horizontal\")\ndev.off()" % (
        "%s/matrices/plots" % data_matrix_dir, os.path.basename(os.path.normpath(args.results_dir)))
        barplot_R_file = open("%s/bargraph.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        keep_logging('Run this R script to generate bargraph plot: %s/bargraph.R' % args.filter2_only_snp_vcf_dir,
                     'Run this R script to generate bargraph plot: %s/bargraph.R' % args.filter2_only_snp_vcf_dir,
                     logger, 'info')

    

    # Commented out for debugging
    """ Methods Steps"""
    Only_ref_variant_positions_for_closely = generate_position_label_data_matrix_All_label()
    
    # temp_generate_position_label_data_matrix_All_label()
    # keep_logging('Running: Generating Barplot statistics data matrices...',
    #              'Running: Generating Barplot statistics data matrices...', logger, 'info')
    # barplot_stats()

    method_time_taken = datetime.now() - method_start_time

    # keep_logging('- Time taken to complete the generate_position_label_data_matrix method: {}'.format(method_time_taken),
    #              '- Time taken to complete the generate_position_label_data_matrix method: {}'.format(method_time_taken), logger, 'info')
    return Only_ref_variant_positions_for_closely

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
    method_start_time = datetime.now()

    print_string_header = "\t"
    if args.outgroup:
        for i in vcf_filenames:
            if "%s_filter2_final.vcf_no_proximate_snp.vcf" % outgroup not in i:
                print_string_header = print_string_header + os.path.basename(i) + "\t"
    else:
        for i in vcf_filenames:
            print_string_header = print_string_header + os.path.basename(i) + "\t"

    def generate_indel_position_label_data_matrix_All_label():
        othercodes = ['0', '2', '3', '4', '5', '6', '7']
        ref_var = ['1', '1TRUE']
        f33 = open("%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f33.write('\t' + print_string_header.strip() + '\n')

        if args.outgroup:
            All_indel_label_final_sorted_header_outgroup = pd.read_csv("All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir, sep='\t', header=0)
            Only_ref_indel_positions_for_closely_outgroup = []
            Only_filtered_indel_positions_for_closely_outgroup = []
            for x in All_indel_label_final_sorted_header_outgroup.itertuples():
                if set(ref_var) & set(x[2:]):
                    if set(othercodes) & set(x[2:]):
                        if x[1] not in outgroup_indel_specific_positions:
                            Only_filtered_indel_positions_for_closely_outgroup.append(x[1])
                            for i in row[1:]:
                                print_string = print_string + "\t" + i
                            STRR2 = row[0] + print_string + "\n"
                            f33.write(STRR2)
                    else:
                        if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
                            if x[1] not in functional_filter_pos_array:
                                Only_ref_indel_positions_for_closely_outgroup.append(x[1])
                        else:
                            Only_ref_indel_positions_for_closely_outgroup.append(x[1])
        else:
            All_indel_label_final_sorted_header = pd.read_csv("%s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, sep='\t', header=0)
            Only_ref_indel_positions_for_closely = []
            Only_filtered_indel_positions_for_closely = []
            for x in All_indel_label_final_sorted_header.itertuples():
                if set(ref_var) & set(x[2:]):
                    if set(othercodes) & set(x[2:]):
                        Only_filtered_indel_positions_for_closely.append(x[1])
                        print_string = ""
                        for i in x[2:]:
                            print_string = print_string + "\t" + i
                        STRR2 = str(x[1]) + print_string + "\n"
                        f33.write(STRR2)
                    else:
                        if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
                            if x[1] not in functional_filter_pos_array:
                                Only_ref_indel_positions_for_closely.append(x[1])
                        else:
                            Only_ref_indel_positions_for_closely.append(x[1])
        f33.close()
        
        Only_ref_indel_variant_positions_for_closely_without_functional_filtered_positions = open('%s/Only_ref_indel_variant_positions_for_closely_without_functional_filtered_positions' % args.filter2_only_snp_vcf_dir, 'w')
        for pos in Only_ref_indel_positions_for_closely:
            Only_ref_indel_variant_positions_for_closely_without_functional_filtered_positions.write(str(pos) + '\n')
        Only_ref_indel_variant_positions_for_closely_without_functional_filtered_positions.close()

        return Only_ref_indel_positions_for_closely
    
    def barplot_indel_stats():
        keep_logging(
            'Read each Sample columns and calculate the percentage of each label to generate barplot statistics.',
            'Read each Sample columns and calculate the percentage of each label to generate barplot statistics.',
            logger, 'info')
        """
        Read each Sample columns and calculate the percentage of each label to generate barplot statistics.
        This will give a visual explanation of how many positions in each samples were filtered out because of different reason
        """

        c_reader = csv.reader(
            open('%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt' % args.filter2_only_snp_vcf_dir,
                 'r'), delimiter='\t')
        columns = list(zip(*c_reader))

        keep_logging('Finished reading columns...', 'Finished reading columns...', logger, 'info')
        counts = 1

        if args.outgroup:
            end = len(vcf_filenames) + 1
            end = end - 1
        else:
            end = len(vcf_filenames) + 1
        # print end

        f_bar_count = open("%s/bargraph_indel_counts.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_perc = open("%s/bargraph_indel_percentage.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_count.write(
            "Sample\tunmapped_positions\treference_allele\ttrue_variant\tOnly_low_AF\tOnly_DP\tOnly_low_MQ\tother\n")
        f_bar_perc.write(
            "Sample\tunmapped_positions_perc\ttrue_variant_perc\tOnly_low_AF_perc\tOnly_DP_perc\tOnly_low_MQ_perc\tother_perc\n")
        for i in xrange(1, end, 1):
            """ Bar Count Statistics: Variant Position Count Statistics """
            true_variant = columns[i].count('VARIANT')
            unmapped_positions = columns[i].count('reference_unmapped_position')
            # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
            reference_allele = columns[i].count('reference_allele') + columns[i].count('reference_allele_indel_proximate')
            Only_low_AF = columns[i].count('LowAF')
            Only_DP = columns[i].count('HighAF_DP')
            Only_low_MQ = columns[i].count('HighAF')
            low_AF_other_parameters = columns[i].count('LowAF_QUAL_DP_proximate_SNP') + columns[i].count(
                'LowAF_DP_QUAL_proximate_SNP') + columns[i].count('LowAF_QUAL_proximate_SNP') + columns[i].count(
                'LowAF_DP_proximate_SNP') + columns[i].count('LowAF_proximate_SNP') + columns[i].count(
                'LowAF_QUAL_DP') + columns[i].count('LowAF_DP_QUAL') + columns[i].count('LowAF_QUAL') + columns[
                                          i].count('LowAF_DP')
            high_AF_other_parameters = columns[i].count('HighAF_QUAL_DP_proximate_SNP') + columns[i].count(
                'HighAF_DP_QUAL_proximate_SNP') + columns[i].count('HighAF_QUAL_proximate_SNP') + columns[i].count(
                'HighAF_DP_proximate_SNP') + columns[i].count('HighAF_proximate_SNP') + columns[i].count(
                'HighAF_QUAL_DP') + columns[i].count('HighAF_DP_QUAL') + columns[i].count('HighAF_QUAL')
            other = low_AF_other_parameters + high_AF_other_parameters
            total = true_variant + unmapped_positions + reference_allele + Only_low_AF + Only_DP + low_AF_other_parameters + high_AF_other_parameters + Only_low_MQ
            filename_count = i - 1
            # bar_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(vcf_filenames_outgroup[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), unmapped_positions, reference_allele, true_variant, Only_low_AF, Only_DP, Only_low_MQ, other)
            if args.outgroup:
                ###

                bar_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(
                    vcf_filenames_outgroup[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                                                                   unmapped_positions, reference_allele, true_variant,
                                                                   Only_low_AF, Only_DP, Only_low_MQ, other)
            else:
                bar_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(
                    vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                                                                   unmapped_positions, reference_allele, true_variant,
                                                                   Only_low_AF, Only_DP, Only_low_MQ, other)

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
                # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
                reference_allele_perc = float(((columns[i].count('reference_allele') + columns[i].count('reference_allele_indel_proximate')) * 100) / total)
            except ZeroDivisionError:
                reference_allele_perc = 0
            try:
                Only_low_AF_perc = float((columns[i].count('LowAF') * 100) / total)
            except ZeroDivisionError:
                Only_low_AF_perc = 0
            try:
                Only_DP_perc = float((columns[i].count('HighAF_DP') * 100) / total)
            except ZeroDivisionError:
                Only_DP_perc = 0
            try:
                Only_low_MQ_perc = float((columns[i].count('HighAF') * 100) / total)
            except ZeroDivisionError:
                Only_low_MQ_perc = 0
            try:
                low_AF_other_parameters_perc = float(((columns[i].count('LowAF_QUAL_DP_proximate_SNP') + columns[
                    i].count('LowAF_DP_QUAL_proximate_SNP') + columns[i].count('LowAF_QUAL_proximate_SNP') + columns[
                                                           i].count('LowAF_DP_proximate_SNP') + columns[i].count(
                    'LowAF_proximate_SNP') + columns[i].count('LowAF_QUAL_DP') + columns[i].count('LowAF_DP_QUAL') +
                                                       columns[i].count('LowAF_QUAL') + columns[i].count(
                            'LowAF_DP')) * 100) / total)
            except ZeroDivisionError:
                low_AF_other_parameters_perc = 0
            try:
                high_AF_other_parameters_perc = float(((columns[i].count('HighAF_QUAL_DP_proximate_SNP') + columns[
                    i].count('HighAF_DP_QUAL_proximate_SNP') + columns[i].count('HighAF_QUAL_proximate_SNP') + columns[
                                                            i].count('HighAF_DP_proximate_SNP') + columns[i].count(
                    'HighAF_proximate_SNP') + columns[i].count('HighAF_QUAL_DP') + columns[i].count('HighAF_DP_QUAL') +
                                                        columns[i].count('HighAF_QUAL')) * 100) / total)
            except ZeroDivisionError:
                high_AF_other_parameters_perc = 0

            other_perc = float(low_AF_other_parameters_perc + high_AF_other_parameters_perc)
            if args.outgroup:
                ###
                bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    os.path.basename(
                        vcf_filenames_outgroup[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                    unmapped_positions_perc, true_variant_perc, Only_low_AF_perc, Only_DP_perc, Only_low_MQ_perc,
                    other_perc)
                f_bar_perc.write(bar_perc_string)
            else:
                bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    os.path.basename(
                        vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                    unmapped_positions_perc, true_variant_perc, Only_low_AF_perc, Only_DP_perc, Only_low_MQ_perc,
                    other_perc)
                f_bar_perc.write(bar_perc_string)

        f_bar_count.close()
        f_bar_perc.close()
        bargraph_R_script = "library(ggplot2)\nlibrary(reshape)\nx1 <- read.table(\"bargraph_indel_percentage.txt\", header=TRUE)\nx1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))\nmdf1=melt(x1,id.vars=\"Sample\")\npdf(\"%s/%s_barplot_indel.pdf\", width = 30, height = 30)\nggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat=\"identity\") + ylab(\"Percentage of Filtered Positions\") + xlab(\"Samples\") + theme(text = element_text(size=9)) + scale_fill_manual(name=\"Reason for filtered out positions\", values=c(\"#08306b\", \"black\", \"orange\", \"darkgrey\", \"#fdd0a2\", \"#7f2704\")) + ggtitle(\"Title Here\") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = \"horizontal\")\ndev.off()" % (
        "%s/matrices/plots" % data_matrix_dir, os.path.basename(os.path.normpath(args.results_dir)))
        barplot_R_file = open("%s/bargraph_indel.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        keep_logging('Run this R script to generate bargraph plot: %s/bargraph_indel.R' % args.filter2_only_snp_vcf_dir,
                     'Run this R script to generate bargraph plot: %s/bargraph_indel.R' % args.filter2_only_snp_vcf_dir,
                     logger, 'info')

    """ Methods Steps"""
    keep_logging('- Generating data matrices.', '- Generating data matrices.', logger, 'info')
    Only_ref_indel_positions_for_closely = generate_indel_position_label_data_matrix_All_label()
    
    # Turning off and removing Only AF and Only DP heatmap plot matrices.
    # Not super useful - Run it as an additional script.
    #temp_generate_indel_position_label_data_matrix_All_label()
    # keep_logging('Running: Generating Barplot statistics data matrices...',
    #              'Running: Generating Barplot statistics data matrices...', logger, 'info')
    # barplot_indel_stats()
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('- Time taken to complete the generate_indel_position_label_data_matrix method: {}'.format(method_time_taken),
    #              '- Time taken to complete the generate_indel_position_label_data_matrix method: {}'.format(method_time_taken), logger, 'info')
    return Only_ref_indel_positions_for_closely

def create_job_fasta(jobrun, vcf_filenames, core_vcf_fasta_dir, functional_filter, script_Directive, job_name_flag):
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
        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (
            os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
            core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/temp_jobs" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            call("qsub %s" % i, logger)


    elif jobrun == "parallel-local" or jobrun == "cluster":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_array = []
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
                core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/temp_jobs" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()
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
            num_cores = int(args.numcores)
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

    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_array = []
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
                core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/temp_jobs" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

        # os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
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
        # os.system("bash command_file")
        call("bash %s" % command_file, logger)

def create_job_allele_variant_fasta(jobrun, vcf_filenames, core_vcf_fasta_dir, config_file, script_Directive,
                                    job_name_flag):
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
        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
            os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
            core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

        # os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_fasta.pbs"
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
        command_file = "%s/commands_list_ref_allele_variants_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
            os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
            core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

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
            num_cores = int(args.numcores)
        else:
            # Slurm Changes here.
            if args.scheduler == "SLURM":
                proc = subprocess.Popen(["echo $SLURM_CPUS_PER_TASK"], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                num_cores = int(out.strip())
            elif args.scheduler == "PBS":
                num_cores = multiprocessing.cpu_count()

        
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)

    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_array = []
        command_file = "%s/commands_list_ref_allele_variants_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
                core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()
        # os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
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
        # os.system("bash command_file")
        call("bash %s" % command_file, logger)

def create_job_DP(jobrun, vcf_filenames, script_Directive, job_name_flag):
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
        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
                core_vcf_fasta_dir, config_file)
            command = "python %s/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (
            os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/temp_jobs" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

        # os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
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
        command_file = "%s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')

        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
                core_vcf_fasta_dir, config_file)
            command = "python %s/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (
            os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/temp_jobs" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()

        # os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
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
        # print len(command_array)
        if args.numcores:
            num_cores = int(args.numcores)
        else:
            # Slurm Changes here.
            if args.scheduler == "SLURM":
                proc = subprocess.Popen(["echo $SLURM_CPUS_PER_TASK"], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                num_cores = int(out.strip())
            elif args.scheduler == "PBS":
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
        ### Great Lakes changes
        for i in vcf_filenames:
            command = "python %s/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i, args.reference,
                core_vcf_fasta_dir, config_file)
            command = "python %s/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (
                os.path.dirname(os.path.abspath(__file__)), args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(i))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/temp_jobs" % args.filter2_only_snp_vcf_dir + '\n')
                out.write(command + '\n')
            out.close()
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        # os.system("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir)
        call("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir, logger)

def generate_vcf_files(Only_ref_variant_positions_for_closely):
    method_start_time = datetime.now()

    

    functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % args.filter2_only_snp_vcf_dir
    functional_filter_pos_array = pd.read_csv(functional_class_filter_positions, sep='\n', header=None)
    functional_filter_pos_array = functional_filter_pos_array.squeeze()
    exclude_ref_var_functional = pd.Series(np.intersect1d(Only_ref_variant_positions_for_closely.values,functional_filter_pos_array.values))
    Only_ref_variant_positions_for_closely_without_functional_filtered_positions = Only_ref_variant_positions_for_closely[~Only_ref_variant_positions_for_closely.isin(exclude_ref_var_functional)]
    Only_ref_variant_positions_for_closely_without_functional_filtered_positions.to_csv('Only_ref_variant_positions_for_closely_without_functional_filtered_positions', index=False, sep='\n')

    print "- No. of core SNPs: %s" % len(Only_ref_variant_positions_for_closely_without_functional_filtered_positions)
    filter2_files_array = []
    for i in vcf_filenames:
        filter2_file = i.replace('_no_proximate_snp.vcf', '')
        filter2_files_array.append(filter2_file)

    filtered_out_vcf_files = []
    for i in filter2_files_array:
        print_array = []
        with open(i) as file_open:
            for line in file_open:
                line = line.strip()
                if line.startswith("#"):
                    print_array.append(line)
                else:
                    split_array = re.split(r'\t+', line)
                    if int(split_array[1]) in Only_ref_variant_positions_for_closely_without_functional_filtered_positions and 'INDEL' not in split_array[7]:
                        print_array.append(line)
        file_open.close()
        file_name = i + "_core.vcf"
        filtered_out_vcf_files.append(file_name)
        f1 = open(file_name, 'w+')
        for ios in print_array:
            print_string = str(ios) + "\n"
            f1.write(print_string)
        f1.close()

    # Turning off generating core fasta alignemnets. No longer used in pipeline
    filename = "%s/consensus.sh" % args.filter2_only_snp_vcf_dir
    # keep_logging('- Extracting Core Genome.', '- Extracting Core Genome.', logger, 'info')
    for file in filtered_out_vcf_files:
        f1 = open(filename, 'a+')
        bgzip_cmd = "bgzip -f %s\n" % (file)
        f1.write(bgzip_cmd)
        subprocess.call([bgzip_cmd], shell=True)
        tabix_cmd = "tabix -f -p vcf %s.gz\n" % (file)
        f1.write(tabix_cmd)
        subprocess.call([tabix_cmd], shell=True)
        fasta_cmd = "cat %s | vcf-consensus %s.gz > %s.fa\n" % (
        args.reference, file, file.replace('_filter2_final.vcf_core.vcf', ''))
        f1.write(fasta_cmd)
        subprocess.call([fasta_cmd], shell=True)
        base = os.path.basename(file)
        header = base.replace('_filter2_final.vcf_core.vcf', '')
        sed_command = "sed -i 's/>.*/>%s/g' %s.fa\n" % (header, file.replace('_filter2_final.vcf_core.vcf', ''))
        subprocess.call([sed_command], shell=True)
        f1.write(sed_command)
    # sequence_lgth_cmd = "for i in %s/*.fa; do %s/%s/bioawk -c fastx \'{ print $name, length($seq) }\' < $i; done" % (args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bioawk", Config)['bioawk_bin'])
    # #os.system(sequence_lgth_cmd)
    # call("%s" % sequence_lgth_cmd, logger)
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('Time taken to complete the core genome method: {}'.format(method_time_taken),
    #              'Time taken to complete the core genome method: {}'.format(method_time_taken), logger, 'info')

def gatk_filter2(final_raw_vcf, out_path, analysis, reference):
    gatk_filter2_parameter_expression = "MQ > 50 && QUAL > 100 && DP > 9"
    gatk_filter2_command = "java -jar %s/%s/GenomeAnalysisTK.jar -T VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s --filterExpression \"%s\" --filterName PASS_filter2" % (
    ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("gatk", Config)['gatk_bin'], reference, out_path,
    analysis, final_raw_vcf, gatk_filter2_parameter_expression)
    keep_logging('Running Command: [%s]' % gatk_filter2_command, 'Running Command: [%s]' % gatk_filter2_command, logger,
                 'info')
    # os.system(gatk_filter2_command)
    call("%s" % gatk_filter2_command, logger)
    filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_gatk.vcf > %s/%s_filter2_final.vcf" % (
    out_path, analysis, out_path, analysis)
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
                # print position + "  " + all_position[next_position_index]
                if position not in remove_proximate_position_array and all_position[
                    next_position_index] not in remove_proximate_position_array:
                    remove_proximate_position_array.append(int(position))
                    remove_proximate_position_array.append(int(all_position[next_position_index]))
    f1 = open(gatk_filter2_final_vcf_file_no_proximate_snp, 'w+')
    with open(gatk_filter2_final_vcf_file, 'rU') as csv_file2:
        for line in csv_file2:
            if line.startswith('gi') or line.startswith('MRSA_8058'):  ##change this!
                line_array = line.split('\t')
                if int(line_array[1]) not in remove_proximate_position_array:
                    print_string = line
                    f1.write(print_string)
            else:
                print_string = line
                f1.write(print_string)
    gatk_filter2_final_vcf_file_no_proximate_snp_positions = gatk_filter2_final_vcf_file + "_no_proximate_snp.vcf_positions_array"
    f2 = open(gatk_filter2_final_vcf_file_no_proximate_snp_positions, 'w+')
    for i in remove_proximate_position_array:
        position_print_string = str(i) + "\n"
        f2.write(position_print_string)
    return gatk_filter2_final_vcf_file_no_proximate_snp

def FQ_analysis():
    for i in vcf_filenames:
        filename_base = os.path.basename(i)
        aln_mpileup_vcf_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                         '_aln_mpileup_raw.vcf_5bp_indel_removed.vcf')
        analysis = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')
        # print aln_mpileup_vcf_file
        grep_reference_file = "grep \'^##reference\' %s" % aln_mpileup_vcf_file
        proc = subprocess.Popen([grep_reference_file], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        reference_file = out.split(':')
        # Change it to multiprocessing
        gatk_filter2_final_vcf_file = gatk_filter2(aln_mpileup_vcf_file, temp_dir, analysis, reference_file[1])
        # print gatk_filter2_final_vcf_file
        gatk_filter2_final_vcf_file_no_proximate_snp = remove_proximate_snps(gatk_filter2_final_vcf_file, temp_dir,
                                                                             analysis, reference_file[1])
        grep_fq_field = "awk -F\'\\t\' \'{print $8}\' %s | grep -o \'FQ=.*\' | sed \'s/FQ=//g\' | awk -F\';\' \'{print $1}\' > %s/%s_FQ_values" % (
        gatk_filter2_final_vcf_file_no_proximate_snp, os.path.dirname(i), analysis)
        # os.system(grep_fq_field)
        call("%s" % grep_fq_field, logger)
        # print grep_fq_field

def DP_analysis():
    create_job_DP(args.jobrun, vcf_filenames, script_Directive, job_name_flag)
    paste_command = "paste %s/extract_DP_positions.txt" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_DP_values')
        paste_command = paste_command + " " + label_file

    paste_file = args.filter2_only_snp_vcf_dir + "/paste_DP_files.sh"
    f2 = open(paste_file, 'w+')
    paste_command = paste_command + " > %s/filtered_DP_values_temp.txt" % args.filter2_only_snp_vcf_dir
    # os.system(paste_command)
    f2.write(paste_command + '\n')
    cat_header = "cat %s/header.txt %s/filtered_DP_values_temp.txt > %s/filtered_DP_values.txt" % (
    args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    # os.system(cat_header)
    f2.write(cat_header + '\n')
    sed_command = "sed -i \'s/_filter2_final.vcf_no_proximate_snp.vcf//g\' %s/filtered_DP_values.txt" % (
        args.filter2_only_snp_vcf_dir)
    # os.system(sed_command)
    f2.write(sed_command + '\n')
    cmd = "bash %s" % paste_file
    # os.system("bash %s/paste_DP_files.sh" % args.filter2_only_snp_vcf_dir)

def DP_analysis_barplot():
    # os.system("bash %s/paste_DP_files.sh" % args.filter2_only_snp_vcf_dir)
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
        bar_string = "%s\t%s\t%s\t%s\t%s\t%s\n" % (
        os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
        reference_position, oneto5, sixto10, elevento14, fifteenorabove)
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
        bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\n" % (
        os.path.basename(vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
        reference_position_perc, oneto5_perc, sixto10_perc, elevento14_perc, fifteenorabove_perc)
        f_bar_perc.write(bar_perc_string)


""" Deprecated Method """


def extract_only_ref_variant_fasta(core_vcf_fasta_dir):
    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes" and \
            ConfigSectionMap("functional_filters", Config)['apply_to_calls'] == "yes":
        functional_filter = "yes"
    create_job_fasta(args.jobrun, vcf_filenames, core_vcf_fasta_dir, functional_filter, script_Directive, job_name_flag)


""" Deprecated Method """


def extract_only_ref_variant_fasta_from_reference():
    method_start_time = datetime.now()
    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes" and \
            ConfigSectionMap("functional_filters", Config)['apply_to_calls'] == "yes":
        ffp = open(
            "%s/Only_ref_variant_positions_for_closely_without_functional_filtered_positions" % args.filter2_only_snp_vcf_dir).readlines()
    else:
        ffp = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir).readlines()
    fasta_string = ""
    # firstLine = ffp.pop(0)
    for lines in ffp:
        lines = lines.strip()
        extract_base = "grep -v \'>\' %s | tr -d \'\\n\'| cut -b%s" % (args.reference, lines)
        proc = subprocess.Popen([extract_base], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        fasta_string = fasta_string + out
        if not out:
            # print lines
            keep_logging('Error extracting reference allele', 'Error extracting reference allele', logger, 'info')
            exit()

    pattern = re.compile(r'\s+')
    fasta_string = re.sub(pattern, '', fasta_string)
    final_fasta_string = ">%s\n" % os.path.basename(
        args.reference.replace('.fasta', '').replace('.fa', '')) + fasta_string + "\n"
    fp = open("%s/%s_variants.fa" % (
    args.filter2_only_snp_vcf_dir, os.path.basename(args.reference.replace('.fasta', '').replace('.fa', ''))), 'w+')
    fp.write(final_fasta_string)
    fp.close()
    method_time_taken = datetime.now() - method_start_time

    keep_logging('Time taken to complete the extract_only_ref_variant_fasta_from_reference method: {}'.format(method_time_taken),
                 'Time taken to complete the extract_only_ref_variant_fasta_from_reference method: {}'.format(method_time_taken), logger, 'info')

def extract_only_ref_variant_fasta_from_reference_allele_variant():
    method_start_time = datetime.now()
    ffp = open("%s/unique_positions_file" % args.filter2_only_snp_vcf_dir).readlines()
    # unique_positions_array = []

    fasta_string = ""
    # firstLine = ffp.pop(0)
    for lines in ffp:
        lines = lines.strip()
        # unique_positions_array.append(lines)
        extract_base = "grep -v \'>\' %s | tr -d \'\\n\'| cut -b%s" % (args.reference, lines)
        proc = subprocess.Popen([extract_base], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.strip()
        fasta_string = fasta_string + out
        if not out:
            # print lines
            keep_logging('Error extracting reference allele', 'Error extracting reference allele', logger, 'info')
            exit()

    pattern = re.compile(r'\s+')
    fasta_string = re.sub(pattern, '', fasta_string)
    final_fasta_string = ">%s\n" % os.path.basename(
        args.reference.replace('.fasta', '').replace('.fa', '')) + fasta_string + "\n"
    fp = open("%s/%s_allele_variants.fa" % (
    args.filter2_only_snp_vcf_dir, os.path.basename(args.reference.replace('.fasta', '').replace('.fa', ''))), 'w+')
    fp.write(final_fasta_string)
    fp.close()
    method_time_taken = datetime.now() - method_start_time

    keep_logging('Time taken to complete the extract_only_ref_variant_fasta_from_reference_allele_variant method: {}'.format(method_time_taken),
                 'Time taken to complete the extract_only_ref_variant_fasta_from_reference_allele_variant method: {}'.format(method_time_taken), logger, 'info')

def core_prep_snp():
    method_start_time = datetime.now()
    """ Generate SNP Filter Label Matrix """
    generate_paste_command()

    generate_paste_command_outgroup()

    """ Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix. """
    Only_ref_variant_positions_for_closely = generate_position_label_data_matrix()

    """ Generate VCF files from final list of variants in Only_ref_variant_positions_for_closely; generate commands for consensus generation """
    generate_vcf_files(Only_ref_variant_positions_for_closely)

    """ Analyze the positions that were filtered out only due to insufficient depth"""
    # DP_analysis()
    method_time_taken = datetime.now() - method_start_time

    keep_logging('- Time taken to parse Single Variant VCFs: {}'.format(method_time_taken),
                 '- Time taken to parse Single Variant VCFs: {}'.format(method_time_taken), logger, 'info')
    return Only_ref_variant_positions_for_closely

def core_prep_indel():
    method_start_time = datetime.now()
    """ Generate SNP Filter Label Matrix """
    generate_indel_paste_command()

    generate_indel_paste_command_outgroup()

    """ Generate different list of Positions from the **All_label_final_sorted_header.txt** SNP position label data matrix. """
    Only_ref_indel_positions_for_closely = generate_indel_position_label_data_matrix()
    method_time_taken = datetime.now() - method_start_time

    keep_logging('- Time taken to parse Indel VCFs: {}'.format(method_time_taken),
                 '- Time taken to parse Indel VCFs: {}'.format(method_time_taken), logger, 'info')
    return Only_ref_indel_positions_for_closely
""" Annotation methods"""

def prepare_snpEff_db(reference_basename):
    method_start_time = datetime.now()
    keep_logging('Preparing snpEff database requirements.', 'Preparing snpEff database requirements.', logger, 'info')
    reference_basename = (os.path.basename(args.reference)).split(".")

    ## Great Lakes Changes
    proc = subprocess.Popen(["find $CONDA_PREFIX/share/ -name snpEff.config"], stdout=subprocess.PIPE, shell=True)
    (out2, err2) = proc.communicate()
    if out2:
        snpeff_config = (str(out2)).strip()
    else:
        print "Unable to find snpEff config file in conda Environment share directory"
        exit()

    # os.system("cp %s $CONDA_PREFIX/bin/" % snpeff_config)
    os.system("cp %s %s" % (snpeff_config, bin_dir))

    if os.path.isfile("%s/snpEff.config" % bin_dir):
        # os.system("cp %s/%s/snpEff.config %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], args.filter2_only_snp_vcf_dir))
        keep_logging("cp %s/snpEff.config %s" % (bin_dir, args.filter2_only_snp_vcf_dir),
                     "cp %s/snpEff.config %s" % (bin_dir, args.filter2_only_snp_vcf_dir), logger, 'debug')
        call("cp %s/snpEff.config %s" % (bin_dir, args.filter2_only_snp_vcf_dir), logger)
    else:
        keep_logging("Error: %s/snpEff.config doesn't exists.\nExiting..." % bin_dir,
                     "Error: %s/snpEff.config doesn't exists.\nExiting..." % bin_dir, logger, 'exception')
        exit()
    make_sure_path_exists("%s/data/%s" % (bin_dir, reference_basename[0]))
    make_sure_path_exists("%s/data/genomes/" % bin_dir)
    # os.system("cp %s %s/%s/data/genomes/" % (args.reference, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']))
    keep_logging("cp %s %s/data/genomes/%s.fa" % (args.reference, bin_dir, reference_basename[0]),
                 "cp %s %s/data/genomes/" % (args.reference, bin_dir), logger, 'debug')
    call("cp %s %s/data/genomes/%s.fa" % (args.reference, bin_dir, reference_basename[0]), logger)
    with open("%s/snpEff.config" % args.filter2_only_snp_vcf_dir, "a") as conf_file:
        conf_file.write(
            "\n\n##Building Custom Database###\n%s.genome\t: %s\n\n" % (reference_basename[0], reference_basename[0]))
    conf_file.close()
    # get the gff name from config file
    if os.path.isfile("%s/%s.gff" % (os.path.dirname(args.reference), reference_basename[0])):
        keep_logging("cp %s/%s.gff %s/data/%s/genes.gff" % (
            os.path.dirname(args.reference), reference_basename[0], bin_dir, reference_basename[0]),
                     "cp %s/%s.gff %s/data/%s/genes.gff" % (os.path.dirname(args.reference), reference_basename[0],
                                                            bin_dir,
                                                            reference_basename[0]), logger, 'debug')
        keep_logging("cp %s/%s.gb* %s/data/%s/genes.gbk" % (
            os.path.dirname(args.reference), reference_basename[0], bin_dir, reference_basename[0]),
                     "cp %s/%s.gff %s/data/%s/genes.gff" % (os.path.dirname(args.reference), reference_basename[0],
                                                            bin_dir,
                                                            reference_basename[0]), logger, 'debug')
        call("cp %s/%s.gff %s/data/%s/genes.gff" % (
            os.path.dirname(args.reference), reference_basename[0], bin_dir, reference_basename[0]), logger)
        call("cp %s/%s.gb* %s/data/%s/genes.gbk" % (
            os.path.dirname(args.reference), reference_basename[0], bin_dir, reference_basename[0]), logger)
    else:
        keep_logging(
            "Error: %s/%s.gff file doesn't exists. Make sure the GFF file has the same prefix as reference fasta file\nExiting..." % (
            os.path.dirname(args.reference), reference_basename[0]),
            "Error: %s/%s.gff file doesn't exists. Make sure the GFF file has the same prefix as reference fasta file\nExiting..." % (
            os.path.dirname(args.reference), reference_basename[0]), logger, 'exception')
        exit()
    # keep_logging("java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), "java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger, 'debug')
    # keep_logging("java -jar %s/%s/%s build -genbank -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), "java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger, 'debug')
    ## Great Lakes Changes
    keep_logging("%s build -genbank -v %s -c %s/snpEff.config -dataDir %s/data" % (
    ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, bin_dir),
                 "%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/data" % (
                 ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir,
                 bin_dir), logger, 'debug')

    # call("java -jar %s/%s/%s build -gff3 -v %s -c %s/snpEff.config -dataDir %s/%s/data" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin']), logger)
    ## Great Lakes Changes
    call("%s build -genbank -v %s -c %s/snpEff.config -dataDir %s/data" % (
    ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], args.filter2_only_snp_vcf_dir, bin_dir),
         logger)
    keep_logging('Finished Preparing snpEff database requirements.', 'Finished Preparing snpEff database requirements.',
                 logger, 'info')
    method_time_taken = datetime.now() - method_start_time

    keep_logging('Time taken to complete the prepare_snpEff_db method: {}'.format(method_time_taken),
                 'Time taken to complete the prepare_snpEff_db method: {}'.format(method_time_taken), logger, 'info')

def variant_annotation():
    method_start_time = datetime.now()
    keep_logging('Annotating Variants using snpEff.', 'Annotating Variants using snpEff.', logger, 'info')

    if ConfigSectionMap("snpeff", Config)['prebuild'] == "yes":
        if ConfigSectionMap("snpeff", Config)['db']:
            print "Using pre-built snpEff database: %s" % ConfigSectionMap("snpeff", Config)['db']
            ## Great Lakes Changes
            proc = subprocess.Popen(["%s databases | grep %s" % (
            ConfigSectionMap("snpeff", Config)['base_cmd'], ConfigSectionMap("snpeff", Config)['db'])],
                                    stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            if out2:
                snpeffdb = ConfigSectionMap("snpeff", Config)['db']
            else:
                print "The database name %s provided was not found. Check the name and try again" % \
                      ConfigSectionMap("snpeff", Config)['db']
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
        annotate_vcf_cmd = "%s -csvStats %s_ANN.csv -dataDir %s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("snpeff", Config)['base_cmd'], raw_vcf, bin_dir,
                            ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir,
                            snpeffdb, raw_vcf, raw_vcf)
        # print annotate_vcf_cmd
        annotate_vcf_cmd_array.append(annotate_vcf_cmd)
        final_vcf = i
        annotate_final_vcf_cmd = "%s -csvStats %s_ANN.csv -dataDir %s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                                 (ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, bin_dir,
                                  ConfigSectionMap("snpeff", Config)['snpeff_parameters'],
                                  args.filter2_only_snp_vcf_dir, snpeffdb, final_vcf, final_vcf)
        annotate_final_vcf_cmd_array.append(annotate_final_vcf_cmd)

    # print annotate_vcf_cmd_array
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_vcf_cmd_array)
    results_2 = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_final_vcf_cmd_array)
    method_time_taken = datetime.now() - method_start_time
    keep_logging('Time taken to complete the variant_annotation method: {}'.format(method_time_taken),
                 'Time taken to complete the variant_annotation method: {}'.format(method_time_taken), logger, 'info')

def indel_annotation():
    method_start_time = datetime.now()
    keep_logging('Annotating indels using snpEff.', 'Annotating indels using snpEff.', logger, 'info')

    if ConfigSectionMap("snpeff", Config)['prebuild'] == "yes":
        if ConfigSectionMap("snpeff", Config)['db']:
            print "Using pre-built snpEff database: %s" % ConfigSectionMap("snpeff", Config)['db']
            proc = subprocess.Popen(["%s databases | grep %s" % (
            ConfigSectionMap("snpeff", Config)['base_cmd'], ConfigSectionMap("snpeff", Config)['db'])],
                                    stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            if out2:
                snpeffdb = ConfigSectionMap("snpeff", Config)['db']
            else:
                print "The database name %s provided was not found. Check the name and try again" % \
                      ConfigSectionMap("snpeff", Config)['db']
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
        annotate_vcf_cmd = "%s -csvStats %s_ANN.csv -dataDir %s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                           (ConfigSectionMap("snpeff", Config)['base_cmd'], raw_vcf, bin_dir,
                            ConfigSectionMap("snpeff", Config)['snpeff_parameters'], args.filter2_only_snp_vcf_dir,
                            snpeffdb, raw_vcf, raw_vcf)
        annotate_vcf_cmd_array.append(annotate_vcf_cmd)
        final_vcf = i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf')
        annotate_final_vcf_cmd = "%s -csvStats %s_ANN.csv -dataDir %s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                                 (ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, bin_dir,
                                  ConfigSectionMap("snpeff", Config)['snpeff_parameters'],
                                  args.filter2_only_snp_vcf_dir, snpeffdb, final_vcf, final_vcf)
        annotate_final_vcf_cmd_array.append(annotate_final_vcf_cmd)
    if args.numcores:
        num_cores = int(args.numcores)
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

    
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_vcf_cmd_array)
    results_2 = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in annotate_final_vcf_cmd_array)
    method_time_taken = datetime.now() - method_start_time

    keep_logging('Time taken to complete the indel_annotation method: {}'.format(method_time_taken),
                 'Time taken to complete the indel_annotation method: {}'.format(method_time_taken), logger, 'info')

def gatk_combine_variants(files_gatk, reference, out_path, merged_file_suffix, logger, Config):
    method_start_time = datetime.now()
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)[
        'gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    # files_gatk = "--variant " + ' --variant '.join(vcf_files_array)
    # keep_logging("java -jar %s/GenomeAnalysisTK.jar -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s" % (
    # os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), reference, files_gatk, out_path,
    # merged_file_suffix), "java -jar %s -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s" % (
    #              ConfigSectionMap("gatk", Config)['base_cmd'], reference, files_gatk, out_path, merged_file_suffix),
    #              logger, 'debug')
    merge_gatk_commands_file = "%s/gatk_merge.sh" % args.filter2_only_snp_vcf_dir
    with open(merge_gatk_commands_file, 'w+') as fopen:
        fopen.write("java -jar %s/GenomeAnalysisTK.jar -T CombineVariants -R %s %s -o %s/Final_vcf_gatk%s --log_to_file %s/gatk_combinevariants.log 2> /dev/null 1> %s/gatkoutput.txt" % (
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), reference, files_gatk, out_path,
        merged_file_suffix, out_path, out_path) + '\n')
    fopen.close()
    # Commenting out calling gatk combine variants with a custom logging call method, problem with python subprocess, OSError: [Errno 7] Argument list too long
    os.system("bash %s" % merge_gatk_commands_file)
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('- Time taken to complete GATK combine variants method: {}'.format(method_time_taken),
    #              '- Time taken to complete GATK combine variants method: {}'.format(method_time_taken), logger, 'info')
    return "%s/Final_vcf_gatk%s" % (out_path, merged_file_suffix)

def extract_locus_tag_from_genbank():
    """ Start of Extract Annotation information from Genbank file

        Extract Annotation information from Genbank file

        - Check if Reference genome Genbank file exists.
        - Initiate dictionaries that maps locus tag to gene name and product. This information will be used for annotating SNP/Indel Matrix
        - Read the locus tag and gene annotations into a dictionary that maps locus tags to gene name/product name

        """
    keep_logging('- Extracting Locus Tag information from Genbank.', '- Extracting Locus Tag information from Genbank.',
                 logger, 'info')
    reference_basename = (os.path.basename(args.reference)).split(".")
    if os.path.isfile("%s/%s.gbf" % (os.path.dirname(args.reference), reference_basename[0])):
        handle = open("%s/%s.gbf" % (os.path.dirname(args.reference), reference_basename[0]), 'rU')
    else:
        raise IOError('%s/%s.gbf does not exist.' % (os.path.dirname(args.reference), reference_basename[0]))
        exit()

    global locus_tag_to_gene_name, locus_tag_to_product, locus_tag_to_strand

    locus_tag_to_gene_name = {}
    locus_tag_to_product = {}
    locus_tag_to_strand = {}
    # locus_tag_to_uniprot = {}
    # locus_tag_to_ec_number = {}

    keep_logging(
        '- Reading annotations from Reference genome genbank file: %s/%s.gbf' % (
            os.path.dirname(args.reference), reference_basename[0]),
        '- Reading annotations from Reference genome genbank file: %s/%s.gbf' % (
            os.path.dirname(args.reference), reference_basename[0]),
        logger, 'info')
    for record in SeqIO.parse(handle, 'genbank'):
        for feature in record.features:
            location = str(feature.location)
            strand = location.split('(')[1].replace(')', '')
            if 'locus_tag' in feature.qualifiers:
                locus_tag_to_strand[str(feature.qualifiers['locus_tag'][0])] = strand
                if 'gene' in feature.qualifiers:
                    locus_tag_to_gene_name[str(feature.qualifiers['locus_tag'][0])] = str(feature.qualifiers['gene'][0])
                else:
                    locus_tag_to_gene_name[str(feature.qualifiers['locus_tag'][0])] = "null or hypothetical protein"
                if 'product' in feature.qualifiers:
                    locus_tag_to_product[str(feature.qualifiers['locus_tag'][0])] = str(
                        feature.qualifiers['product'][0])
                else:
                    locus_tag_to_product[str(feature.qualifiers['locus_tag'][0])] = "null or hypothetical protein"
            else:
                continue
                # keep_logging(
                #     '',
                #     'Warning: locus_tag specifications for the below feature doesnt exists. Please check the format of genbank file\n%s' % str(
                #         feature),
                #     logger, 'warning')

    global first_locus_tag, last_element, last_locus_tag

    first_locus_tag = record.features[1].qualifiers['locus_tag'][0]
    last_element = len(record.features) - 1
    last_locus_tag = record.features[last_element].qualifiers['locus_tag'][0]

    # #Debugging prints
    # print first_locus_tag
    # print locus_tag_to_gene_name[first_locus_tag]
    # print last_locus_tag
    # print locus_tag_to_gene_name[last_locus_tag]

    """ End of Extract Annotation information from Genbank file 

        Extract Annotation information from Genbank file 

        - Check if Reference genome Genbank file exists.
        - Initiate dictionaries that maps locus tag to gene name and product. This information will be used for annotating SNP/Indel Matrix
        - Read the locus tag and gene annotations into a dictionary that maps locus tags to gene name/product name

    """

    return locus_tag_to_gene_name, locus_tag_to_product, locus_tag_to_strand, first_locus_tag, last_element, last_locus_tag

def merge_vcf():
    """ Start of Merging Step:

        - Merge Individual Annotated raw and filtered vcf files to generate a Final merged vcf file using Gatk combine variants method.
        - Parse this merged Final_vcf* file and generate a SNP/Indel matrix

        """
    method_start_time = datetime.now()
    keep_logging(
        '- Merging Final Annotated VCF files into Final_vcf_no_proximate_snp.vcf using bcftools',
        '- Merging Final Annotated VCF files into Final_vcf_no_proximate_snp.vcf using bcftools',
        logger, 'info')

    # Commented for SNP Matrix debugging
    # files_for_tabix = glob.glob("%s/*.vcf_no_proximate_snp.vcf_ANN.vcf" % args.filter2_only_snp_vcf_dir)
    # tabix(files_for_tabix, "vcf", logger, Config)
    # files_for_tabix = glob.glob("%s/*_filter2_indel_final.vcf_ANN.vcf" % args.filter2_only_snp_vcf_dir)
    # tabix(files_for_tabix, "vcf", logger, Config)

    files = ' '.join(vcf_filenames)

    """ bcftools merging is deprecated. Replaced with GATK combinevariants """
    merge_commands_file = "%s/bcftools_merge.sh" % args.filter2_only_snp_vcf_dir

    with open(merge_commands_file, 'w+') as fopen:
        fopen.write("bcftools merge -i ANN:join -m both -o %s/Final_vcf_no_proximate_snp.vcf -O v %s" % (
        args.filter2_only_snp_vcf_dir, files.replace("_filter2_final.vcf_no_proximate_snp.vcf",
                                                     "_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz")) + '\n')
        fopen.write(
            "bcftools merge -i ANN:join -m both -o %s/Final_vcf_indel.vcf -O v %s" % (args.filter2_only_snp_vcf_dir,
                                                                                      files.replace(
                                                                                          "_filter2_final.vcf_no_proximate_snp.vcf",
                                                                                          "_filter2_indel_final.vcf_ANN.vcf.gz")) + '\n')

    fopen.close()

    os.system("bash %s" % merge_commands_file)

    """ Merge with Gatk combine variants method """
    # Commented for SNP Matrix debugging
    merged_file_suffix = "_no_proximate_snp.vcf"

    annotated_no_proximate_snp_file = "%s/annotated_no_proximate_snp_list.txt" % args.filter2_only_snp_vcf_dir
    annotated_no_proximate_snp_indel_file = "%s/annotated_no_proximate_snp_indel_list.txt" % args.filter2_only_snp_vcf_dir

    with open(annotated_no_proximate_snp_file, 'w+') as fopen:
        for i in vcf_filenames:
            fopen.write(i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                  '_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz') + '\n')
    fopen.close()

    with open(annotated_no_proximate_snp_indel_file, 'w+') as fopen:
        for i in vcf_filenames:
            fopen.write(
                i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf_ANN.vcf.gz') + '\n')
    fopen.close()

    # files_gatk = "--variant " + ' --variant '.join(vcf_filenames)
    files_gatk = ""
    for i in vcf_filenames:
        files_gatk = files_gatk + " --variant " + i
    final_gatk_snp_merged_vcf = gatk_combine_variants(files_gatk.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                                                         '_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz'),
                                                      args.reference, args.filter2_only_snp_vcf_dir, merged_file_suffix,
                                                      logger, Config)

    # Test this merge and annotate this merged file - Testing Mode Right now.
    # merged_file_suffix = "_no_proximate_snp_1.vcf"
    # final_gatk_snp_merged_vcf_1 = gatk_combine_variants(files_gatk,args.reference, args.filter2_only_snp_vcf_dir, merged_file_suffix, logger, Config)
    merged_file_suffix = "_indel.vcf"
    final_gatk_indel_merged_vcf = gatk_combine_variants(files_gatk.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                                                           '_filter2_indel_final.vcf_ANN.vcf.gz'),
                                                        args.reference, args.filter2_only_snp_vcf_dir,
                                                        merged_file_suffix,
                                                        logger, Config)

    """ Tabix index the combined GATK Final vcf file """
    files_for_tabix = glob.glob("%s/Final_vcf_*.vcf" % args.filter2_only_snp_vcf_dir)
    tabix(files_for_tabix, "vcf", logger, Config)

    """ End of Merging Step. """

    annotated_no_proximate_snp_file = "%s/annotated_no_proximate_snp_list.txt" % args.filter2_only_snp_vcf_dir
    annotated_no_proximate_snp_indel_file = "%s/annotated_no_proximate_snp_indel_list.txt" % args.filter2_only_snp_vcf_dir
    final_gatk_snp_merged_vcf = "Final_vcf_gatk_no_proximate_snp.vcf"
    final_gatk_indel_merged_vcf = "Final_vcf_gatk_indel.vcf"
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('- Time taken to complete the merge_vcf method: {}'.format(method_time_taken),
    #              '- Time taken to complete the merge_vcf method: {}'.format(method_time_taken), logger, 'info')
    return annotated_no_proximate_snp_file, annotated_no_proximate_snp_indel_file, final_gatk_snp_merged_vcf, final_gatk_indel_merged_vcf


def extract_annotations_from_multivcf():
    method_start_time = datetime.now()
    """ Extract ANN information from bcftools Final vcf file. (There is a reason why i am using bcftools merged file to extract ANN information) """
    snp_var_ann_dict = {}
    indel_var_ann_dict = {}

    for variants in VCF("%s/Final_vcf_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir):
        snp_var_ann_dict[variants.POS] = variants.INFO.get('ANN')

    # Commented out for legionella bug
    for variants in VCF("%s/Final_vcf_indel.vcf.gz" % args.filter2_only_snp_vcf_dir):
        indel_var_ann_dict[variants.POS] = variants.INFO.get('ANN')
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('- Time taken to complete the extract_annotations_from_multivcf method: {}'.format(method_time_taken),
    #              '- Time taken to complete the extract_annotations_from_multivcf method: {}'.format(method_time_taken), logger, 'info')
    """ End of Extract ANN information from bcftools Final vcf file"""
    return snp_var_ann_dict, indel_var_ann_dict


def extract_core_positions():
    method_start_time = datetime.now()
    """ Generate an array of core positions. Read Only_ref_variant_positions_for_closely* to get final core variant positions into core_positions array"""
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
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('Time taken to complete the extract_core_positions method: {}'.format(method_time_taken),
    #              'Time taken to complete the extract_core_positions method: {}'.format(method_time_taken), logger, 'info')
    """ End: Generate an array of core positions. """
    return core_positions, indel_core_positions


def extract_functional_class_positions():
    """ Generate a list of functional class positions from Phaster, Mummer and Custom Masking results/files"""
    method_start_time = datetime.now()
    """ Read in functional class filter positions. """
    functional_filter_pos_array = []
    with open(functional_class_filter_positions, 'rU') as f_functional:
        for line_func in f_functional:
            functional_filter_pos_array.append(line_func.strip())

    """ GET individual PHAGE/Repetitive/masked region positions to assign functional class group string """
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
            else:
                raise IOError('%s/phage_region_positions.txt does not exist.' % args.filter2_only_snp_vcf_dir)
                exit()
        # GET REPETITIVE REGIONS
        if ConfigSectionMap("functional_filters", Config)['find_repetitive_region'] == "yes":
            repetitive_positions_file = "%s/repeat_region_positions.txt" % args.filter2_only_snp_vcf_dir
            if os.path.isfile(repetitive_positions_file):
                with open(repetitive_positions_file, 'rU') as frep:
                    for line in frep:
                        repetitive_positions.append(line.strip())
                frep.close()
            else:
                raise IOError('%s/repeat_region_positions.txt does not exist.' % args.filter2_only_snp_vcf_dir)
                exit()
        # GET MASK REGIONS
        if ConfigSectionMap("functional_filters", Config)['mask_region'] == "yes":
            mask_positions_file = "%s/mask_positions.txt" % args.filter2_only_snp_vcf_dir
            if os.path.isfile(mask_positions_file):
                with open(mask_positions_file, 'rU') as fmask:
                    for line in fmask:
                        mask_positions.append(line.strip())
                fmask.close()
            else:
                raise IOError('%s/mask_positions.txt does not exist.' % args.filter2_only_snp_vcf_dir)
                exit()

    """ End: Generate a list of functional class positions from Phaster, Mummer and Custom Masking results/files"""
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('Time taken to complete the extract_functional_class_positions method: {}'.format(method_time_taken),
    #              'Time taken to complete the extract_functional_class_positions method: {}'.format(method_time_taken), logger, 'info')
    return functional_filter_pos_array, phage_positions, repetitive_positions, mask_positions

def generate_position_label_dict(final_merge_anno_file):
    method_start_time = datetime.now()
    global position_label, position_indel_label

    """ Prepare a *final_ordered_sorted.txt file with sorted unique variant positions. """
    paste_label_command = "paste %s/unique_positions_file " % args.filter2_only_snp_vcf_dir
    paste_indel_label_command = "paste %s/unique_indel_positions_file " % args.filter2_only_snp_vcf_dir
    paste_label_command_exclude_outgroup = "paste %s/unique_positions_file " % args.filter2_only_snp_vcf_dir
    paste_indel_label_command_exclude_outgroup = "paste %s/unique_indel_positions_file " % args.filter2_only_snp_vcf_dir

    for filename_base in final_merge_anno_file.samples:
        if re.search('R1_001_final.fastq.gz', filename_base):
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
        elif re.search('_R1.fastq.gz', filename_base):
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
            # Changed on 03/15/2019
        elif re.search('R1.fastq.gz', filename_base):
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
            # Changed on 03/15/2019
            # first_part = re.sub("_S.*", "", first_part)
        elif re.search('1_combine.fastq.gz', filename_base):
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
        elif re.search('1_sequence.fastq.gz', filename_base):
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
        elif re.search('_forward.fastq.gz', filename_base):
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
        elif re.search('R1_001.fastq.gz', filename_base):
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
        elif re.search('_1.fastq.gz', filename_base):
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
        elif re.search('.1.fastq.gz', filename_base):
            second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
            first_part_split = filename_base.split('.1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
        # print filename_base
        #first_part = re.sub("_S.*_", "", first_part)
        # sample_label_file = "%s/%s_filter2_final.vcf_no_proximate_snp.vcf_positions_label" % (
        # args.filter2_only_snp_vcf_dir, first_part)
        # sample_indel_label_file = "%s/%s_filter2_indel_final.vcf_indel_positions_label" % (
        #     args.filter2_only_snp_vcf_dir, first_part)
        
        sample_label_file = "%s/%s_filter2_final.vcf_no_proximate_snp.vcf_positions_label" % ((args.filter2_only_snp_vcf_dir).replace('core_temp_dir', '%s/%s_vcf_results' % (first_part, first_part)), first_part)
        sample_indel_label_file = "%s/%s_filter2_indel_final.vcf_indel_positions_label" % ((args.filter2_only_snp_vcf_dir).replace('core_temp_dir', '%s/%s_vcf_results' % (first_part, first_part)), first_part)
        
        paste_label_command = paste_label_command + sample_label_file + " "
        paste_indel_label_command = paste_indel_label_command + sample_indel_label_file + " "
        if args.outgroup:
            if outgroup not in sample_label_file:
                paste_label_command_exclude_outgroup = paste_label_command_exclude_outgroup + sample_label_file + " "
                paste_indel_label_command_exclude_outgroup = paste_indel_label_command_exclude_outgroup + sample_indel_label_file + " "

    paste_label_command = paste_label_command + " > %s/All_label_final_ordered.txt" % args.filter2_only_snp_vcf_dir
    paste_indel_label_command = paste_indel_label_command + " > %s/All_indel_label_final_ordered.txt" % args.filter2_only_snp_vcf_dir
    sort_ordered_label_cmd = "sort -n -k1,1 %s/All_label_final_ordered.txt > %s/All_label_final_ordered_sorted.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    sort_ordered_indel_label_cmd = "sort -n -k1,1 %s/All_indel_label_final_ordered.txt > %s/All_indel_label_final_ordered_sorted.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    if args.outgroup:
        paste_label_command_exclude_outgroup = paste_label_command_exclude_outgroup + " > %s/All_label_final_ordered_exclude_outgroup.txt" % args.filter2_only_snp_vcf_dir
        paste_indel_label_command_exclude_outgroup = paste_indel_label_command_exclude_outgroup + " > %s/All_indel_label_final_ordered_exclude_outgroup.txt" % args.filter2_only_snp_vcf_dir
        sort_ordered_label_cmd_exclude_outgroup = "sort -n -k1,1 %s/All_label_final_ordered_exclude_outgroup.txt > %s/All_label_final_ordered_exclude_outgroup_sorted.txt" % (
            args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
        sort_ordered_indel_label_cmd_exclude_outgroup = "sort -n -k1,1 %s/All_indel_label_final_ordered_exclude_outgroup.txt > %s/All_indel_label_final_ordered_exclude_outgroup_sorted.txt" % (
            args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    with open('%s/All_label_final_ordered.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
        outfile.write(paste_label_command + '\n')
        outfile.write(sort_ordered_label_cmd + '\n')
        outfile.write(paste_indel_label_command + '\n')
        outfile.write(sort_ordered_indel_label_cmd + '\n')
    outfile.close()

    os.system("bash %s/All_label_final_ordered.sh" % args.filter2_only_snp_vcf_dir)

    if args.outgroup:
        # Just in case if os.system past commands doesn't work
        with open('%s/All_label_final_ordered_exclude_outgroup.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
            outfile.write(paste_label_command_exclude_outgroup + '\n')
            outfile.write(sort_ordered_label_cmd_exclude_outgroup + '\n')
            outfile.write(paste_indel_label_command_exclude_outgroup + '\n')
            outfile.write(sort_ordered_indel_label_cmd_exclude_outgroup + '\n')
        outfile.close()

        # Changed: Uncomment this
        os.system("bash %s/All_label_final_ordered_exclude_outgroup.sh" % args.filter2_only_snp_vcf_dir)

    """ End: Prepare a All_indel_label_final_ordered_sorted.txt file with sorted unique variant positions. """

    """ Generate a position_label and position_indel_label dictionary that will contain information about each unique variant position that passed variant filters in any sample and reasons for being filtered out in any sample """
    position_label = OrderedDict()
    with open("%s/All_label_final_ordered_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
        keep_logging(
            '- Loading %s/All_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
            '- Loading %s/All_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
            logger, 'info')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            position_label[row[0]] = ','.join(row[1:])
    csv_file.close()

    # #Commented for debugging
    position_indel_label = OrderedDict()
    with open("%s/All_indel_label_final_ordered_sorted.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
        keep_logging(
            '- Loading %s/All_indel_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
            '- Loading %s/All_indel_label_final_ordered_sorted.txt' % args.filter2_only_snp_vcf_dir,
            logger, 'info')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if row[0] not in position_label.keys():
                position_indel_label[row[0]] = ','.join(row[1:])
            else:
                position_indel_label[row[0]] = ','.join(row[1:])
                # keep_logging('Warning: position %s already present as a SNP' % row[0],
                #              'Warning: position %s already present as a SNP' % row[0], logger, 'info')
    csv_file.close()

    # All_label_final_ordered_sorted = pd.read_csv("%s/All_label_final_ordered_sorted.txt" % args.filter2_only_snp_vcf_dir, sep='\t', header=None)
    # All_indel_label_final_ordered_sorted = pd.read_csv("%s/All_indel_label_final_ordered_sorted.txt" % args.filter2_only_snp_vcf_dir, sep='\t', header=None)

    # position_label = All_label_final_ordered_sorted.set_index(0).T.to_dict('list')
    # position_indel_label = All_indel_label_final_ordered_sorted.set_index(0).T.to_dict('list')

    """ End: Generate a position_label and position_indel_label dictionary """
    method_time_taken = datetime.now() - method_start_time

    # keep_logging('- Time taken to complete the generate_position_label_dict method: {}'.format(method_time_taken),
    #              '- Time taken to complete the generate_position_label_dict method: {}'.format(method_time_taken), logger, 'info')
    return position_label, position_indel_label

def get_low_fq_mq_positions(position_label):
    """ Generate mask_fq_mq_positions array with positions where a variant was filtered because of LowFQ or LowMQ """
    mask_fq_mq_positions = []
    mask_fq_mq_positions_outgroup_specific = []
    if args.outgroup:
        position_label_exclude_outgroup = OrderedDict()
        with open("%s/All_label_final_ordered_exclude_outgroup_sorted.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            # keep_logging(
            #     '- Reading All label positions file: %s/All_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     '- Reading All label positions file: %s/All_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                position_label_exclude_outgroup[row[0]] = ','.join(row[1:])
        csv_file.close()

        # Commented for debugging
        position_indel_label_exclude_outgroup = OrderedDict()
        with open("%s/All_indel_label_final_ordered_exclude_outgroup_sorted.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            # keep_logging(
            #     '- Reading All label Indel positions file: %s/All_indel_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     '- Reading All label Indel positions file: %s/All_indel_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                if row[0] not in position_label_exclude_outgroup.keys():
                    position_indel_label_exclude_outgroup[row[0]] = ','.join(row[1:])
                else:
                    position_indel_label_exclude_outgroup[row[0]] = ','.join(row[1:])
                    keep_logging('Warning: position %s already present as a SNP' % row[0],
                                 'Warning: position %s already present as a SNP' % row[0], logger, 'info')
        csv_file.close()

        for key in position_label_exclude_outgroup.keys():
            label_sep_array = position_label_exclude_outgroup[key].split(',')
            for i in label_sep_array:
                if "LowFQ" in str(i):
                    if key not in mask_fq_mq_positions:
                        if int(key) not in outgroup_specific_positions:
                            mask_fq_mq_positions.append(key)
                        elif int(key) in outgroup_specific_positions:
                            mask_fq_mq_positions_outgroup_specific.append(key)
                if i == "HighFQ":
                    if key not in mask_fq_mq_positions:
                        if int(key) not in outgroup_specific_positions:
                            mask_fq_mq_positions.append(key)
                        elif int(key) in outgroup_specific_positions:
                            mask_fq_mq_positions_outgroup_specific.append(key)
    else:
        for key in position_label.keys():
            label_sep_array = position_label[key].split(',')
            for i in label_sep_array:
                if "LowFQ" in str(i):
                    if key not in mask_fq_mq_positions:
                        mask_fq_mq_positions.append(key)
                if i == "HighFQ":
                    if key not in mask_fq_mq_positions:
                        mask_fq_mq_positions.append(key)

    fp = open("%s/mask_fq_mq_positions.txt" % (args.filter2_only_snp_vcf_dir), 'w+')
    for i in mask_fq_mq_positions:
        fp.write(i + '\n')
    fp.close()

    fp = open("%s/mask_fq_mq_positions_outgroup_specific.txt" % (args.filter2_only_snp_vcf_dir), 'w+')
    for i in mask_fq_mq_positions_outgroup_specific:
        fp.write(i + '\n')
    fp.close()

    print "Length of mask_fq_mq_positions:%s" % len(mask_fq_mq_positions)
    print "Length of mask_fq_mq_positions specific to outgroup:%s" % len(mask_fq_mq_positions_outgroup_specific)

    """ End: Generate mask_fq_mq_positions array """
    return mask_fq_mq_positions, mask_fq_mq_positions_outgroup_specific

def get_low_fq_mq_positions_indel(position_label, position_indel_label):
    """ Generate mask_fq_mq_positions array with positions where a variant was filtered because of LowFQ or LowMQ"""
    mask_fq_mq_positions = []
    mask_fq_mq_positions_outgroup_specific = []

    if args.outgroup:
        position_label_exclude_outgroup = OrderedDict()
        with open("%s/All_label_final_ordered_exclude_outgroup_sorted.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            # keep_logging(
            #     'Reading All label positions file: %s/All_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     'Reading All label positions file: %s/All_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                position_label_exclude_outgroup[row[0]] = ','.join(row[1:])
        csv_file.close()

        position_indel_label_exclude_outgroup = OrderedDict()
        with open("%s/All_indel_label_final_ordered_exclude_outgroup_sorted.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            # keep_logging(
            #     'Reading All label positions file: %s/All_indel_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     'Reading All label positions file: %s/All_indel_label_final_ordered_exclude_outgroup_sorted.txt' % args.filter2_only_snp_vcf_dir,
            #     logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                if row[0] not in position_label_exclude_outgroup.keys():
                    position_indel_label_exclude_outgroup[row[0]] = ','.join(row[1:])
                else:
                    position_indel_label_exclude_outgroup[row[0]] = ','.join(row[1:])
                    keep_logging('Warning: position %s already present as a SNP' % row[0],
                                 'Warning: position %s already present as a SNP' % row[0], logger, 'info')
        csv_file.close()
        for key in position_label_exclude_outgroup.keys():
            label_sep_array = position_label_exclude_outgroup[key].split(',')
            for i in label_sep_array:
                if "LowFQ" in str(i):
                    if key not in mask_fq_mq_positions:
                        if int(key) not in outgroup_specific_positions:
                            mask_fq_mq_positions.append(key)
                        elif int(key) in outgroup_specific_positions:
                            mask_fq_mq_positions_outgroup_specific.append(key)
                if i == "HighFQ":
                    if key not in mask_fq_mq_positions:
                        if int(key) not in outgroup_specific_positions:
                            mask_fq_mq_positions.append(key)
                        elif int(key) in outgroup_specific_positions:
                            mask_fq_mq_positions_outgroup_specific.append(key)
    else:
        for key in position_label.keys():
            label_sep_array = position_label[key].split(',')
            for i in label_sep_array:
                if "LowFQ" in str(i):
                    if key not in mask_fq_mq_positions:
                        mask_fq_mq_positions.append(key)
                if i == "HighFQ":
                    if key not in mask_fq_mq_positions:
                        mask_fq_mq_positions.append(key)

    # print "Length of Indel mask_fq_mq_positions:%s" % len(mask_fq_mq_positions)
    # print "Length of Indel mask_fq_mq_positions specific to outgroup:%s" % len(mask_fq_mq_positions_outgroup_specific)
    return mask_fq_mq_positions, mask_fq_mq_positions_outgroup_specific

def generate_SNP_matrix(final_merge_anno_file, functional_filter_pos_array, phage_positions, repetitive_positions,
                        mask_positions, position_label, core_positions, snp_var_ann_dict):
    method_start_time = datetime.now()
    """ Prepare SNP/Indel Matrix print strings and add matrix row information subsequently """

    header_print_string = "Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product"
    for sample in final_merge_anno_file.samples:
        header_print_string = header_print_string + "\t" + sample
    header_print_string = header_print_string + "\n"

    """ Main: Generate SNP Matrix """

    """ Open Matrix files to write strings """
    fp_code = open("%s/SNP_matrix_code.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_code_unmasked = open("%s/SNP_matrix_code_unmasked.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    # fp_allele = open("%s/SNP_matrix_allele_outdated.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele_new = open("%s/SNP_matrix_allele_new.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele_phage = open("%s/SNP_matrix_allele_unmasked.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_code.write(header_print_string)
    fp_code_unmasked.write(header_print_string)
    #fp_allele.write(header_print_string)
    fp_allele_new.write(header_print_string)
    fp_allele_phage.write(header_print_string)
    
    """ Parse variant positions from the loaded cyvcf VCF object/ merged gatk multivcf file and generate matrix row information """
    for variants in VCF("%s/Final_vcf_gatk_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir):
        # Initiate print_string variable to add matrix row information.
        # print_string generator no. 1
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

        # Initiate variant code string where the code means:
        # REF allele = 0,
        # core = 1,
        # Filtered = 2,
        # unmapped = -1,
        # True but non-core = 3

        code_string = position_label[str(variants.POS)]
        
        # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
        code_string = code_string.replace('reference_allele_indel_proximate', '4')
        code_string = code_string.replace('reference_allele', '0')
        code_string = code_string.replace('reference_unmapped_position', '-1')
        # Changing LowFQ code from 2 to -3
        # Changing HighFQ but LowMQ code from 2 to -4
        code_string = code_string.replace('LowFQ_QUAL_DP_proximate_SNP', '-3')
        code_string = code_string.replace('LowFQ_DP_QUAL_proximate_SNP', '-3')
        code_string = code_string.replace('LowFQ_QUAL_proximate_SNP', '-3')
        code_string = code_string.replace('LowFQ_DP_proximate_SNP', '-3')
        code_string = code_string.replace('LowFQ_proximate_SNP', '-3')
        code_string = code_string.replace('LowFQ_QUAL_DP', '-3')
        code_string = code_string.replace('LowFQ_DP_QUAL', '-3')
        code_string = code_string.replace('LowFQ_QUAL', '-3')
        code_string = code_string.replace('LowFQ_DP', '-3')
        code_string = code_string.replace('HighFQ_QUAL_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_DP_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_proximate_SNP', '2')
        code_string = code_string.replace('HighFQ_QUAL_DP', '2')
        code_string = code_string.replace('HighFQ_DP_QUAL', '2')
        code_string = code_string.replace('HighFQ_QUAL', '2')
        code_string = code_string.replace('HighFQ_DP', '2')
        code_string = code_string.replace('LowFQ', '-3')
        code_string = code_string.replace('HighFQ', '-4')
        code_string_unmasked =code_string

        if str(variants.POS) in core_positions:
            code_string = code_string.replace('VARIANT', '1')
            code_string_unmasked = code_string_unmasked.replace('VARIANT', '1')
        # Adding functional class status code to SNP matrix: 2018-07-24
        elif str(variants.POS) in functional_filter_pos_array:
            # Changing Functional class filter code to -2 from 2: 2018-12-04
            code_string = code_string.replace('VARIANT', '-2')
            code_string_unmasked = code_string_unmasked.replace('VARIANT', '3')
            # code_string_unmasked = code_string_unmasked.replace('VARIANT', 'PHAGE')
        else:
            code_string = code_string.replace('VARIANT', '3')
            code_string_unmasked = code_string_unmasked.replace('VARIANT', '3')

        # Annotation Bug fix 2
        # Changing SNP type: Date 28/05/2019
        if variants.POS in snp_var_ann_dict.keys():
            if snp_var_ann_dict[variants.POS] is not None:
                if "protein_coding" in set(
                        snp_var_ann_dict[variants.POS].split('|')) and "intergenic_region" not in set(
                        snp_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Coding SNP"
                elif "protein_coding" in set(snp_var_ann_dict[variants.POS].split('|')) and "intergenic_region" in set(
                        snp_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Coding and Non-coding SNP"
                elif "protein_coding" not in set(
                        snp_var_ann_dict[variants.POS].split('|')) and "intergenic_region" in set(
                        snp_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Non-Coding SNP"
                elif "protein_coding" not in set(
                        snp_var_ann_dict[variants.POS].split('|')) and "intragenic_variant" in set(
                        snp_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Non-Coding SNP"
                else:
                    # print set((snp_var_ann_dict[variants.POS].split('|')))
                    snp_type = "No_protein_coding/intergenic_region_field_in_ANN SNP"
            # print snp_type
        else:
            keep_logging(
                'Warning: position %s not found in snp_var_ann_dict dictionary. Assigning Not found as SNP type.' % variants.POS,
                'Warning: position %s not found in snp_var_ann_dict dictionary. Assigning Not found as SNP type.' % variants.POS,
                logger, 'info')
            # print set((snp_var_ann_dict[variants.POS].split('|')))
            snp_type = "Not Found in Annotated VCF file"

        # print_string generator no. 2
        print_string = print_string + snp_type + " at %s > " % str(variants.POS) + str(
            ",".join(variants.ALT)) + " functional=%s" % functional_field

        # Get ANN field from variant INFO column and save it as an array. Split and Go through each elements, add bells and whistles
        if variants.INFO.get('ANN'):

            ann_array = (variants.INFO.get('ANN')).split(',')

            # Generate tag string before generating ann_string
            if len(ann_array) > 1:
                # print variants.INFO.get('ANN')
                # print list(set(ann_array))
                tag_list = []

                for i_again in set(snp_var_ann_dict[variants.POS].split(',')):
                    i_split_again = i_again.split('|')

                    if "-" not in i_split_again[4]:
                        if i_split_again[4] not in tag_list:
                            tag_list.append(i_split_again[4])

                    else:
                        split_tags = i_split_again[4].split('-')
                        for splittagsindividual in split_tags:
                            if splittagsindividual not in tag_list:
                                tag_list.append(splittagsindividual)

                if len(tag_list) == 1:
                    tag = tag_list[0]

                elif len(tag_list) == 2:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1])

                elif len(tag_list) == 3:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1]) + "-" + str(tag_list[2])
                    print "Error: More than 2 Locus Tags were found at %s - %s" % (variants.POS, tag_list)
                    # exit()
                elif len(tag_list) > 3:
                    print "Error: More than 3 Locus Tags were found at %s - %s" % (variants.POS, tag_list)
                    tag = str(tag_list[0]) + "-" + str(tag_list[1]) + "-" + str(tag_list[2]) + "-" + str(tag_list[3])
                tag = tag.replace('CHR_START-', '')
                tag = tag.replace('-CHR_END', '')
            else:
                for i in list(set(ann_array)):
                    i_split = i.split('|')
                    tag = str(i_split[4]).replace('CHR_START-', '')
                    tag = str(tag).replace('-CHR_END', '')

            # Generate ann_string variable
            ann_string = ";"
            for i in list(set(ann_array)):
                i_split = i.split('|')
                # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"

                # MOve this tag before this for loop because of multiple tags associated.
                # tag = str(i_split[4]).replace('CHR_START-', '')
                # tag = str(tag).replace('-CHR_END', '')

                if "-" in tag:
                    # print tag
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
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags, extra_tags_prot]) + ";"
                    # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                    ann_string = ann_string + '|'.join(
                        [i_split[0], i_split[1], i_split[2], i_split[4], i_split[9], i_split[10], i_split[11],
                         i_split[13], extra_tags, extra_tags_prot]) + ";"

                # Changing SNP type: Date 28/05/2019
                elif tag == "":
                    print "ERROR: Issues with this locus tag. Check this tag in genbank file"
                    print list(set(ann_array))
                    # Adding this so that Ann string is not empty: 30/05/2019
                    if tag in locus_tag_to_gene_name.keys() and tag in locus_tag_to_product.keys():
                        extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                    else:
                        print "tag key not found: %s" % tag
                        extra_tags = "NULL" + "|" + "NULL"
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                    # Added 2019-31-05
                    if "ERROR_OUT_OF_CHROMOSOME_RANGE" in i:
                        ann_string = ann_string + '|'.join(
                            [i_split[0], "intergenic_region", i_split[2], "ERROR_OUT_OF_CHROMOSOME_RANGE", i_split[9],
                             i_split[10], i_split[11],
                             i_split[13], extra_tags]) + ";"
                        # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                        # No changes here
                    else:
                        # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                        # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                        ann_string = ann_string + '|'.join(
                            [i_split[0], i_split[1], i_split[2], i_split[4], i_split[9], i_split[10], i_split[11],
                             i_split[13], extra_tags]) + ";"
                # Changing SNP type: Date 28/05/2019
                else:
                    if tag in locus_tag_to_gene_name.keys() and tag in locus_tag_to_product.keys():
                        extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                    else:
                        print "tag key not found: %s" % tag
                        extra_tags = "NULL" + "|" + "NULL"
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                    # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                    ann_string = ann_string + '|'.join(
                        [i_split[0], i_split[1], i_split[2], i_split[4], i_split[9], i_split[10], i_split[11],
                         i_split[13], extra_tags]) + ";"
        # Annotation Bug fix 4
        # Changing SNP type: Date 28/05/2019
        # Working/Testing
        else:
            if len(variants.ALT) > 1 and snp_var_ann_dict[variants.POS]:
                # print variants.ALT
                # print ';'.join(set(snp_var_ann_dict[variants.POS].split(',')))

                ann_string = ";%s" % ';'.join(set(snp_var_ann_dict[variants.POS].split(',')))
                # Get Tag here; Multiple tag names.
                tag_list = []

                for i in set(snp_var_ann_dict[variants.POS].split(',')):
                    i_split = i.split('|')
                    if i_split[4] not in tag_list:
                        tag_list.append(i_split[4])
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # Changing >1 to ==2
                if len(tag_list) == 1:
                    tag = tag_list[0]
                elif len(tag_list) == 2:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1])
                elif len(tag_list) == 3:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1]) + "-" + str(tag_list[2])
                    print "Error: More than 2 Locus Tags were found at %s - %s" % (variants.POS, tag_list)
                    # exit()
                elif len(tag_list) > 3:
                    print "Error: More than 3 Locus Tags were found at %s - %s" % (variants.POS, tag_list)
                    tag = str(tag_list[0]) + "-" + str(tag_list[1]) + "-" + str(tag_list[2]) + "-" + str(tag_list[3])

            else:
                ann_string = ";None"

        # Annotation Bug fix 5
        # Changing SNP type: Date 28/05/2019
        ann_string = ann_string.replace('ERROR_OUT_OF_CHROMOSOME_RANGE', '%s-%s' % (
        locus_tag_to_gene_name[last_locus_tag], locus_tag_to_gene_name[first_locus_tag]))
        ann_string = ann_string.replace('CHR_END', '%s' % locus_tag_to_gene_name[first_locus_tag])

        # SNP Matrix Bug
        # No changes here: 28/05/2019
        ann_string_split = ann_string.split(';')
        # print len(ann_string_split)
        if len(ann_string_split) == 3:
            first_allele_ann_string_split = ann_string_split[1].split('|')
            second_allele_ann_string_split = ann_string_split[2].split('|')
            if len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) == 10:
                ann_string = ann_string
            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) == 10:
                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + "|" + first_allele_ann_string_split[15]
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[2] + "|" + first_allele_ann_string_split[4] + "|" + first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[10] + "|" + first_allele_ann_string_split[11] + "|" + first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[
                                                  2] + "|" + first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[
                                                  10] + "|" + first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                ann_string = new_first_allele_ann_string + str(ann_string_split[2])

            elif len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) > 10:

                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + "|" + second_allele_ann_string_split[15]
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[1] + "|" + second_allele_ann_string_split[2] + "|" + \
                # second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[9] + "|" + \
                # second_allele_ann_string_split[10] + "|" + second_allele_ann_string_split[11] + "|" + \
                # second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[
                    1] + "|" + second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[
                                                   9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                ann_string = str(ann_string_split[1]) + new_second_allele_ann_string
            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) > 10:

                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + "|" + first_allele_ann_string_split[15]
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[2] + "|" + first_allele_ann_string_split[4] + "|" + first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[10] + "|" + first_allele_ann_string_split[11] + "|" + first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[
                                                  2] + "|" + first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[
                                                  10] + "|" + first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + "|" + second_allele_ann_string_split[15]
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[1] + "|" + second_allele_ann_string_split[2] + "|" + \
                # second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[9] + "|" + \
                # second_allele_ann_string_split[10] + "|" + second_allele_ann_string_split[11] + "|" + \
                # second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[
                    1] + "|" + second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[
                                                   9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                ann_string = new_first_allele_ann_string + new_second_allele_ann_string

        # Changed from >3 to ==4 Issue #29 26-03-2020
        if len(ann_string_split) > 3:
            first_allele_ann_string_split = ann_string_split[1].split('|')
            second_allele_ann_string_split = ann_string_split[2].split('|')
            third_allele_ann_string_split = ann_string_split[3].split('|')
            if len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) == 10 and len(
                    third_allele_ann_string_split) == 10:
                ann_string = ann_string

            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) == 10 and len(
                    third_allele_ann_string_split) == 10:
                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + "|" + first_allele_ann_string_split[15]
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[2] + "|" + first_allele_ann_string_split[4] + "|" + first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[10] + "|" + first_allele_ann_string_split[11] + "|" + first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[
                                                  2] + "|" + first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[
                                                  10] + "|" + first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                ann_string = new_first_allele_ann_string + str(ann_string_split[2]) + str(ann_string_split[3])

            elif len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) > 10 and len(
                    third_allele_ann_string_split) == 10:

                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + "|" + second_allele_ann_string_split[15]
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[1] + "|" + second_allele_ann_string_split[2] + "|" + \
                # second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[9] + "|" + \
                # second_allele_ann_string_split[10] + "|" + second_allele_ann_string_split[11] + "|" + \
                # second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[
                    1] + "|" + second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[
                                                   9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                ann_string = str(ann_string_split[1]) + new_second_allele_ann_string + str(ann_string_split[3])

            elif len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) == 10 and len(
                    third_allele_ann_string_split) > 10:

                if third_allele_ann_string_split[14] == "" and third_allele_ann_string_split[15] == "":
                    prod = third_allele_ann_string_split[3] + third_allele_ann_string_split[15]
                else:
                    prod = third_allele_ann_string_split[14] + "|" + third_allele_ann_string_split[15]
                # new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + third_allele_ann_string_split[1] + "|" + third_allele_ann_string_split[2] + "|" + \
                #                               third_allele_ann_string_split[4] + "|" + third_allele_ann_string_split[9] + "|" + \
                #                               third_allele_ann_string_split[10] + "|" + third_allele_ann_string_split[11] + "|" + \
                #                               third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + third_allele_ann_string_split[
                    1] + "|" + third_allele_ann_string_split[2] + "|" + \
                                              third_allele_ann_string_split[4] + "|" + third_allele_ann_string_split[
                                                  9] + "|" + \
                                              third_allele_ann_string_split[10] + "|" + third_allele_ann_string_split[
                                                  11] + "|" + \
                                              third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                ann_string = str(ann_string_split[1]) + str(ann_string_split[2]) + new_third_allele_ann_string

            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) > 10 and len(
                    third_allele_ann_string_split) > 10:
                # print ann_string
                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + "|" + first_allele_ann_string_split[15]
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[2] + "|" + first_allele_ann_string_split[4] + "|" + first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[10] + "|" + first_allele_ann_string_split[11] + "|" + first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + first_allele_ann_string_split[
                                                  2] + "|" + first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + first_allele_ann_string_split[
                                                  10] + "|" + first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + "|" + second_allele_ann_string_split[15]
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[1] + "|" + second_allele_ann_string_split[2] + "|" + \
                # second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[9] + "|" + \
                # second_allele_ann_string_split[10] + "|" + second_allele_ann_string_split[11] + "|" + \
                # second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + second_allele_ann_string_split[
                    1] + "|" + second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + second_allele_ann_string_split[
                                                   9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                if third_allele_ann_string_split[14] == "" and third_allele_ann_string_split[15] == "":
                    prod = third_allele_ann_string_split[3] + third_allele_ann_string_split[15]
                else:
                    prod = third_allele_ann_string_split[14] + "|" + third_allele_ann_string_split[15]
                # new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + third_allele_ann_string_split[1] + "|" + third_allele_ann_string_split[2] + "|" + \
                #                               third_allele_ann_string_split[4] + "|" + third_allele_ann_string_split[9] + "|" + \
                #                               third_allele_ann_string_split[10] + "|" + third_allele_ann_string_split[11] + "|" + \
                #                               third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + third_allele_ann_string_split[
                    1] + "|" + third_allele_ann_string_split[2] + "|" + \
                                              third_allele_ann_string_split[4] + "|" + third_allele_ann_string_split[
                                                  9] + "|" + \
                                              third_allele_ann_string_split[10] + "|" + third_allele_ann_string_split[
                                                  11] + "|" + \
                                              third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                ann_string = new_first_allele_ann_string + new_second_allele_ann_string + new_third_allele_ann_string

        # Added this extra check Issue #29 26-03-2020
        if len(ann_string_split) > 4:
            # print ann_string_split
            # print variants.POS
            new_allele_string_array = []
            for i_ann_string_split in ann_string_split[1:]:
                if len(i_ann_string_split.split('|')) == 10:
                    ann_string = ann_string
                elif len(i_ann_string_split.split('|')) > 10:
                    ann_string = ";"
                    i_ann_string_split_array = i_ann_string_split.split('|')
                    if i_ann_string_split_array[14] == "" and i_ann_string_split_array[15] == "":
                        prod = i_ann_string_split_array[3] + i_ann_string_split_array[15]
                    else:
                        prod = i_ann_string_split_array[14] + i_ann_string_split_array[15]

                    # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                    # new_allele_string = i_ann_string_split_array[0] + "|" + i_ann_string_split_array[1] + "|" + \
                    #                     i_ann_string_split_array[2] + "|" + \
                    #                     i_ann_string_split_array[4] + "|" + \
                    #                     i_ann_string_split_array[9] + "|" + \
                    #                     i_ann_string_split_array[10] + "|" + \
                    #                     i_ann_string_split_array[11] + "|" + \
                    #                     i_ann_string_split_array[13] + "|" + prod + "|" + prod
                    new_allele_string = i_ann_string_split_array[0] + "|" + i_ann_string_split_array[1] + "|" + \
                                        i_ann_string_split_array[2] + "|" + \
                                        i_ann_string_split_array[4] + "|" + \
                                        i_ann_string_split_array[9] + "|" + \
                                        i_ann_string_split_array[10] + "|" + \
                                        i_ann_string_split_array[11] + "|" + \
                                        i_ann_string_split_array[13] + "|" + prod + "|" + prod
                    new_allele_string_array.append(new_allele_string)

                    ann_string = ann_string + ";".join(new_allele_string_array)

        # print_string generator no. 3

        # Annotation Bug fix 6
        # Changing Strandness string: Date 28/05/2019
        # Each Locus ID with a strand information

        strandness = " Strand Information: "
        if "-" in tag:
            tagsplit = tag.split('-')
            for i in tagsplit:
                if i in locus_tag_to_strand.keys():
                    if "," in locus_tag_to_strand[i]:
                        locus_tag_to_strand_split = locus_tag_to_strand[i].split(',')
                        strand = locus_tag_to_strand_split[0]
                    else:
                        strand = locus_tag_to_strand[i]
                    strandness = strandness + i + "=" + strand + "/"
                else:
                    if i == "" or i == "None":
                        strandness = strandness + "NULL=" + "No Strand Information found" + "/"
                    else:
                        strandness = strandness + i + "=" + "No Strand Information found" + "/"
        else:
            if tag in locus_tag_to_strand.keys():
                # strandness = strandness + locus_tag_to_strand[tag]
                if "," in locus_tag_to_strand[tag]:
                    locus_tag_to_strand_split = locus_tag_to_strand[tag].split(',')
                    strand = locus_tag_to_strand_split[0]
                else:
                    strand = locus_tag_to_strand[tag]
                strandness = strandness + tag + "=" + strand
            else:
                if tag == "" or tag == "None":
                    strandness = strandness + "NULL=" + "No Strand Information found"
                else:
                    strandness = strandness + tag + "=" + "No Strand Information found"

        # Annotation Bug fix 7
        # Changing tag equals NULL: Date 30/05/2019
        if tag == "" or tag == "None":
            tag = "NULL"

        ann_string_old = ann_string
        # Test No. of Fields
        ann_string_split = ann_string.split(';')
        for i in ann_string_split:
            ann_string_array_split_columns = i.split('|')
            if len(ann_string_array_split_columns) > 10:
                #print variants.POS
                # ann_string = ann_string.replace('||WARNING_TRANSCRIPT_NO_START_CODON||WARNING_TRANSCRIPT_NO_START_CODON', '|WARNING_TRANSCRIPT_NO_START_CODON|WARNING_TRANSCRIPT_NO_START_CODON')
                ann_string = ";" + str('|'.join(ann_string_array_split_columns[: -2 or None])) + ";"
                #print ann_string
            else:
                ann_string = ann_string_old
                # print ann_string
        ann_string = ann_string.replace(';;', ';')
        print_string = print_string + " locus_tag=" + tag + strandness + ann_string
        print_string_phage = print_string

        """ Go over each genotype for a variant and generate a gt_string variable """
        gt_string = ""
        for gt in variants.gt_bases:
            gt_modified = gt.replace("./.", variants.REF)
            gt_string = gt_string + "," + gt_modified
        gt_string = gt_string.replace("A/A", "A")
        gt_string = gt_string.replace("G/G", "G")
        gt_string = gt_string.replace("C/C", "C")
        gt_string = gt_string.replace("T/T", "T")

        # Extra Check
        if "A/A" in gt_string:
            gt_string = gt_string.replace("A/A", "A")
        elif "G/G" in gt_string:
            gt_string = gt_string.replace("G/G", "G")
        elif "C/C" in gt_string:
            gt_string = gt_string.replace("C/C", "C")
        elif "T/T" in gt_string:
            gt_string = gt_string.replace("T/T", "T")
        # Not recognizing "." as a string. Its also redundant to replace "./." with "." and again with Reference allele
        # elif "." in gt_string:
        #     gt_string = gt_string.replace(".", variants.REF)
        #     print "Go over each genotype for a variant and generate a gt_string variable"
        # else:
        #     print "Warning: Doesn't recognize GT code at %s - %s" % (variants.POS, gt_string)
        #     exit()
        
        
        # print_string generator no. 4
        # Replace various seperators that were used in old matrix. Clean up this block of code
        final_allele_string = print_string + gt_string.replace(',', '\t') + '\n'
        # Replace code at Phage Positions with -2
        if str(variants.POS) in functional_filter_pos_array:
            code_string_array = code_string.split(',')
            for (i, item) in enumerate(code_string_array):
                if item == "0":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "1":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "2":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "3":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "4":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "-1":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "-2":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "-3":
                    code_string_array[i] = "-2"
            for (i, item) in enumerate(code_string_array):
                if item == "-4":
                    code_string_array[i] = "-2"
            code_string = ','.join(code_string_array)

        final_code_string = print_string + "\t" + code_string.replace(',', '\t') + '\n'
        final_code_string_unmasked = print_string + "\t" + code_string_unmasked.replace(',', '\t') + '\n'
        final_allele_string = final_allele_string.replace(',|', '|')

        final_allele_string = final_allele_string.replace(',;,', ':::')
        final_allele_string = final_allele_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(',|', '|')
        final_code_string_unmasked = final_code_string_unmasked.replace(',|', '|')

        final_code_string_unmasked = final_code_string_unmasked.replace(',;,', ':::')
        final_code_string_unmasked = final_code_string_unmasked.replace(';,', ':::')
        final_code_string_unmasked = final_code_string_unmasked.replace(';\t\t', ';\t')
        final_code_string_unmasked = final_code_string_unmasked.replace('\t\t', '\t')

        final_code_string = final_code_string.replace(',;,', ':::')
        final_code_string = final_code_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(';\t\t', ';\t')
        final_code_string = final_code_string.replace('\t\t', '\t')

        final_allele_string = final_allele_string.replace('\t\t', '\t')
        # fp_allele.write(final_allele_string)
        fp_code.write(final_code_string)
        fp_code_unmasked.write(final_code_string_unmasked)

        ntd_string = ""
        ntd_string_phage = ""
        count = 0
        code_string_array = code_string.split(',')
        gt_string_array = gt_string[1:].split(',')

        for i in gt_string_array:
            if str(code_string_array[count]) == "0" or str(code_string_array[count]) == "1" or str(
                    code_string_array[count]) == "3":
                ntd_string = ntd_string + "\t" + str(i)
                ntd_string_phage = ntd_string_phage + "\t" + str(i)
            if code_string_array[count] == "-1":
                ntd_string = ntd_string + "\t" + "-"
                ntd_string_phage = ntd_string_phage + "\t" + "-"
            # Changing Functional class filter code to -2 from 2 and replacing variant allele with N: 2018-12-04
            if str(code_string_array[count]) == "2" or str(code_string_array[count]) == "-2" or str(
                    code_string_array[count]) == "-3" or str(code_string_array[count]) == "-4" or str(code_string_array[count]) == "4":
                ntd_string = ntd_string + "\t" + "N"
            if str(code_string_array[count]) == "2":
                ntd_string_phage = ntd_string_phage + "\t" + "N"
            if str(code_string_array[count]) == "-2":
                ntd_string_phage = ntd_string_phage + "\t" + str(i)
            count += 1

        # Annotation Bug fix 8
        """ Mask Phage positions in SNP_matrix_allele_new.tsv. This is the default matrix. """
        if str(variants.POS) in functional_filter_pos_array:
            ntd_string_array = ntd_string.split('\t')
            # print ntd_string_array
            ntd_string = ""
            for i in ntd_string_array[1:]:
                ntd_string = ntd_string + "\t" + "N"
            ntd_string_array = ntd_string.split('\t')
            # print ntd_string_array

        """ Generate a print_string for each of the matrix - SNP_matrix_allele_new.tsv and SNP_matrix_allele_phage.tsv """
        print_string = print_string + ntd_string + "\n"

        print_string_phage = print_string_phage + ntd_string_phage + "\n"

        """ This is a hardcoded solution. Find the root cause of these strings getting into the print_strint variable """
        print_string.replace(',;,', '\t')
        print_string.replace(';,', '\t')
        fp_allele_new.write(print_string)

        print_string_phage.replace(',;,', '\t')
        print_string_phage.replace(';,', '\t')
        fp_allele_phage.write(print_string_phage)
    
    fp_code.close()
    fp_allele_new.close()
    fp_allele_phage.close()
    fp_code_unmasked.close()
    # fp_allele.close()
    
    
    method_time_taken = datetime.now() - method_start_time

    keep_logging('- Time taken to complete the Generate SNP matrix method: {}'.format(method_time_taken),
                 '- Time taken to complete the Generate SNP matrix method: {}'.format(method_time_taken), logger, 'info')

def generate_Indel_matrix(final_merge_anno_file, functional_filter_pos_array, phage_positions, repetitive_positions,
                          mask_positions, position_indel_label, indel_core_positions, indel_var_ann_dict):
    method_start_time = datetime.now()
    """ Prepare Indel Matrix header print strings and add matrix row information subsequently """
    header_print_string = "Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product"
    for sample in final_merge_anno_file.samples:
        # header_print_string = header_print_string + "," + sample
        header_print_string = header_print_string + "\t" + sample
    header_print_string = header_print_string + "\n"

    fp_code = open("%s/Indel_matrix_code.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_code_unmasked = open("%s/Indel_matrix_code_unmasked.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    fp_allele = open("%s/Indel_matrix_allele.tsv" % args.filter2_only_snp_vcf_dir, 'w+')
    
    fp_code.write(header_print_string)
    fp_code_unmasked.write(header_print_string)
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
        # Extract positions filtered by Indel Proximate filters and assign N instead of reference allele - 2020-05-20
        code_string = code_string.replace('reference_allele_indel_proximate', '4')
        code_string = code_string.replace('reference_allele', '0')
        code_string = code_string.replace('reference_unmapped_position', '-1')
        code_string = code_string.replace('LowAF_QUAL_DP_proximate_SNP', '2')
        code_string = code_string.replace('LowAF_DP_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('LowAF_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('LowAF_DP_proximate_SNP', '2')
        code_string = code_string.replace('LowAF_proximate_SNP', '2')
        code_string = code_string.replace('LowAF_QUAL_DP', '2')
        code_string = code_string.replace('LowAF_DP_QUAL', '2')
        code_string = code_string.replace('LowAF_QUAL', '2')
        code_string = code_string.replace('LowAF_DP', '2')
        code_string = code_string.replace('HighAF_QUAL_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighAF_DP_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighAF_QUAL_proximate_SNP', '2')
        code_string = code_string.replace('HighAF_DP_proximate_SNP', '2')
        code_string = code_string.replace('HighAF_proximate_SNP', '2')
        code_string = code_string.replace('HighAF_QUAL_DP', '2')
        code_string = code_string.replace('HighAF_DP_QUAL', '2')
        code_string = code_string.replace('HighAF_QUAL', '2')
        code_string = code_string.replace('HighAF_DP', '2')
        code_string = code_string.replace('LowAF', '-3')
        code_string = code_string.replace('HighAF', '-4')
        code_string_unmasked =code_string

        if str(variants.POS) in indel_core_positions:
            code_string = code_string.replace('VARIANT', '1')
            code_string_unmasked = code_string_unmasked.replace('VARIANT', '1')
        # Adding functional class status code to SNP matrix: 2018-07-24
        elif str(variants.POS) in functional_filter_pos_array:
            # Changing Functional class filter code to -2 from 2: 2018-12-04
            code_string = code_string.replace('VARIANT', '-2')
            code_string_unmasked = code_string_unmasked.replace('VARIANT', '3')
        else:
            code_string = code_string.replace('VARIANT', '3')
            code_string_unmasked = code_string_unmasked.replace('VARIANT', '3')

        # Changing SNP type: Date 28/05/2019
        # Assign type of snp: coding / non-coding
        if variants.POS in indel_var_ann_dict.keys():
            if indel_var_ann_dict[variants.POS] is not None:
                if "protein_coding" in set(
                        indel_var_ann_dict[variants.POS].split('|')) and "intergenic_region" not in set(
                    indel_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Coding Indel"
                elif "protein_coding" in set(
                        indel_var_ann_dict[variants.POS].split('|')) and "intergenic_region" in set(
                    indel_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Coding and Non-coding Indel"
                elif "protein_coding" not in set(
                        indel_var_ann_dict[variants.POS].split('|')) and "intergenic_region" in set(
                    indel_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Non-Coding Indel"
                elif "protein_coding" not in set(
                        indel_var_ann_dict[variants.POS].split('|')) and "intragenic_variant" in set(
                    indel_var_ann_dict[variants.POS].split('|')):
                    snp_type = "Non-Coding Indel"
                else:
                    # print set((indel_var_ann_dict[variants.POS].split('|')))
                    snp_type = "No_protein_coding/intergenic_region_field_in_ANN SNP"
            # print snp_type
        else:
            keep_logging(
                'Warning: position %s not found in snp_var_ann_dict dictionary. Assigning Not found as SNP type.' % variants.POS,
                'Warning: position %s not found in snp_var_ann_dict dictionary. Assigning Not found as SNP type.' % variants.POS,
                logger, 'info')
            # print set((indel_var_ann_dict[variants.POS].split('|')))
            snp_type = "Not Found in Annotated VCF file"

        print_string = print_string + snp_type + " at %s > " % str(variants.POS) + str(
            ",".join(variants.ALT)) + " functional=%s" % functional_field

        # Get ANN field from variant INFO column and save it as an array. Split and Go through each elements, add bells and whistles
        if variants.INFO.get('ANN'):

            ann_array = (variants.INFO.get('ANN')).split(',')

            # Generate tag string before generating ann_string
            if len(ann_array) > 1:
                # print variants.INFO.get('ANN')
                # print list(set(ann_array))
                tag_list = []

                for i_again in set(indel_var_ann_dict[variants.POS].split(',')):
                    i_split_again = i_again.split('|')

                    if "-" not in i_split_again[4]:
                        if i_split_again[4] not in tag_list:
                            tag_list.append(i_split_again[4])

                    else:
                        split_tags = i_split_again[4].split('-')
                        for splittagsindividual in split_tags:
                            if splittagsindividual not in tag_list:
                                tag_list.append(splittagsindividual)

                if len(tag_list) == 1:
                    tag = tag_list[0]

                elif len(tag_list) == 2:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1])

                elif len(tag_list) == 3:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1]) + "-" + str(tag_list[2])
                    print "Error: More than two locus tags were found at %s - %s" % (variants.POS, tag_list)
                    
                elif len(tag_list) > 3:
                    print "Error: More than three locus tags were found at %s - %s" % (variants.POS, tag_list)
                    
                    exit()
                tag = tag.replace('CHR_START-', '')
                tag = tag.replace('-CHR_END', '')
            else:
                for i in list(set(ann_array)):
                    i_split = i.split('|')
                    tag = str(i_split[4]).replace('CHR_START-', '')
                    tag = str(tag).replace('-CHR_END', '')

            # Generating Annotation string
            ann_string = ";"
            for i in list(set(ann_array)):
                i_split = i.split('|')

                # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13]]) + ";"

                # MOve this tag before this for loop because of multiple tags associated.
                # tag = str(i_split[4]).replace('CHR_START-', '')
                # tag = str(tag).replace('-CHR_END', '')

                if "-" in tag:
                    # print tag
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
                    # ann_string = ann_string + '|'.join(
                    #     [i_split[0], i_split[1], i_split[2], i_split[3], i_split[9], i_split[10], i_split[11],
                    #      i_split[13], extra_tags, extra_tags_prot]) + ";"
                    # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                    ann_string = ann_string + '|'.join(
                        [i_split[0], i_split[1], i_split[2], i_split[4], i_split[9], i_split[10], i_split[11],
                         i_split[13], extra_tags, extra_tags_prot]) + ";"

                # Changing SNP type: Date 28/05/2019
                elif tag == "":
                    print "ERROR: Issues with this locus tag. Check this tag in genbank file"
                    print list(set(ann_array))
                    # Adding this so that Ann string is not empty: 30/05/2019
                    if tag in locus_tag_to_gene_name.keys() and tag in locus_tag_to_product.keys():
                        extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                    else:
                        print "tag key not found: %s" % tag
                        extra_tags = "NULL" + "|" + "NULL"
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                    # Added 2019-31-05
                    if "ERROR_OUT_OF_CHROMOSOME_RANGE" in i:
                        ann_string = ann_string + '|'.join(
                            [i_split[0], "intergenic_region", i_split[2], "ERROR_OUT_OF_CHROMOSOME_RANGE", i_split[9],
                             i_split[10], i_split[11],
                             i_split[13], extra_tags]) + ";"
                    else:
                        # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                        # ann_string = ann_string + '|'.join(
                        #     [i_split[0], i_split[1], i_split[2], i_split[3], i_split[9], i_split[10], i_split[11],
                        #      i_split[13], extra_tags]) + ";"
                        ann_string = ann_string + '|'.join(
                            [i_split[0], i_split[1], i_split[2], i_split[4], i_split[9], i_split[10], i_split[11],
                             i_split[13], extra_tags]) + ";"

                # Changing SNP type: Date 28/05/2019
                else:
                    if tag in locus_tag_to_gene_name.keys() and tag in locus_tag_to_product.keys():
                        extra_tags = str(locus_tag_to_gene_name[tag]) + "|" + str(locus_tag_to_product[tag])
                    else:
                        print "tag key not found: %s" % tag
                        extra_tags = "NULL" + "|" + "NULL"
                    # ann_string = ann_string + '|'.join([i_split[0],i_split[1],i_split[2],i_split[3],i_split[9], i_split[10], i_split[11], i_split[13], extra_tags]) + ";"
                    # ann_string = ann_string + '|'.join(
                    #     [i_split[0], i_split[1], i_split[2], i_split[3], i_split[9], i_split[10], i_split[11],
                    #      i_split[13], extra_tags]) + ";"
                    # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                    ann_string = ann_string + '|'.join(
                        [i_split[0], i_split[1], i_split[2], i_split[4], i_split[9], i_split[10], i_split[11],
                         i_split[13], extra_tags]) + ";"


        # Changing SNP type: Date 28/05/2019
        # Working/Testing
        else:

            if len(variants.ALT) > 1 and indel_var_ann_dict[variants.POS]:

                # print variants.ALT
                # print ';'.join(set(snp_var_ann_dict[variants.POS].split(',')))

                ann_string = ";%s" % ';'.join(set(indel_var_ann_dict[variants.POS].split(',')))
                # Get Tag here; Multiple tag names.
                tag_list = []

                for i in set(indel_var_ann_dict[variants.POS].split(',')):
                    i_split = i.split('|')
                    if i_split[4] not in tag_list:
                        tag_list.append(i_split[4])

                if len(tag_list) == 1:
                    tag = tag_list[0]

                elif len(tag_list) == 2:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1])

                elif len(tag_list) == 3:
                    tag = str(tag_list[0]) + "-" + str(tag_list[1]) + "-" + str(tag_list[2])
                    print "Error: More than two locus tags were found at %s - %s" % (variants.POS, tag_list)
                    print tag_list
                elif len(tag_list) > 3:
                    print "Error: More than three locus tags were found at %s - %s" % (variants.POS, tag_list)
                    print tag_list
                    exit()

                # if len(set(snp_var_ann_dict[variants.POS].split(','))) > 2:
                #     print tag
                #     print set(snp_var_ann_dict[variants.POS].split(','))

            else:
                ann_string = ";None"

        # Changing SNP type: Date 28/05/2019
        ann_string = ann_string.replace('ERROR_OUT_OF_CHROMOSOME_RANGE', '%s-%s' % (
            locus_tag_to_gene_name[last_locus_tag], locus_tag_to_gene_name[first_locus_tag]))
        ann_string = ann_string.replace('CHR_END', '%s' % locus_tag_to_gene_name[first_locus_tag])

        # SNP Matrix Bug
        ann_string_split = ann_string.split(';')
        if len(ann_string_split) == 3:
            first_allele_ann_string_split = ann_string_split[1].split('|')
            second_allele_ann_string_split = ann_string_split[2].split('|')
            if len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) == 10:
                ann_string = ann_string
            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) == 10:
                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + first_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                #                               first_allele_ann_string_split[1] + "|" + \
                #                               first_allele_ann_string_split[2] + "|" + \
                #                               first_allele_ann_string_split[4] + "|" + \
                #                               first_allele_ann_string_split[9] + "|" + \
                #                               first_allele_ann_string_split[10] + "|" + \
                #                               first_allele_ann_string_split[11] + "|" + \
                #                               first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + \
                                              first_allele_ann_string_split[2] + "|" + \
                                              first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + \
                                              first_allele_ann_string_split[10] + "|" + \
                                              first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"

                ann_string = new_first_allele_ann_string + str(ann_string_split[2])

            elif len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) > 10:

                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + second_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                #                                second_allele_ann_string_split[1] + "|" + \
                #                                second_allele_ann_string_split[2] + "|" + \
                #                                second_allele_ann_string_split[4] + "|" + \
                #                                second_allele_ann_string_split[9] + "|" + \
                #                                second_allele_ann_string_split[10] + "|" + \
                #                                second_allele_ann_string_split[11] + "|" + \
                #                                second_allele_ann_string_split[
                #                                    13] + "|" + prod + "|" + prod + ";"
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                                               second_allele_ann_string_split[1] + "|" + \
                                               second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + \
                                               second_allele_ann_string_split[9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[
                                                   13] + "|" + prod + "|" + prod + ";"

                ann_string = str(ann_string_split[1]) + new_second_allele_ann_string
            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) > 10:

                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + first_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                #                               first_allele_ann_string_split[1] + "|" + \
                #                               first_allele_ann_string_split[2] + "|" + \
                #                               first_allele_ann_string_split[4] + "|" + \
                #                               first_allele_ann_string_split[9] + "|" + \
                #                               first_allele_ann_string_split[10] + "|" + \
                #                               first_allele_ann_string_split[11] + "|" + \
                #                               first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + \
                                              first_allele_ann_string_split[2] + "|" + \
                                              first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + \
                                              first_allele_ann_string_split[10] + "|" + \
                                              first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"

                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + second_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                #                                second_allele_ann_string_split[1] + "|" + \
                #                                second_allele_ann_string_split[2] + "|" + \
                #                                second_allele_ann_string_split[4] + "|" + \
                #                                second_allele_ann_string_split[9] + "|" + \
                #                                second_allele_ann_string_split[10] + "|" + \
                #                                second_allele_ann_string_split[11] + "|" + \
                #                                second_allele_ann_string_split[
                #                                    13] + "|" + prod + "|" + prod + ";"
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                                               second_allele_ann_string_split[1] + "|" + \
                                               second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + \
                                               second_allele_ann_string_split[9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[
                                                   13] + "|" + prod + "|" + prod + ";"

                ann_string = new_first_allele_ann_string + new_second_allele_ann_string

        # Changed >3 to ==4 Issue #29 26-03-2020
        if len(ann_string_split) == 4:

            first_allele_ann_string_split = ann_string_split[1].split('|')
            second_allele_ann_string_split = ann_string_split[2].split('|')
            third_allele_ann_string_split = ann_string_split[3].split('|')

            if len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) == 10 and len(
                    third_allele_ann_string_split) == 10:
                ann_string = ann_string

            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) == 10 and len(
                    third_allele_ann_string_split) == 10:
                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + first_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                #                               first_allele_ann_string_split[1] + "|" + \
                #                               first_allele_ann_string_split[2] + "|" + \
                #                               first_allele_ann_string_split[4] + "|" + \
                #                               first_allele_ann_string_split[9] + "|" + \
                #                               first_allele_ann_string_split[10] + "|" + \
                #                               first_allele_ann_string_split[11] + "|" + \
                #                               first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + \
                                              first_allele_ann_string_split[2] + "|" + \
                                              first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + \
                                              first_allele_ann_string_split[10] + "|" + \
                                              first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"

                ann_string = new_first_allele_ann_string + str(ann_string_split[2]) + str(ann_string_split[3])

            elif len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) > 10 and len(
                    third_allele_ann_string_split) == 10:

                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + second_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                #                                second_allele_ann_string_split[1] + "|" + \
                #                                second_allele_ann_string_split[2] + "|" + \
                #                                second_allele_ann_string_split[4] + "|" + \
                #                                second_allele_ann_string_split[9] + "|" + \
                #                                second_allele_ann_string_split[10] + "|" + \
                #                                second_allele_ann_string_split[11] + "|" + \
                #                                second_allele_ann_string_split[
                #                                    13] + "|" + prod + "|" + prod + ";"
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                                               second_allele_ann_string_split[1] + "|" + \
                                               second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + \
                                               second_allele_ann_string_split[9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[
                                                   13] + "|" + prod + "|" + prod + ";"

                ann_string = str(ann_string_split[1]) + new_second_allele_ann_string + str(ann_string_split[3])

            elif len(first_allele_ann_string_split) == 10 and len(second_allele_ann_string_split) == 10 and len(
                    third_allele_ann_string_split) > 10:

                if third_allele_ann_string_split[14] == "" and third_allele_ann_string_split[15] == "":
                    prod = third_allele_ann_string_split[3] + third_allele_ann_string_split[15]
                else:
                    prod = third_allele_ann_string_split[14] + third_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + \
                #                               third_allele_ann_string_split[1] + "|" + \
                #                               third_allele_ann_string_split[2] + "|" + \
                #                               third_allele_ann_string_split[4] + "|" + \
                #                               third_allele_ann_string_split[9] + "|" + \
                #                               third_allele_ann_string_split[10] + "|" + \
                #                               third_allele_ann_string_split[11] + "|" + \
                #                               third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + \
                                              third_allele_ann_string_split[1] + "|" + \
                                              third_allele_ann_string_split[2] + "|" + \
                                              third_allele_ann_string_split[4] + "|" + \
                                              third_allele_ann_string_split[9] + "|" + \
                                              third_allele_ann_string_split[10] + "|" + \
                                              third_allele_ann_string_split[11] + "|" + \
                                              third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"

                ann_string = str(ann_string_split[1]) + str(ann_string_split[2]) + new_third_allele_ann_string

            elif len(first_allele_ann_string_split) > 10 and len(second_allele_ann_string_split) > 10 and len(
                    third_allele_ann_string_split) > 10:
                # print ann_string
                if first_allele_ann_string_split[14] == "" and first_allele_ann_string_split[15] == "":
                    prod = first_allele_ann_string_split[3] + first_allele_ann_string_split[15]
                else:
                    prod = first_allele_ann_string_split[14] + first_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                #                               first_allele_ann_string_split[1] + "|" + \
                #                               first_allele_ann_string_split[2] + "|" + \
                #                               first_allele_ann_string_split[4] + "|" + \
                #                               first_allele_ann_string_split[9] + "|" + \
                #                               first_allele_ann_string_split[10] + "|" + \
                #                               first_allele_ann_string_split[11] + "|" + \
                #                               first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                new_first_allele_ann_string = ";" + first_allele_ann_string_split[0] + "|" + \
                                              first_allele_ann_string_split[1] + "|" + \
                                              first_allele_ann_string_split[2] + "|" + \
                                              first_allele_ann_string_split[4] + "|" + \
                                              first_allele_ann_string_split[9] + "|" + \
                                              first_allele_ann_string_split[10] + "|" + \
                                              first_allele_ann_string_split[11] + "|" + \
                                              first_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"

                if second_allele_ann_string_split[14] == "" and second_allele_ann_string_split[15] == "":
                    prod = second_allele_ann_string_split[3] + second_allele_ann_string_split[15]
                else:
                    prod = second_allele_ann_string_split[14] + second_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                #                                second_allele_ann_string_split[1] + "|" + \
                #                                second_allele_ann_string_split[2] + "|" + \
                #                                second_allele_ann_string_split[4] + "|" + \
                #                                second_allele_ann_string_split[9] + "|" + \
                #                                second_allele_ann_string_split[10] + "|" + \
                #                                second_allele_ann_string_split[11] + "|" + \
                #                                second_allele_ann_string_split[
                #                                    13] + "|" + prod + "|" + prod + ";"
                new_second_allele_ann_string = second_allele_ann_string_split[0] + "|" + \
                                               second_allele_ann_string_split[1] + "|" + \
                                               second_allele_ann_string_split[2] + "|" + \
                                               second_allele_ann_string_split[4] + "|" + \
                                               second_allele_ann_string_split[9] + "|" + \
                                               second_allele_ann_string_split[10] + "|" + \
                                               second_allele_ann_string_split[11] + "|" + \
                                               second_allele_ann_string_split[
                                                   13] + "|" + prod + "|" + prod + ";"

                if third_allele_ann_string_split[14] == "" and third_allele_ann_string_split[15] == "":
                    prod = third_allele_ann_string_split[3] + third_allele_ann_string_split[15]
                else:
                    prod = third_allele_ann_string_split[14] + third_allele_ann_string_split[15]
                # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                # new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + \
                #                               third_allele_ann_string_split[1] + "|" + \
                #                               third_allele_ann_string_split[2] + "|" + \
                #                               third_allele_ann_string_split[4] + "|" + \
                #                               third_allele_ann_string_split[9] + "|" + \
                #                               third_allele_ann_string_split[10] + "|" + \
                #                               third_allele_ann_string_split[11] + "|" + \
                #                               third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"
                new_third_allele_ann_string = third_allele_ann_string_split[0] + "|" + \
                                              third_allele_ann_string_split[1] + "|" + \
                                              third_allele_ann_string_split[2] + "|" + \
                                              third_allele_ann_string_split[4] + "|" + \
                                              third_allele_ann_string_split[9] + "|" + \
                                              third_allele_ann_string_split[10] + "|" + \
                                              third_allele_ann_string_split[11] + "|" + \
                                              third_allele_ann_string_split[13] + "|" + prod + "|" + prod + ";"

                ann_string = new_first_allele_ann_string + new_second_allele_ann_string + new_third_allele_ann_string

        # Added this extra check Issue #29 26-03-2020
        if len(ann_string_split) > 4:
            new_allele_string_array = []
            for i_ann_string_split in ann_string_split[1:]:
                if len(i_ann_string_split.split('|')) == 10:
                    ann_string = ann_string
                elif len(i_ann_string_split.split('|')) > 10:
                    ann_string = ";"
                    i_ann_string_split_array = i_ann_string_split.split('|')
                    if i_ann_string_split_array[14] == "" and i_ann_string_split_array[15] == "":
                        prod = i_ann_string_split_array[3] + i_ann_string_split_array[15]
                    else:
                        prod = i_ann_string_split_array[14] + i_ann_string_split_array[15]

                    # Remove gene symbol and insert Locus Tag in field 4 #30 - 2020-03-27
                    # new_allele_string = i_ann_string_split_array[0] + "|" + i_ann_string_split_array[1] + "|" + \
                    #                     i_ann_string_split_array[2] + "|" + \
                    #                     i_ann_string_split_array[4] + "|" + \
                    #                     i_ann_string_split_array[9] + "|" + \
                    #                     i_ann_string_split_array[10] + "|" + \
                    #                     i_ann_string_split_array[11] + "|" + \
                    #                     i_ann_string_split_array[13] + "|" + prod + "|" + prod
                    new_allele_string = i_ann_string_split_array[0] + "|" + i_ann_string_split_array[1] + "|" + \
                                        i_ann_string_split_array[2] + "|" + \
                                        i_ann_string_split_array[4] + "|" + \
                                        i_ann_string_split_array[9] + "|" + \
                                        i_ann_string_split_array[10] + "|" + \
                                        i_ann_string_split_array[11] + "|" + \
                                        i_ann_string_split_array[13] + "|" + prod + "|" + prod
                    new_allele_string_array.append(new_allele_string)

                    ann_string = ann_string + ";".join(new_allele_string_array)
        # # JUST FOR THE SAKE OF DEBUGGING
        # ann_string_split = ann_string.split(';')
        # for i in ann_string_split:
        #     if len(i.split('|')) != 10 and len(i.split('|')) != 1:
        #         print ann_string

        # Changing Strandness string: Date 28/05/2019
        # Each Locus ID with a strand information
        strandness = " Strand Information: "
        if "-" in tag:
            tagsplit = tag.split('-')
            for i in tagsplit:
                if i in locus_tag_to_strand.keys():
                    if "," in locus_tag_to_strand[i]:
                        locus_tag_to_strand_split = locus_tag_to_strand[i].split(',')
                        strand = locus_tag_to_strand_split[0]
                    else:
                        strand = locus_tag_to_strand[i]
                    strandness = strandness + i + "=" + strand + "/"
                else:
                    if i == "" or i == "None":
                        strandness = strandness + "NULL=" + "No Strand Information found" + "/"
                    else:
                        strandness = strandness + i + "=" + "No Strand Information found" + "/"
        else:
            if tag in locus_tag_to_strand.keys():
                # strandness = strandness + locus_tag_to_strand[tag]
                if "," in locus_tag_to_strand[tag]:
                    locus_tag_to_strand_split = locus_tag_to_strand[tag].split(',')
                    strand = locus_tag_to_strand_split[0]
                else:
                    strand = locus_tag_to_strand[tag]
                strandness = strandness + tag + "=" + strand
            else:
                if tag == "" or tag == "None":
                    strandness = strandness + "NULL=" + "No Strand Information found"
                else:
                    strandness = strandness + tag + "=" + "No Strand Information found"

        # Changing tag equals NULL: Date 30/05/2019
        if tag == "" or tag == "None":
            tag = "NULL"

        # # Test No. of Fields
        ann_string_split = ann_string.split(';')
        # ann_string = ";"
        for i in ann_string_split:
            ann_string_array_split_columns = i.split('|')
            if len(ann_string_array_split_columns) > 10:
                print "Warning: More than 10 field - %s %s" % (variants.POS, i)
                # print variants.POS
                # ann_string = ann_string.replace('||WARNING_TRANSCRIPT_NO_START_CODON||WARNING_TRANSCRIPT_NO_START_CODON', '|WARNING_TRANSCRIPT_NO_START_CODON|WARNING_TRANSCRIPT_NO_START_CODON')
                # ann_string = ann_string + str('|'.join(ann_string_array_split_columns[: -2 or None])) + ";"
                # print ann_string
        print_string = print_string + " locus_tag=" + tag + strandness + ann_string
        print_string_phage = print_string
        
        gt_string = ""
        for gt in variants.gt_bases:
            gt = gt.replace('./.', '.')
            if "/" in gt:
                gt_split = gt.split('/')
                gt = gt_split[1]
            gt_string = gt_string + "," + gt
        gt_string = gt_string.replace('.', variants.REF)

        """Replacing Phage/Functional filter position code"""
        if str(variants.POS) in functional_filter_pos_array:
            code_string_array = code_string.split(',')
            code_string = ""
            for i in code_string_array:
                code_string = code_string + "," + "-2"

        final_allele_string = print_string + gt_string.replace(',', '\t') + '\n'
        final_code_string = print_string + "\t" + code_string.replace(',', '\t') + '\n'
        final_code_string_unmasked = print_string + "\t" + code_string_unmasked.replace(',', '\t') + '\n'
        final_allele_string = final_allele_string.replace(',|', '|')
        # final_allele_string = final_allele_string.replace(',;,', ':::')
        # final_allele_string = final_allele_string.replace(';,', ':::')
        final_allele_string = final_allele_string.replace(',;,', ':::')
        final_allele_string = final_allele_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(',|', '|')
        final_code_string_unmasked = final_code_string_unmasked.replace(',|', '|')
        # final_code_string = final_code_string.replace(',;,', ':::')
        # final_code_string = final_code_string.replace(';,', ':::')
        final_code_string = final_code_string.replace(',;,', ':::')
        final_code_string = final_code_string.replace(';,', ':::')
        final_code_string = final_code_string.replace('\t\t', '\t')
        final_code_string_unmasked = final_code_string_unmasked.replace(',;,', ':::')
        final_code_string_unmasked = final_code_string_unmasked.replace(';,', ':::')
        final_code_string_unmasked = final_code_string_unmasked.replace('\t\t', '\t')

        final_allele_string = final_allele_string.replace('\t\t', '\t')
        fp_allele.write(final_allele_string)
        fp_code.write(final_code_string)
        fp_code_unmasked.write(final_code_string_unmasked)

    fp_code.close()
    fp_code_unmasked.close()
    fp_allele.close()
    method_time_taken = datetime.now() - method_start_time

    keep_logging('Time taken to complete the Generate Indel matrix method: {}'.format(method_time_taken),
                 'Time taken to complete the Generate Indel matrix method: {}'.format(method_time_taken), logger, 'info')

def annotated_snp_matrix():
    """
    :return: Annotate core vcf files generated at core_prep steps.
    Read Genbank file and return a dictionary of Prokka ID mapped to Gene Name, Prokka ID mapped to Product Name.
    This dictionary will then be used to insert annotation into SNP/Indel matrix
    """
    method_start_time = datetime.now()
    """Annotate all VCF file formats with SNPeff"""
    # Removing snpEff steps from this script. 
    # snpEff will be individually run as a part of variant calling steps.
    # variant_annotation()
    # indel_annotation()

    keep_logging('- Adding snpEff annotations to SNP/Indel Matrix.', '- Adding snpEff annotations to SNP/Indel Matrix.',
                 logger, 'info')

    locus_tag_to_gene_name, locus_tag_to_product, locus_tag_to_strand, first_locus_tag, last_element, last_locus_tag = extract_locus_tag_from_genbank()

    # Commented for debugging
    annotated_no_proximate_snp_file, annotated_no_proximate_snp_indel_file, final_gatk_snp_merged_vcf, final_gatk_indel_merged_vcf = merge_vcf()

    annotated_no_proximate_snp_file = "%s/annotated_no_proximate_snp_list.txt" % args.filter2_only_snp_vcf_dir
    annotated_no_proximate_snp_indel_file = "%s/annotated_no_proximate_snp_indel_list.txt" % args.filter2_only_snp_vcf_dir
    final_gatk_snp_merged_vcf = "Final_vcf_gatk_no_proximate_snp.vcf"
    final_gatk_indel_merged_vcf = "Final_vcf_gatk_indel.vcf"
    
    snp_var_ann_dict, indel_var_ann_dict = extract_annotations_from_multivcf()

    core_positions, indel_core_positions = extract_core_positions()

    functional_filter_pos_array, phage_positions, repetitive_positions, mask_positions = extract_functional_class_positions()

    """ Read and parse final GATK merged vcf file cyvcf library; Generate a header string from the sample list of this merged vcf file"""

    final_merge_anno_file = VCF("%s/Final_vcf_gatk_no_proximate_snp.vcf.gz" % args.filter2_only_snp_vcf_dir)

    """ End """

    position_label, position_indel_label = generate_position_label_dict(final_merge_anno_file)

    # Commented for debugging
    #mask_fq_mq_positions, mask_fq_mq_positions_outgroup_specific = get_low_fq_mq_positions(position_label)

    generate_SNP_matrix(final_merge_anno_file, functional_filter_pos_array, phage_positions, repetitive_positions,
                        mask_positions, position_label, core_positions, snp_var_ann_dict)

    mask_fq_mq_positions, mask_fq_mq_positions_outgroup_specific = get_low_fq_mq_positions_indel(position_label,
                                                                                                 position_indel_label)

    final_merge_anno_file = VCF("%s/Final_vcf_gatk_indel.vcf.gz" % args.filter2_only_snp_vcf_dir)

    # print "Generating Indel Matrix"

    generate_Indel_matrix(final_merge_anno_file, functional_filter_pos_array, phage_positions, repetitive_positions,
                          mask_positions, position_indel_label, indel_core_positions, indel_var_ann_dict)
    method_time_taken = datetime.now() - method_start_time

    keep_logging('Time taken to complete the Annotated SNP/Indel Matrix method: {}'.format(method_time_taken),
                 'Time taken to complete the Annotated SNP/Indel Matrix method: {}'.format(method_time_taken), logger, 'info')

""" Report methods """


def alignment_report(data_matrix_dir):
    keep_logging('Generating Alignment report...', 'Generating Alignment report...', logger, 'info')
    varcall_dir = os.path.dirname(args.results_dir)
    # print varcall_dir
    report_string = ""
    header = "Sample,QC-passed reads,Mapped reads,% mapped reads,mean depth,%_bases_above_5,%_bases_above_10,%_bases_above_15,unmapped_positions,READ_PAIR_DUPLICATES,READ_PAIR_OPTICAL_DUPLICATES,unmapped reads,% unmapped reads"
    fp = open("%s/Report_alignment.txt" % (data_matrix_dir), 'w+')
    fp.write(header + '\n')
    for vcf in vcf_filenames:
        sample = os.path.basename(vcf.replace('_filter2_final.vcf_no_proximate_snp.vcf', ''))
        # print sample
        report_string = sample + ","
        qc = (subprocess.check_output(
            "grep \'QC-passed\' %s/%s/%s_alignment_stats | sed \'s/ + 0 in total (QC-passed reads + QC-failed reads)//g\'" % (
            varcall_dir, sample, sample), shell=True)).strip()
        mapped = (subprocess.check_output(
            "grep \'mapped (\' %s/%s/%s_alignment_stats | awk -F\' \' \'{print $1}\'" % (varcall_dir, sample, sample),
            shell=True)).strip()
        replace = "%:-nan%)"
        perc_mapped = (subprocess.check_output(
            "grep \'mapped (\' %s/%s/%s_alignment_stats | awk -F\' \' \'{print $5}\' | sed \'s/%s//g\' | sed \'s/(//g\'" % (
            varcall_dir, sample, sample, replace), shell=True)).strip()
        depth_of_coverage = (subprocess.check_output(
            "awk -F\'\\t\' \'{OFS=\",\"};FNR==2{print $3,$7,$8,$9}\' %s/%s/%s_depth_of_coverage.sample_summary" % (
            varcall_dir, sample, sample), shell=True)).strip()
        unmapped_positions = (subprocess.check_output(
            "wc -l %s/%s/%s_unmapped.bed_positions | cut -d\' \' -f1" % (varcall_dir, sample, sample),
            shell=True)).strip()
        opt_dup = (subprocess.check_output(
            "awk -F\'\\t\' \'{OFS=\",\"};FNR==8{print $7,$8,$5}\' %s/%s/%s_markduplicates_metrics" % (
            varcall_dir, sample, sample), shell=True)).strip()
        perc_unmapped = str(100 - float(perc_mapped))
        myList = ','.join(
            map(str, (sample, qc, mapped, perc_mapped, depth_of_coverage, unmapped_positions, opt_dup, perc_unmapped)))
        # print myList
        fp.write(myList + '\n')
    fp.close()
    keep_logging('Alignment report can be found in %s/Report_alignment.txt' % data_matrix_dir,
                 'Alignment report can be found in %s/Report_alignment.txt' % data_matrix_dir, logger, 'info')


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
        unmapped_positions = (
            subprocess.check_output("wc -l %s/core_temp_dir/unique_positions_file | cut -d\' \' -f1" % (varcall_dir),
                                    shell=True)).strip()
        core_snps = (subprocess.check_output(
            "wc -l %s/core_temp_dir/Only_ref_variant_positions_for_closely | cut -d\' \' -f1" % (varcall_dir),
            shell=True)).strip()
        filtered_snp_count = (subprocess.check_output(
            "grep -w \'^%s\' %s/core_temp_dir/bargraph_counts.txt | awk -F\'\\t\' \'{OFS=\",\"};{print $2,$3,$4,$5,$6,$7}\'" % (
            sample, varcall_dir), shell=True)).strip()
        filtered_snp_perc = (subprocess.check_output(
            "grep -w \'^%s\' %s/core_temp_dir/bargraph_percentage.txt | awk -F\'\\t\' \'{OFS=\",\"};{print $2,$3,$4,$5,$6,$7}\'" % (
            sample, varcall_dir), shell=True)).strip()
        myList = ','.join(map(str, (sample, unmapped_positions, core_snps, filtered_snp_count, filtered_snp_perc)))
        fp.write(myList + '\n')
    fp.close()
    keep_logging('Variant call report can be found in %s/Report_variants.txt' % data_matrix_dir,
                 'Variant call report can be found in %s/Report_variants.txt' % data_matrix_dir, logger, 'info')


""" Gubbins/Iqtree methods"""


def gubbins(gubbins_dir, input_fasta, jobrun, logger, Config):
    keep_logging('\nRunning Gubbins on input: %s\n' % input_fasta, '\nRunning Gubbins on input: %s\n' % input_fasta,
                 logger,
                 'info')

    call("module load bioperl python-anaconda2/201607 biopython dendropy reportlab fasttree RAxML fastml/gub gubbins",
         logger)
    # os.system("module load bioperl python-anaconda2/201607 biopython dendropy reportlab fasttree RAxML fastml/gub gubbins")
    # gubbins_cmd = "%s/%s --prefix %s/%s %s" % (
    # ConfigSectionMap("gubbins", Config)['gubbins_bin'], ConfigSectionMap("gubbins", Config)['base_cmd'], gubbins_dir,
    # (os.path.basename(input_fasta)).replace('.fa', ''), input_fasta)

    load_module = "module load bioperl python-anaconda2/201607 biopython dendropy reportlab fasttree RAxML fastml/gub gubbins"
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
        job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=12,mem=47000mb,walltime=250:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\ncd %s\n%s\n%s" % (
        job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'],
        ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
        gubbins_dir, load_module, gubbins_cmd)
        f1 = open(job_file_name, 'w+')
        f1.write(job_print_string)
        f1.close()
        # os.system("qsub %s" % job_file_name)
        call("qsub %s" % job_file_name, logger)


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

    # Updated - Initialize Global variables
    global logger

    analysis_name_log = "step_" + str(args.steps)
    logger = generate_logger(args.filter2_only_snp_vcf_dir, analysis_name_log, log_unique_time)
    print_details = "This step will parse final vcf files(*_no_proximate_snp.vcf) generated at the end of Variant Calling Pipeline. At the end of this step, the following results will be generated and placed in output directory:\n\n" \
                    "1. Final Core SNP Positions list(Variant positions that were not filtered out in any of the samples and passed all the filters)\n" \
                    "2. SNP Positions that were filtered out with labels indicating the reason (Depth, FQ, MQ, Unmapped in one or other samples, Proximate SNPS, Quality of Variant) why they were filtered out.\n" \
                    "3. Barplot Statistics about the filtered variants and their reason for getting filtered.\n" \
                    "4. Final Consensus fasta file using only Core SNP Positions\n"
    #keep_logging('%s' % print_details, '%s' % print_details, logger, 'info')
    keep_logging('- Parsing VCF files generated with SNPKIT', '- Parsing VCF files generated with SNPKIT', logger, 'info')

    # Create temporary Directory /tmp/snpkit_temp for storing temporary intermediate files. Check if core_temp_dir contains all the required files to run these pipeline.
    global temp_dir
    temp_dir = "/tmp/snpkit_temp"
    make_sure_path_exists(temp_dir)

    # Read Config file into Config object that will be used to extract configuration settings set up in config file.
    global config_file
    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
    global Config
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    keep_logging('- Path to config file: %s' % config_file, '- Path to config file: %s' % config_file, logger, 'info')

    global num_cores
    if args.numcores:
        num_cores = int(args.numcores)
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

    global bin_dir
    proc = subprocess.Popen(["which gatk"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    bin_dir = os.path.dirname(out)

    global outgroup_specific_positions
    global outgroup_indel_specific_positions

    # Get outgroup_Sample name
    outgroup = get_outgroup()
    outgroup_vcf_filename = str(outgroup) + "_filter2_final.vcf_no_proximate_snp.vcf"
    outgroup_indel_vcf_filename = str(outgroup) + "_filter2_indel_final.vcf"

    # Read filenames. Core variants and final results will be extracted considering only these files.
    filter2_only_snp_vcf_filenames = args.filter2_only_snp_vcf_filenames
    vcf_filenames_temp = []
    vcf_filenames_temp_outgroup = []

    with open(filter2_only_snp_vcf_filenames) as fp:
        for line in fp:
            line = line.strip()
            #line = args.filter2_only_snp_vcf_dir + line
            line = (args.filter2_only_snp_vcf_dir).replace('core_temp_dir', '*/*_vcf_results') + line
            line = glob.glob(line)
            vcf_filenames_temp.append(line[0])
            if args.outgroup:
                if "%s_filter2_final.vcf_no_proximate_snp.vcf" % outgroup not in line:
                    vcf_filenames_temp_outgroup.append(line)
        fp.close()
    vcf_filenames = sorted(vcf_filenames_temp)
    vcf_filenames_outgroup = sorted(vcf_filenames_temp_outgroup)

    make_sure_files_exists(vcf_filenames, Config, logger)

    log_file_handle = "%s/%s_%s.log.txt" % (args.filter2_only_snp_vcf_dir, log_unique_time, analysis_name_log)

    scheduler_directives, script_Directive, job_name_flag = get_scheduler_directive(args.scheduler, Config)

    global functional_filter_pos_array
    functional_filter_pos_array = []
    functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % args.filter2_only_snp_vcf_dir

    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
        keep_logging(
            '- Extracting Functional Class filter regions such as Phage (Phaster), Repeat (MUMmer), Custom Mask to %s' % os.path.basename(functional_class_filter_positions),
            '- Extracting Functional Class filter regions such as Phage (Phaster), Repeat (MUMmer), Custom Mask to %s' % os.path.basename(functional_class_filter_positions),
            logger,
            'info')

        
        with open(functional_class_filter_positions, 'rU') as f_functional:
            for line_func in f_functional:
                functional_filter_pos_array.append(line_func.strip())

        keep_logging(
            '- Number of Functional Class filter positions: %s' % len(functional_filter_pos_array),
            '- Number of Functional Class filter positions: %s' % len(functional_filter_pos_array),
            logger,
            'info')
    # Start Variant Calling Core Pipeline steps based on steps argument supplied.
    if "1" in args.steps:
        """ 
        core_prep step
        """
        method_start_time = datetime.now()
        # Gather SNP positions from each final *_no_proximate_snp.vcf file (that passed the variant filter parameters from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files for use in downstream methods
        keep_logging('- Gathering Filtered SNP positions information.',
                     '- Gathering Filtered SNP positions information.', logger,
                     'info')

        # Extract All the unique SNO and Indel position list from final filtered *_no_proximate_snp.vcf files.
        unique_position_file = create_positions_filestep(vcf_filenames, temp_dir, args.outgroup, logger,
                                                         args.filter2_only_snp_vcf_dir, num_cores)
        unique_indel_position_file = create_indel_positions_filestep(vcf_filenames, temp_dir, args.outgroup, logger,
                                                                     args.filter2_only_snp_vcf_dir, num_cores)

        # bgzip and tabix all the vcf files in core_temp_dir.
        files_for_tabix = glob.glob((args.filter2_only_snp_vcf_dir).replace('core_temp_dir/', '*/*_vcf_results/*.vcf'))
        tabix(files_for_tabix, "vcf", logger, Config)

        # Get the cluster option; create and run jobs based on given parameter. The jobs will parse all the intermediate vcf file to extract information such as if any unique variant position was unmapped in a sample, if it was filtered out dur to DP,MQ, FQ, proximity to indel, proximity to other SNPs and other variant filter parameters set in config file.
        tmp_dir = "/tmp/temp_%s/" % log_unique_time
        make_sure_path_exists(tmp_dir)

        create_job(args.jobrun, vcf_filenames, unique_position_file, temp_dir, scheduler_directives, script_Directive,
                   job_name_flag, temp_dir, args.outgroup, logger, args.filter2_only_snp_vcf_dir, num_cores)

        create_indel_job(args.jobrun, vcf_filenames, unique_indel_position_file, temp_dir, scheduler_directives,
                         script_Directive, job_name_flag, temp_dir, args.outgroup, logger,
                         args.filter2_only_snp_vcf_dir, num_cores)

        # If Phaster Summary file doesn't exist in reference genome folder
        if not os.path.isfile("%s/summary.txt" % os.path.dirname(args.reference)):
            if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
                keep_logging('Functional class filter is set to yes. Preparing Functional class filters\n',
                             'Functional class filter is set to yes. Preparing Functional class filters\n', logger,
                             'info')
                if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
                    # Submit Phaster jobs to find ProPhage region in reference genome.
                    run_phaster(args.reference, args.filter2_only_snp_vcf_dir, logger, Config)

        call("cp %s %s/Logs/core_prep/" % (
        log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)
        method_time_taken = datetime.now() - method_start_time

        keep_logging('- Time taken to complete the Core Step 1 method: {}'.format(method_time_taken),
                     '- Time taken to complete the Core Step 1 method: {}'.format(method_time_taken), logger, 'info')

    if "2" in args.steps:
        """ 
        core step 
        """
        method_start_time = datetime.now()
        # Set variables; check if the output from core_prep steps (*label files) exists and was completed without any errors.
        snp_unique_positions_file = args.filter2_only_snp_vcf_dir + "/unique_positions_file"
        indel_unique_positions_file = args.filter2_only_snp_vcf_dir + "/unique_indel_positions_file"
        uniq_snp_positions = sum(1 for line in open('%s' % snp_unique_positions_file))
        uniq_indel_positions = sum(1 for line in open('%s' % indel_unique_positions_file))
        if not os.path.isfile(snp_unique_positions_file) and not os.path.isfile(indel_unique_positions_file):
            keep_logging(
                '- Error finding unique_positions_file/unique_indel_positions_file. Please rerun core_prep step.',
                '- Error finding unique_positions_file/unique_indel_positions_file. Please rerun core_prep step.', logger,
                'exception')
            exit()
        make_sure_label_files_exists(vcf_filenames, uniq_snp_positions, uniq_indel_positions, Config, logger)

        # Set up Report and results directories to transfer the final results.
        data_matrix_dir = args.results_dir + '/data_matrix'
        core_vcf_fasta_dir = args.results_dir + '/core_snp_consensus'
        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(core_vcf_fasta_dir)

        functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % args.filter2_only_snp_vcf_dir

        # Get outgroup specific variant positions
        if args.outgroup:
            f_outgroup = open("%s/outgroup_indel_specific_positions.txt" % args.filter2_only_snp_vcf_dir, 'r+')

            outgroup_indel_specific_positions = []
            for i in f_outgroup:
                i = i.strip()
                outgroup_indel_specific_positions.append(int(i))
            f_outgroup.close()

            f_outgroup = open("%s/outgroup_specific_positions.txt" % args.filter2_only_snp_vcf_dir, 'r+')

            outgroup_specific_positions = []
            for i in f_outgroup:
                i = i.strip()
                outgroup_specific_positions.append(int(i))
            f_outgroup.close()

            print "No. of outgroup specific variant positions: %s" % len(outgroup_specific_positions)
            print "No. of outgroup specific Indel variant positions: %s" % len(outgroup_indel_specific_positions)
        else:
            outgroup_indel_specific_positions = []
            outgroup_specific_positions = []
            #print "No. of outgroup specific variant positions: %s" % len(outgroup_specific_positions)
            #print "No. of outgroup specific Indel variant positions: %s" % len(outgroup_indel_specific_positions)

        # Commented out for debugging
        # Run core steps. Generate SNP and data Matrix results. Extract core SNPS and consensus files.
        Only_ref_indel_positions_for_closely = core_prep_indel()

        Only_ref_variant_positions_for_closely = core_prep_snp()

        # Annotate core variants. Generate SNP and Indel matrix.
        # Commented out for debugging
        annotated_snp_matrix()

        # # Read new allele matrix and generate fasta; generate a seperate function
        # keep_logging('Generating Fasta from Variant Alleles...\n', 'Generating Fasta from Variant Alleles...\n', logger,
        #              'info')

        # create_job_allele_variant_fasta(args.jobrun, vcf_filenames, args.filter2_only_snp_vcf_dir, config_file,
        #                                 script_Directive, job_name_flag)

        # extract_only_ref_variant_fasta_from_reference_allele_variant()

        # # mask_fq_mq_positions_specific_to_outgroup()

        # call("cp %s %s/Logs/core/" % (
        #     log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)
        # method_time_taken = datetime.now() - method_start_time

        # keep_logging('Time taken to complete the Core Step 2 method: {}'.format(method_time_taken),
        #              'Time taken to complete the Core Step 2 method: {}'.format(method_time_taken), logger, 'info')

    if "3" in args.steps:
        """ report step """

        # Get outgroup_Sample name
        outgroup = get_outgroup()

        keep_logging('Step 3: Generate Results folder.', 'Step 3: Generate Results folder.',
                     logger, 'info')

        # Set up Report and results directories to transfer the final results.
        keep_logging('Step 3: Setting up Result directory - %s' % args.results_dir, 'Step 3: Setting up Result directory - %s' % args.results_dir,
                     logger, 'info')
        data_matrix_dir = args.results_dir + '/data_matrix'
        plots_dir = "%s/plots" % data_matrix_dir
        matrices_dir = "%s/matrices" % data_matrix_dir
        functional_ann_dir = "%s/Functional_annotation_results" % data_matrix_dir
        logs_dir = "%s/logs" % data_matrix_dir
        data_matrix_snpeff_dir = data_matrix_dir + '/snpEff_results'

        # Create the above directories 
        make_sure_path_exists(args.results_dir)
        make_sure_path_exists(data_matrix_dir)
        make_sure_path_exists(plots_dir)
        make_sure_path_exists(matrices_dir)
        make_sure_path_exists(functional_ann_dir)
        make_sure_path_exists(logs_dir)
        make_sure_path_exists(data_matrix_snpeff_dir)

        # Move results to the results directory
        # List of files to move
        
        keep_logging('Step 3: Copying data to - %s' % args.results_dir, 'Step 3: Copying data to - %s' % args.results_dir,
                     logger, 'info')

        list_of_data_matrix_output_files = ['unique_positions_file', 'unique_indel_positions_file', 
        'All_label_final_ordered_sorted.txt', 'All_indel_label_final_ordered_sorted.txt',
        'SNP_matrix_allele_new.tsv', 'SNP_matrix_allele_unmasked.tsv', 'SNP_matrix_code.tsv', 'SNP_matrix_code_unmasked.tsv',
        'Indel_matrix_allele.tsv', 'Indel_matrix_code.tsv', 'Indel_matrix_code_unmasked.tsv',
        'mask_fq_mq_positions_outgroup_specific.txt', 'mask_fq_mq_positions.txt']

        list_of_snpeff_results_output_files = ['*_aln_mpileup_raw.vcf_ANN.vcf', '*_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf', '*_filter2_indel_final.vcf_ANN.vcf']

        list_of_functional_ann_output_files = ['Functional_class_filter_positions.txt', 'inexact_repeat_region_positions.txt', 'phage_region_positions.txt', 'repeat_region_positions.txt']
                
        for files in list_of_data_matrix_output_files:
            os.system("cp %s/%s %s/matrices" % (args.filter2_only_snp_vcf_dir, files, data_matrix_dir))
        
        for files in list_of_snpeff_results_output_files:
            os.system("cp %s/%s %s" % (args.filter2_only_snp_vcf_dir, files, data_matrix_snpeff_dir))
        
        for files in list_of_functional_ann_output_files:
            os.system("cp %s/%s %s" % (args.filter2_only_snp_vcf_dir, files, functional_ann_dir))

        """ Generating Gubbins MFA files"""
        keep_logging('Step 3: Generating Gubbins MFA files', 'Generating Gubbins MFA files',
                     logger, 'info')
        reference_base = os.path.basename(args.reference).split('.')[0]
        gubbins_dir = args.results_dir + '/gubbins'

        make_sure_path_exists(gubbins_dir)
        
        prepare_ref_allele_unmapped_consensus_input = "%s/gubbins/%s_%s_genome_aln_w_alt_allele_unmapped.fa" % (
            args.results_dir, (os.path.basename(os.path.normpath(args.results_dir))).replace('_core_results', ''),
            reference_base)

        prepare_ref_allele_unmapped_consensus_input_cmd = "cat %s %s/*_ref_allele_unmapped_variants.fa > %s" % (
            args.reference, args.filter2_only_snp_vcf_dir, prepare_ref_allele_unmapped_consensus_input)
        
        call("%s" % prepare_ref_allele_unmapped_consensus_input_cmd, logger)

        os.chdir(gubbins_dir)

        if args.scheduler == "SLURM":
            job_file_name = "%s" % (prepare_ref_allele_unmapped_consensus_input.replace('.fa', '.sbat'))
        elif args.scheduler == "PBS":
            job_file_name = "%s" % (prepare_ref_allele_unmapped_consensus_input.replace('.fa', '.pbs'))
        else:
            job_file_name = "%s" % (prepare_ref_allele_unmapped_consensus_input.replace('.fa', '.sbat'))

        if args.mask:
            gubbins_iqtree_script = "python %s/scripts/gubbins_iqtree.py -w %s" % (
            os.path.dirname(os.path.abspath(__file__)), prepare_ref_allele_unmapped_consensus_input)
            if args.outgroup:
                # Get outgroup_Sample name
                outgroup = get_outgroup()
                gubbins_iqtree_script = gubbins_iqtree_script + " -o %s" % outgroup
            print gubbins_iqtree_script
            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(job_file_name))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s" % os.path.dirname(prepare_ref_allele_unmapped_consensus_input) + '\n')
                out.write(gubbins_iqtree_script + '\n')
            out.close()

        else:
            iqtree_results_dir = args.results_dir + '/gubbins/iqtree_results'
            gubbins_command = "run_gubbins.py --prefix %s --threads %s %s" % (
            os.path.basename(prepare_ref_allele_unmapped_consensus_input).replace('.fa', ''), num_cores,
            prepare_ref_allele_unmapped_consensus_input)
            if args.outgroup:
                # Get outgroup_Sample name
                outgroup = get_outgroup()
                gubbins_command = gubbins_command + " --outgroup %s" % outgroup
            iqtree_command = "iqtree -s %s/%s.filtered_polymorphic_sites.fasta -nt AUTO -bb 1000 -m MFP -pre %s/%s" % (
            os.path.dirname(prepare_ref_allele_unmapped_consensus_input),
            os.path.basename(prepare_ref_allele_unmapped_consensus_input).replace('.fa', ''), iqtree_results_dir,
            os.path.basename(prepare_ref_allele_unmapped_consensus_input.replace('.fa', '')))

            with open(job_file_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(job_file_name))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s" % os.path.dirname(prepare_ref_allele_unmapped_consensus_input) + '\n')
                out.write(gubbins_command + '\n')
                out.write(iqtree_command + '\n')
            out.close()

        keep_logging('sbatch %s' % job_file_name, 'sbatch %s' % job_file_name, logger, 'info')

        call("cp %s %s/Logs/tree/" % (
            log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)


        call("cp %s/%s*.log.txt %s" % (args.filter2_only_snp_vcf_dir, datetime.now().strftime('%Y_%m_%d'), logs_dir), logger)

        

        call("cp %s %s/Logs/report/" % (
            log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)

    """ The below steps are for debugging purpose only."""
    if "5" in args.steps:
        """ 
        Debugging Purposes only: Run only SNP matrix annotation step 
        """

        keep_logging('Step 5: Running SNP matrix annotation step.', 'Step 5: Running SNP matrix annotation step.',
                     logger, 'info')

        functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % args.filter2_only_snp_vcf_dir

        
        functional_filter_pos_array = pd.read_csv(functional_class_filter_positions, sep='\n', header=None)
        functional_filter_pos_array = functional_filter_pos_array.squeeze()
        Only_ref_variant_positions_for_closely = pd.read_csv("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, sep='\n', header=None)
        Only_ref_variant_positions_for_closely = Only_ref_variant_positions_for_closely.squeeze()
        exclude_ref_var_functional = pd.Series(np.intersect1d(Only_ref_variant_positions_for_closely.values,functional_filter_pos_array.values))
        Only_ref_variant_positions_for_closely_without_functional_filtered_positions = Only_ref_variant_positions_for_closely[~Only_ref_variant_positions_for_closely.isin(exclude_ref_var_functional)]
        Only_ref_variant_positions_for_closely_without_functional_filtered_positions.to_csv('%s/Only_ref_variant_positions_for_closely_without_functional_filtered_positions' % args.filter2_only_snp_vcf_dir, index=False, sep='\n', header=None)
        
        # Get outgroup specific variant positions
        if args.outgroup:
            f_outgroup = open("%s/outgroup_indel_specific_positions.txt" % args.filter2_only_snp_vcf_dir, 'r+')

            outgroup_indel_specific_positions = []
            for i in f_outgroup:
                i = i.strip()
                outgroup_indel_specific_positions.append(int(i))
            f_outgroup.close()

            f_outgroup = open("%s/outgroup_specific_positions.txt" % args.filter2_only_snp_vcf_dir, 'r+')

            outgroup_specific_positions = []
            for i in f_outgroup:
                i = i.strip()
                outgroup_specific_positions.append(int(i))
            f_outgroup.close()

            print "No. of outgroup specific variant positions: %s" % len(outgroup_specific_positions)
            print "No. of outgroup specific Indel variant positions: %s" % len(outgroup_indel_specific_positions)
        else:

            outgroup_indel_specific_positions = []
            outgroup_specific_positions = []

        # Annotate core variants. Generate SNP and Indel matrix.
        annotated_snp_matrix()

        # # Read new allele matrix and generate fasta; generate a seperate function
        keep_logging('Generating Fasta from Variant Alleles...\n', 'Generating Fasta from Variant Alleles...\n', logger,
                     'info')

        create_job_allele_variant_fasta(args.jobrun, vcf_filenames, args.filter2_only_snp_vcf_dir, config_file,
                                        script_Directive, job_name_flag)

        extract_only_ref_variant_fasta_from_reference_allele_variant()

        # mask_fq_mq_positions_specific_to_outgroup()

        call("cp %s %s/Logs/core/" % (
            log_file_handle, os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))), logger)

    time_taken = datetime.now() - start_time_2
    if args.remove_temp:
        del_command = "rm -r %s" % temp_dir
        os.system(del_command)
