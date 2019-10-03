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
from tabix import *
from Bio import SeqIO
from core_prep_sanity_checks import *
from PBS_generate_jobs import *
from core_pipeline_core_prep_main import *

def core_prep_label(vcf_filenames, filter2_only_snp_vcf_dir, outgroup, reference, log_unique_time, log_file_handle, logger, jobrun, Config):
    # Create temporary Directory core_temp_dir/temp for storing temporary intermediate files. Check if core_temp_dir contains all the required files to run these pipeline.
    global temp_dir
    temp_dir = filter2_only_snp_vcf_dir + "/temp"

    # # Extract All the unique SNO and Indel position list from final filtered *_no_proximate_snp.vcf files.
    unique_position_file = create_positions_filestep(vcf_filenames, filter2_only_snp_vcf_dir, outgroup, logger)
    unique_indel_position_file = create_indel_positions_filestep(vcf_filenames, filter2_only_snp_vcf_dir, outgroup, logger)

    # bgzip and tabix all the vcf files in core_temp_dir.
    files_for_tabix = glob.glob("%s/*.vcf" % filter2_only_snp_vcf_dir)
    tabix(files_for_tabix, "vcf", logger, Config)

    # Get the cluster option; create and run jobs based on given parameter. The jobs will parse all the intermediate vcf file to extract information such as if any unique variant position was unmapped in a sample, if it was filtered out dur to DP,MQ, FQ, proximity to indel, proximity to other SNPs and other variant filter parameters set in config file.
    tmp_dir = "/tmp/temp_%s/" % log_unique_time

    create_job(filter2_only_snp_vcf_dir, jobrun, vcf_filenames, unique_position_file, tmp_dir, Config)

    create_indel_job(filter2_only_snp_vcf_dir, jobrun, vcf_filenames, unique_indel_position_file, tmp_dir, Config)

    # If Phaster Summary file doesn't exist in reference genome folder
    if not os.path.isfile("%s/summary.txt" % os.path.dirname(reference)):
        if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
            keep_logging('Functional class filter is set to yes. Preparing Functional class filters\n',
                         'Functional class filter is set to yes. Preparing Functional class filters\n', logger,
                         'info')
            if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
                # Submit Phaster jobs to find ProPhage region in reference genome.
                run_phaster(reference, filter2_only_snp_vcf_dir, logger, Config)

    call(
        "cp %s %s/Logs/core_prep/" % (log_file_handle, os.path.dirname(os.path.dirname(filter2_only_snp_vcf_dir))),
        logger)


"""core_prep methods 

    This block contains methods that are respnsible for running the first part of core_All step of the pipeline.
    This methods generates all the necessary intermediate files required for the second part of core_All step.
    Example of intermediate files: various diagnostics files/matrices where it decides why a variant was filtered out.

"""

def create_positions_filestep(vcf_filenames, filter2_only_snp_vcf_dir, outgroup, logger):

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
    if outgroup is not None:
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
        #print len(position_array_sort_excluding_outgroup)
        outgroup_specific_positions = []
        f_outgroup = open("%s/outgroup_specific_positions.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        for i in outgroup_position_array:
            if i not in position_array_sort_excluding_outgroup:
                f_outgroup.write(str(i) + '\n')
                outgroup_specific_positions.append(int(i))
                # outgroup_indel_specific_positions.append(int(i))
        f_outgroup.close()
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
        # Check why python sorting is not working
        keep_logging('Sorting unique variant positions.\n', 'Sorting unique variant positions.\n', logger, 'info')
        position_array_unique = set(position_array)
        position_array_sort = sorted(position_array_unique)
        keep_logging('\nThe number of unique variant positions:%s' % len(position_array_sort), '\nThe number of unique variant positions:%s' % len(position_array_sort), logger, 'info')
        unique_position_file = "%s/unique_positions_file" % filter2_only_snp_vcf_dir
        f=open(unique_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()

        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
            exit()

        return unique_position_file

    else:

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
        unique_position_file = "%s/unique_positions_file" % filter2_only_snp_vcf_dir
        f=open(unique_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()

        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
            exit()
        return unique_position_file

def create_indel_positions_filestep(vcf_filenames, filter2_only_snp_vcf_dir, outgroup, logger):

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
    if outgroup is not None:
        outgroup_position_indel_file_name = temp_dir + "/" + outgroup_indel_vcf_filename + "_positions"
        print outgroup_position_indel_file_name
        outgroup_position_indel_array = []
        f1 = open(outgroup_position_indel_file_name, 'r+')
        for lines in f1:
            lines = lines.strip()
            outgroup_position_indel_array.append(int(lines))
        f1.close()
        #print len(outgroup_position_indel_array)

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
        keep_logging('\nThe number of unique indel positions:%s' % len(position_array_sort), '\nThe number of unique indel positions:%s' % len(position_array_sort), logger, 'info')
        unique_indel_position_file = "%s/unique_indel_positions_file" % filter2_only_snp_vcf_dir
        f=open(unique_indel_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()
        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
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
        keep_logging('\nThe number of unique indel positions:%s' % len(position_array_sort), '\nThe number of unique indel positions:%s' % len(position_array_sort), logger, 'info')
        unique_indel_position_file = "%s/unique_indel_positions_file" % filter2_only_snp_vcf_dir
        f=open(unique_indel_position_file, 'w+')
        for i in position_array_sort:
            # Changed variable to suit sorting: 25-07-2018
            f.write(str(i) + "\n")
        f.close()
        if len(position_array_sort) == 0:
            keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
            exit()
        return unique_indel_position_file