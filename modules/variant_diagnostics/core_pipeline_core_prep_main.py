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
from modules.logging_subprocess import *
from modules.log_modules import *
from modules.tabix import *
from Bio import SeqIO
from modules.core_prep_sanity_checks import *
from PBS_generate_jobs import *


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
    paste_file = args.filter2_only_snp_vcf_dir + "/paste_label_files.sh"
    f4 = open(paste_file, 'w+')
    paste_command = "paste %s/unique_positions_file" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                               '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
        paste_command = paste_command + " " + label_file
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (
    args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

    call("%s" % header_awk_cmd, logger)
    call("%s" % sed_header, logger)
    call("%s" % sed_header_2, logger)

    temp_paste_command = paste_command + " > %s/temp_label_final_raw.txt" % args.filter2_only_snp_vcf_dir
    paste_command = paste_command + " > %s/All_label_final_raw" % args.filter2_only_snp_vcf_dir
    f4.write(paste_command)
    f4.close()
    sort_All_label_cmd = "sort -n -k1,1 %s/All_label_final_raw > %s/All_label_final_sorted.txt" % (
    args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    paste_command_header = "cat %s/header.txt %s/All_label_final_sorted.txt > %s/All_label_final_sorted_header.txt" % (
    args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    ls = []
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                               '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
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
    call("%s" % sort_All_label_cmd, logger)
    call("%s" % paste_command_header, logger)

    """ Assign numeric code to each variant filter reason"""
    subprocess.call([
                        "sed -i 's/reference_unmapped_position/0/g' %s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
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
        paste_file = args.filter2_only_snp_vcf_dir + "/paste_label_files_outgroup.sh"
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

        temp_paste_command = paste_command + " > %s/temp_label_final_raw_outgroup.txt" % args.filter2_only_snp_vcf_dir
        paste_command = paste_command + " > %s/All_label_final_raw_outgroup" % args.filter2_only_snp_vcf_dir
        f4.write(paste_command)
        f4.close()
        sort_All_label_cmd = "sort -n -k1,1 %s/All_label_final_raw_outgroup > %s/All_label_final_sorted_outgroup.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
        paste_command_header = "cat %s/header_outgroup.txt %s/All_label_final_sorted_outgroup.txt > %s/All_label_final_sorted_header_outgroup.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

        ls = []
        for i in vcf_filenames:
            label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                   '_filter2_final.vcf_no_proximate_snp.vcf_positions_label')
            ls.append(label_file)
        ls.insert(0, "%s/unique_positions_file" % args.filter2_only_snp_vcf_dir)

        with open('%s/All_label_final_raw_outgroup.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
            outfile.write(paste_command)
        outfile.close()

        with open('%s/temp_label_final_raw_outgroup.txt.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile:
            outfile.write(temp_paste_command)
        outfile.close()
        call("bash %s/All_label_final_raw_outgroup.sh" % args.filter2_only_snp_vcf_dir, logger)
        call("bash %s/temp_label_final_raw_outgroup.txt.sh" % args.filter2_only_snp_vcf_dir, logger)

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
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
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

    else:
        print "Skip generating seperate intermediate files for outgroup"


def generate_indel_paste_command():
    """
    This Function will take all the *label file and generate/paste it column wise to generate a matrix. These matrix will be used in downstream analysis.
    :param: null
    :return: null
    """

    """ Paste/Generate and sort SNP Filter Label Matrix """
    paste_file = args.filter2_only_snp_vcf_dir + "/paste_indel_label_files.sh"
    f4 = open(paste_file, 'w+')
    paste_command = "paste %s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                               '_filter2_indel_final.vcf_indel_positions_label')
        paste_command = paste_command + " " + label_file
    header_awk_cmd = "awk \'{ORS=\"\t\";}{print $1}\' %s > %s/header.txt" % (
    args.filter2_only_snp_vcf_filenames, args.filter2_only_snp_vcf_dir)
    sed_header = "sed -i \'s/^/\t/\' %s/header.txt" % args.filter2_only_snp_vcf_dir
    sed_header_2 = "sed -i -e \'$a\\' %s/header.txt" % args.filter2_only_snp_vcf_dir

    # os.system(header_awk_cmd)
    # os.system(sed_header)
    # os.system(sed_header_2)

    call("%s" % header_awk_cmd, logger)
    call("%s" % sed_header, logger)
    call("%s" % sed_header_2, logger)

    temp_paste_command = paste_command + " > %s/temp_indel_label_final_raw.txt" % args.filter2_only_snp_vcf_dir
    paste_command = paste_command + " > %s/All_indel_label_final_raw" % args.filter2_only_snp_vcf_dir
    f4.write(paste_command)
    f4.close()

    call("bash %s" % paste_file, logger)

    sort_All_label_cmd = "sort -n -k1,1 %s/All_indel_label_final_raw > %s/All_indel_label_final_sorted.txt" % (
    args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
    paste_command_header = "cat %s/header.txt %s/All_indel_label_final_sorted.txt > %s/All_indel_label_final_sorted_header.txt" % (
    args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

    ls = []
    for i in vcf_filenames:
        label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                               '_filter2_indel_final.vcf_indel_positions_label')
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
    subprocess.call([
                        "sed -i 's/reference_unmapped_position/0/g' %s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir],
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
        paste_file = args.filter2_only_snp_vcf_dir + "/paste_indel_label_files_outgroup.sh"
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

        temp_paste_command = paste_command + " > %s/temp_indel_label_final_raw_outgroup.txt" % args.filter2_only_snp_vcf_dir
        paste_command = paste_command + " > %s/All_indel_label_final_raw_outgroup" % args.filter2_only_snp_vcf_dir
        f4.write(paste_command)
        f4.close()

        call("bash %s" % paste_file, logger)

        sort_All_label_cmd = "sort -n -k1,1 %s/All_indel_label_final_raw_outgroup > %s/All_indel_label_final_sorted_outgroup.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)
        paste_command_header = "cat %s/header_outgroup.txt %s/All_indel_label_final_sorted_outgroup.txt > %s/All_indel_label_final_sorted_header_outgroup.txt" % (
        args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir)

        ls = []
        for i in vcf_filenames:
            label_file = i.replace('_filter2_final.vcf_no_proximate_snp.vcf',
                                   '_filter2_indel_final.vcf_indel_positions_label')
            ls.append(label_file)
        ls.insert(0, "%s/unique_indel_positions_file" % args.filter2_only_snp_vcf_dir)

        with open('%s/All_indel_label_final_raw_outgroup.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile2:
            outfile2.write(paste_command)
        outfile2.close()

        with open('%s/temp_indel_label_final_raw_outgroup.txt.sh' % args.filter2_only_snp_vcf_dir, 'w') as outfile2:
            outfile2.write(temp_paste_command)
        outfile2.close()

        # Why is this not working?
        call("bash %s/All_indel_label_final_raw_outgroup.sh" % args.filter2_only_snp_vcf_dir, logger)
        call("bash %s/temp_indel_label_final_raw_outgroup.txt.sh" % args.filter2_only_snp_vcf_dir, logger)
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
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir],
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
    else:
        print "Skip generating seperate intermediate files for outgroup"


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
        f1 = open("%s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'w+')
        f2 = open("%s/Only_ref_variant_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f3 = open("%s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f4 = open(
            "%s/Only_filtered_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir,
            'w+')
        if args.outgroup:
            with open("%s/All_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir,
                      'rU') as csv_file:
                keep_logging(
                    'Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir,
                    'Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir,
                    logger, 'info')
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    position_label[row[0]] = row[1:]
                keep_logging('Generating different list of Positions and heatmap data matrix... \n',
                             'Generating different list of Positions and heatmap data matrix... \n', logger, 'info')
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
                            if int(value) not in outgroup_specific_positions:
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
                            if int(value) not in outgroup_specific_positions:
                                strr = value + "\n"
                                f1.write(strr)
                                STRR3 = value + "\t" + str(position_label[value]) + "\n"
                                f2.write(STRR3)
            csv_file.close()
            f1.close()
            f2.close()
            f3.close()
            f4.close()
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/1TRUE/-1/g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)

        else:
            with open("%s/All_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
                keep_logging(
                    'Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir,
                    'Reading All label positions file: %s/All_label_final_sorted_header.txt \n' % args.filter2_only_snp_vcf_dir,
                    logger, 'info')
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    position_label[row[0]] = row[1:]
                keep_logging('Generating different list of Positions and heatmap data matrix... \n',
                             'Generating different list of Positions and heatmap data matrix... \n', logger, 'info')
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
                            STRR3 = value + "\t" + str(position_label[value]) + "\n"
                            f2.write(STRR3)
            csv_file.close()
            f1.close()
            f2.close()
            f3.close()
            f4.close()
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_variant_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/1TRUE/-1/g' %s/Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)

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
        keep_logging(
            'Reading temporary label positions file: %s/temp_label_final_raw.txt \n' % args.filter2_only_snp_vcf_dir,
            'Reading temporary label positions file: %s/temp_label_final_raw.txt \n' % args.filter2_only_snp_vcf_dir,
            logger, 'info')
        lll = ['reference_unmapped_position', 'LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP',
               'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ', 'LowFQ_proximate_SNP',
               'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP',
               'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP',
               'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
        ref_var = ['reference_allele', 'VARIANT']

        if args.outgroup:
            print "here"
            with open("%s/temp_label_final_raw_outgroup.txt" % args.filter2_only_snp_vcf_dir, 'r') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    if set(ref_var) & set(row[1:]):
                        if set(lll) & set(row[1:]):
                            if int(row[0]) not in outgroup_specific_positions:

                                print_string = ""
                                for i in row[1:]:
                                    print_string = print_string + "\t" + i
                                STRR2 = row[0] + print_string + "\n"
                                f33.write(STRR2)
            csv_file.close()
            f33.close()

        else:
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
        f44 = open("%s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            keep_logging(
                'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir,
                'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir,
                logger, 'info')
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
        subprocess.call([
                            "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ/3/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_FQ.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)

        """
        Read temp_Only_filtered_positions_for_closely_matrix file and generate a matrix of positions that are being filtered just because of Dp
        """
        temp_position_label_DP = OrderedDict()
        f44 = open("%s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        with open("%s/temp_Only_filtered_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            keep_logging(
                'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir,
                'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_positions_for_closely_matrix.txt \n' % args.filter2_only_snp_vcf_dir,
                logger, 'info')
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
        subprocess.call([
                            "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_proximate_SNP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL_DP/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_QUAL/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ_DP/3/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowFQ/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighFQ/4/g' %s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)

    def barplot_stats():
        keep_logging(
            '\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n',
            '\nRead each Sample columns and calculate the percentage of each label to generate barplot statistics.\n',
            logger, 'info')
        """
        Read each Sample columns and calculate the percentage of each label to generate barplot statistics.
        This will give a visual explanation of how many positions in each samples were filtered out because of different reason
        """

        c_reader = csv.reader(
            open('%s/temp_Only_filtered_positions_for_closely_matrix.txt' % args.filter2_only_snp_vcf_dir, 'r'),
            delimiter='\t')
        columns = list(zip(*c_reader))
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
            true_variant = columns[i].count('VARIANT')
            unmapped_positions = columns[i].count('reference_unmapped_position')
            reference_allele = columns[i].count('reference_allele')
            Only_low_FQ = columns[i].count('LowFQ')
            Only_DP = columns[i].count('HighFQ_DP')
            Only_low_MQ = columns[i].count('HighFQ')
            low_FQ_other_parameters = columns[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns[i].count(
                'LowFQ_DP_QUAL_proximate_SNP') + columns[i].count('LowFQ_QUAL_proximate_SNP') + columns[i].count(
                'LowFQ_DP_proximate_SNP') + columns[i].count('LowFQ_proximate_SNP') + columns[i].count(
                'LowFQ_QUAL_DP') + columns[i].count('LowFQ_DP_QUAL') + columns[i].count('LowFQ_QUAL') + columns[
                                          i].count('LowFQ_DP')
            high_FQ_other_parameters = columns[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns[i].count(
                'HighFQ_DP_QUAL_proximate_SNP') + columns[i].count('HighFQ_QUAL_proximate_SNP') + columns[i].count(
                'HighFQ_DP_proximate_SNP') + columns[i].count('HighFQ_proximate_SNP') + columns[i].count(
                'HighFQ_QUAL_DP') + columns[i].count('HighFQ_DP_QUAL') + columns[i].count('HighFQ_QUAL')
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
                low_FQ_other_parameters_perc = float(((columns[i].count('LowFQ_QUAL_DP_proximate_SNP') + columns[
                    i].count('LowFQ_DP_QUAL_proximate_SNP') + columns[i].count('LowFQ_QUAL_proximate_SNP') + columns[
                                                           i].count('LowFQ_DP_proximate_SNP') + columns[i].count(
                    'LowFQ_proximate_SNP') + columns[i].count('LowFQ_QUAL_DP') + columns[i].count('LowFQ_DP_QUAL') +
                                                       columns[i].count('LowFQ_QUAL') + columns[i].count(
                            'LowFQ_DP')) * 100) / total)
            except ZeroDivisionError:
                low_FQ_other_parameters_perc = 0
            try:
                high_FQ_other_parameters_perc = float(((columns[i].count('HighFQ_QUAL_DP_proximate_SNP') + columns[
                    i].count('HighFQ_DP_QUAL_proximate_SNP') + columns[i].count('HighFQ_QUAL_proximate_SNP') + columns[
                                                            i].count('HighFQ_DP_proximate_SNP') + columns[i].count(
                    'HighFQ_proximate_SNP') + columns[i].count('HighFQ_QUAL_DP') + columns[i].count('HighFQ_DP_QUAL') +
                                                        columns[i].count('HighFQ_QUAL')) * 100) / total)
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
                bar_perc_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (os.path.basename(
                    vcf_filenames[filename_count].replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                                                                        unmapped_positions_perc, reference_allele_perc,
                                                                        true_variant_perc,
                                                                        Only_low_FQ_perc, Only_DP_perc,
                                                                        Only_low_MQ_perc, other_perc)
            f_bar_count.write(bar_string)
            f_bar_perc.write(bar_perc_string)
        f_bar_count.close()
        f_bar_perc.close()
        bargraph_R_script = "library(ggplot2)\nlibrary(reshape)\nx1 <- read.table(\"bargraph_percentage.txt\", header=TRUE)\nx1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))\nmdf1=melt(x1,id.vars=\"Sample\")\npdf(\"%s/%s_barplot.pdf\", width = 30, height = 30)\nggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat=\"identity\") + ylab(\"Percentage of Filtered Positions\") + xlab(\"Samples\") + theme(text = element_text(size=9)) + scale_fill_manual(name=\"Reason for filtered out positions\", values=c(\"#08306b\", \"black\", \"orange\", \"darkgrey\", \"#fdd0a2\", \"#7f2704\")) + ggtitle(\"Title Here\") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = \"horizontal\")\ndev.off()" % (
        args.filter2_only_snp_vcf_dir, os.path.basename(os.path.normpath(args.results_dir)))
        barplot_R_file = open("%s/bargraph.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        keep_logging('Run this R script to generate bargraph plot: %s/bargraph.R' % args.filter2_only_snp_vcf_dir,
                     'Run this R script to generate bargraph plot: %s/bargraph.R' % args.filter2_only_snp_vcf_dir,
                     logger, 'info')

    """ Methods Steps"""
    keep_logging('Running: Generating data matrices...', 'Running: Generating data matrices...', logger, 'info')
    generate_position_label_data_matrix_All_label()
    keep_logging('Running: Changing variables in data matrices to codes for faster processing...',
                 'Running: Changing variables in data matrices to codes for faster processing...', logger, 'info')
    temp_generate_position_label_data_matrix_All_label()
    keep_logging('Running: Generating Barplot statistics data matrices...',
                 'Running: Generating Barplot statistics data matrices...', logger, 'info')
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
        f1 = open("%s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir, 'w+')
        f2 = open("%s/Only_ref_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f3 = open("%s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f4 = open(
            "%s/Only_filtered_indel_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir,
            'w+')

        if args.outgroup:
            with open("%s/All_indel_label_final_sorted_header_outgroup.txt" % args.filter2_only_snp_vcf_dir,
                      'rU') as csv_file:
                keep_logging(
                    'Reading All label positions file: %s/All_indel_label_final_sorted_header.txt' % args.filter2_only_snp_vcf_dir,
                    'Reading All label positions file: %s/All_indel_label_final_sorted_header.txt' % args.filter2_only_snp_vcf_dir,
                    logger, 'info')
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    position_label[row[0]] = row[1:]
                keep_logging('Generating different list of Positions and heatmap data matrix...',
                             'Generating different list of Positions and heatmap data matrix...', logger, 'info')
                print_string_header = "\t"
                for i in vcf_filenames:
                    print_string_header = print_string_header + os.path.basename(i) + "\t"
                # f.write('\t' + print_string_header.strip() + '\n')
                f2.write('\t' + print_string_header.strip() + '\n')
                f3.write('\t' + print_string_header.strip() + '\n')
                f4.write('\t' + print_string_header.strip() + '\n')
                for value in position_label:
                    lll = ['0', '2', '3', '4', '5', '6', '7']
                    ref_var = ['1', '1TRUE']
                    if set(ref_var) & set(position_label[value]):
                        if set(lll) & set(position_label[value]):
                            if int(value) not in outgroup_indel_specific_positions:
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
                            if int(value) not in outgroup_indel_specific_positions:
                                strr = value + "\n"
                                f1.write(strr)
                                STRR3 = value + "\t" + str(position_label[value]) + "\n"
                                f2.write(STRR3)
            csv_file.close()
            f1.close()
            f2.close()
            f3.close()
            f4.close()
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_indel_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/1TRUE/-1/g' %s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
        else:
            with open("%s/All_indel_label_final_sorted_header.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_file:
                keep_logging(
                    'Reading All label positions file: %s/All_indel_label_final_sorted_header.txt' % args.filter2_only_snp_vcf_dir,
                    'Reading All label positions file: %s/All_indel_label_final_sorted_header.txt' % args.filter2_only_snp_vcf_dir,
                    logger, 'info')
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    position_label[row[0]] = row[1:]
                keep_logging('Generating different list of Positions and heatmap data matrix...',
                             'Generating different list of Positions and heatmap data matrix...', logger, 'info')
                print_string_header = "\t"
                for i in vcf_filenames:
                    print_string_header = print_string_header + os.path.basename(i) + "\t"
                # f.write('\t' + print_string_header.strip() + '\n')
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
                            STRR3 = value + "\t" + str(position_label[value]) + "\n"
                            f2.write(STRR3)
            csv_file.close()
            f1.close()
            f2.close()
            f3.close()
            f4.close()
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_indel_positions_for_closely" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_ref_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/Only_filtered_indel_positions_for_closely_matrix_TRUE_variants_filtered_out.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)
            subprocess.call([
                                "sed -i 's/1TRUE/-1/g' %s/Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir],
                            shell=True)

    def temp_generate_indel_position_label_data_matrix_All_label():

        """
        Read **temp_label_final_raw.txt** SNP position label data matrix for generating barplot statistics.
        """
        temp_position_label = OrderedDict()
        f33 = open("%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        print_string_header = "\t"
        if args.outgroup:
            for i in vcf_filenames:

                if "%s_filter2_final.vcf_no_proximate_snp.vcf" % outgroup not in i:
                    print_string_header = print_string_header + os.path.basename(i) + "\t"
        else:
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"

        f33.write('\t' + print_string_header.strip() + '\n')
        keep_logging(
            'Reading temporary label positions file: %s/temp_label_final_raw.txt' % args.filter2_only_snp_vcf_dir,
            'Reading temporary label positions file: %s/temp_label_final_raw.txt' % args.filter2_only_snp_vcf_dir,
            logger, 'info')
        # lll = ['reference_unmapped_position', 'LowFQ', 'LowFQ_DP', 'LowFQ_QUAL', 'LowFQ_DP_QUAL', 'LowFQ_QUAL_DP', 'HighFQ_DP', 'HighFQ_QUAL', 'HighFQ_DP_QUAL', 'HighFQ_QUAL_DP', 'HighFQ', 'LowFQ_proximate_SNP', 'LowFQ_DP_proximate_SNP', 'LowFQ_QUAL_proximate_SNP', 'LowFQ_DP_QUAL_proximate_SNP', 'LowFQ_QUAL_DP_proximate_SNP', 'HighFQ_DP_proximate_SNP', 'HighFQ_QUAL_proximate_SNP', 'HighFQ_DP_QUAL_proximate_SNP', 'HighFQ_QUAL_DP_proximate_SNP', 'HighFQ_proximate_SNP', '_proximate_SNP']
        lll = ['reference_unmapped_position', 'LowAF', 'LowAF_DP', 'LowAF_QUAL', 'LowAF_DP_QUAL', 'LowAF_QUAL_DP',
               'HighAF_DP', 'HighAF_QUAL', 'HighAF_DP_QUAL', 'HighAF_QUAL_DP', 'HighAF', 'LowAF_proximate_SNP',
               'LowAF_DP_proximate_SNP', 'LowAF_QUAL_proximate_SNP', 'LowAF_DP_QUAL_proximate_SNP',
               'LowAF_QUAL_DP_proximate_SNP', 'HighAF_DP_proximate_SNP', 'HighAF_QUAL_proximate_SNP',
               'HighAF_DP_QUAL_proximate_SNP', 'HighAF_QUAL_DP_proximate_SNP', 'HighAF_proximate_SNP', '_proximate_SNP']
        ref_var = ['reference_allele', 'VARIANT']

        if args.outgroup:
            with open("%s/temp_indel_label_final_raw_outgroup.txt" % args.filter2_only_snp_vcf_dir, 'r') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                next(csv_reader, None)
                for row in csv_reader:
                    if set(ref_var) & set(row[1:]):
                        if set(lll) & set(row[1:]):
                            if int(row[0]) not in outgroup_indel_specific_positions:
                                print_string = ""
                                for i in row[1:]:
                                    print_string = print_string + "\t" + i
                                STRR2 = row[0] + print_string + "\n"
                                f33.write(STRR2)
            csv_file.close()
            f33.close()
        else:
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
        Read temp_Only_filtered_positions_for_closely_matrix file and generate a matrix of positions that are being filtered just because of AF
        """
        temp_position_label_AF = OrderedDict()
        f44 = open("%s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir,
                   'w+')
        with open("%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            keep_logging(
                'Reading temporary Only_filtered_indel_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir,
                'Reading temporary Only_filtered_indel_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir,
                logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)

            for row in csv_reader:
                temp_position_label_AF[row[0]] = row[1:]
            print_string_header = "\t"
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"
            f44.write('\t' + print_string_header.strip() + '\n')
            for value in temp_position_label_AF:
                lll = ['LowAF']
                if set(lll) & set(temp_position_label_AF[value]):

                    print_string = ""
                    for i in temp_position_label_AF[value]:
                        print_string = print_string + "\t" + i
                    STRR2 = value + print_string + "\n"
                    f44.write(STRR2)
            f44.close()
            csv_file.close()
            f44.close()

        """
        Perform Sed on temp files. Find a faster way to do this.
        """
        subprocess.call([
                            "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF/3/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_AF.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)

        """
        Read temp_Only_filtered_positions_for_closely_matrix file and generate a matrix of positions that are being filtered just because of Dp
        """
        temp_position_label_DP = OrderedDict()
        f44 = open("%s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir,
                   'w+')
        with open("%s/temp_Only_filtered_indel_positions_for_closely_matrix.txt" % args.filter2_only_snp_vcf_dir,
                  'rU') as csv_file:
            keep_logging(
                'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir,
                'Reading temporary Only_filtered_positions label file: %s/temp_Only_filtered_indel_positions_for_closely_matrix.txt ' % args.filter2_only_snp_vcf_dir,
                logger, 'info')
            csv_reader = csv.reader(csv_file, delimiter='\t')
            next(csv_reader, None)
            for row in csv_reader:
                temp_position_label_DP[row[0]] = row[1:]
            print_string_header = "\t"
            for i in vcf_filenames:
                print_string_header = print_string_header + os.path.basename(i) + "\t"
            f44.write('\t' + print_string_header.strip() + '\n')
            for value in temp_position_label_DP:
                lll = ['HighAF_DP']
                ref_var = ['reference_allele', 'VARIANT']
                if set(lll) & set(temp_position_label_AF[value]):
                    print_string = ""
                    for i in temp_position_label_AF[value]:
                        print_string = print_string + "\t" + i
                    STRR2 = value + print_string + "\n"
                    f44.write(STRR2)
        f44.close()
        csv_file.close()

        """
        Perform Sed on temp files. Find a faster way to do this.
        """
        subprocess.call([
                            "sed -i 's/_filter2_final.vcf_no_proximate_snp.vcf//g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_unmapped_position/0/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/reference_allele/1/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/VARIANT/2/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_proximate_SNP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL_DP/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_QUAL/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF_DP/3/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/LowAF/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)
        subprocess.call([
                            "sed -i 's/HighAF/4/g' %s/temp_Only_filtered_indel_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir],
                        shell=True)

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
        print len(columns)
        keep_logging('Finished reading columns...', 'Finished reading columns...', logger, 'info')
        counts = 1

        if args.outgroup:
            end = len(vcf_filenames) + 1
            end = end - 1
        else:
            end = len(vcf_filenames) + 1
        print end

        f_bar_count = open("%s/bargraph_indel_counts.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_perc = open("%s/bargraph_indel_percentage.txt" % args.filter2_only_snp_vcf_dir, 'w+')
        f_bar_count.write(
            "Sample\tunmapped_positions\treference_allele\ttrue_variant\tOnly_low_AF\tOnly_DP\tOnly_low_MQ\tother\n")
        f_bar_perc.write(
            "Sample\tunmapped_positions_perc\ttrue_variant_perc\tOnly_low_AF_perc\tOnly_DP_perc\tOnly_low_MQ_perc\tother_perc\n")
        for i in xrange(1, end, 1):
            """ Bar Count Statistics: Variant Position Count Statistics """
            print i
            true_variant = columns[i].count('VARIANT')
            unmapped_positions = columns[i].count('reference_unmapped_position')
            reference_allele = columns[i].count('reference_allele')
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
                reference_allele_perc = float((columns[i].count('reference_allele') * 100) / total)
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
        args.filter2_only_snp_vcf_dir, os.path.basename(os.path.normpath(args.results_dir)))
        barplot_R_file = open("%s/bargraph_indel.R" % args.filter2_only_snp_vcf_dir, 'w+')
        barplot_R_file.write(bargraph_R_script)
        keep_logging('Run this R script to generate bargraph plot: %s/bargraph_indel.R' % args.filter2_only_snp_vcf_dir,
                     'Run this R script to generate bargraph plot: %s/bargraph_indel.R' % args.filter2_only_snp_vcf_dir,
                     logger, 'info')

    """ Methods Steps"""
    keep_logging('Running: Generating data matrices...', 'Running: Generating data matrices...', logger, 'info')
    # if args.outgroup:
    #     f_outgroup = open("%s/outgroup_indel_specific_positions.txt" % args.filter2_only_snp_vcf_dir, 'r+')
    #     global outgroup_indel_specific_positions
    #     outgroup_indel_specific_positions = []
    #     for i in f_outgroup:
    #         outgroup_indel_specific_positions.append(i)
    #     f_outgroup.close()
    #
    #     f_outgroup = open("%s/outgroup_specific_positions.txt" % args.filter2_only_snp_vcf_dir, 'r+')
    #     global outgroup_specific_positions
    #     outgroup_specific_positions = []
    #     for i in f_outgroup:
    #         outgroup_specific_positions.append(i)
    #     f_outgroup.close()
    # else:
    #     global outgroup_specific_positions
    #     global outgroup_indel_specific_positions
    #     outgroup_indel_specific_positions = []
    #     outgroup_specific_positions = []
    generate_indel_position_label_data_matrix_All_label()
    keep_logging('Running: Changing variables in data matrices to codes for faster processing...',
                 'Running: Changing variables in data matrices to codes for faster processing...', logger, 'info')
    temp_generate_indel_position_label_data_matrix_All_label()
    keep_logging('Running: Generating Barplot statistics data matrices...',
                 'Running: Generating Barplot statistics data matrices...', logger, 'info')
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
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
        command_file = "%s/commands_list_fasta.sh" % args.filter2_only_snp_vcf_dir
        f3 = open(command_file, 'w+')
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)
            f1 = open(job_file_name, 'w+')
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -functional_filter %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, functional_filter)
            job_file_name = "%s_fasta.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
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
        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            f1 = open(job_file_name, 'w+')
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
            job_print_string = "#PBS -N %s_fasta\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/extract_only_ref_variant_fasta_unique_positions.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -reference %s -out_core %s -config %s\n" % (
            job_name, ConfigSectionMap("scheduler", Config)['email'],
            ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'],
            ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'],
            args.filter2_only_snp_vcf_dir, i, args.reference, core_vcf_fasta_dir, config_file)
            job_file_name = "%s_ref_allele_variants_fasta.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
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
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (
            job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
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

        for i in vcf_filenames:
            job_name = os.path.basename(i)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (
            job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
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
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,mem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/DP_analysis.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s\n" % (
            job_name, args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir, i)
            job_file_name = "%s_DP.pbs" % (i)
            f1 = open(job_file_name, 'w+')
            f1.write(job_print_string)
            f1.close()
        pbs_dir = args.filter2_only_snp_vcf_dir + "/*_DP.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            f3.write("bash %s\n" % i)
        f3.close()
        # os.system("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir)
        call("bash %s/commands_list_DP.sh" % args.filter2_only_snp_vcf_dir, logger)


def generate_vcf_files():
    if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
        keep_logging(
            'Removing Variants falling in Functional filters positions file: %s\n' % functional_class_filter_positions,
            'Removing Variants falling in Functional filters positions file: %s\n' % functional_class_filter_positions,
            logger,
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

    f_file = open(
        "%s/Only_ref_variant_positions_for_closely_without_functional_filtered_positions" % args.filter2_only_snp_vcf_dir,
        'w+')
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

    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)[
        'vcftools_bin']
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
        bgzip_cmd = "%s/%s/bgzip -f %s\n" % (
        ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], file)
        f1.write(bgzip_cmd)
        subprocess.call([bgzip_cmd], shell=True)
        tabix_cmd = "%s/%s/tabix -f -p vcf %s.gz\n" % (
        ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], file)
        f1.write(tabix_cmd)
        subprocess.call([tabix_cmd], shell=True)
        fasta_cmd = "cat %s | %s/vcf-consensus %s.gz > %s.fa\n" % (
        args.reference, base_vcftools_bin, file, file.replace('_filter2_final.vcf_core.vcf', ''))
        f1.write(fasta_cmd)
        subprocess.call([fasta_cmd], shell=True)
        base = os.path.basename(file)
        header = base.replace('_filter2_final.vcf_core.vcf', '')
        sed_command = "sed -i 's/>.*/>%s/g' %s.fa\n" % (header, file.replace('_filter2_final.vcf_core.vcf', ''))
        subprocess.call([sed_command], shell=True)
        f1.write(sed_command)
    keep_logging('The consensus commands are in : %s' % filename, 'The consensus commands are in : %s' % filename,
                 logger, 'info')
    sequence_lgth_cmd = "for i in %s/*.fa; do %s/%s/bioawk -c fastx \'{ print $name, length($seq) }\' < $i; done" % (
    args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'],
    ConfigSectionMap("bioawk", Config)['bioawk_bin'])
    # os.system(sequence_lgth_cmd)
    call("%s" % sequence_lgth_cmd, logger)


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
    create_job_DP(args.jobrun, vcf_filenames)
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
