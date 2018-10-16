from __future__ import division
import sys
import os
import errno
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *

# """ Sanity Check Methods"""
# def make_sure_path_exists(out_path):
#     """
#     Fuction to make sure output folder exists. If not, create it.
#     :param: out_path
#     :return: null/exception
#     """
#     try:
#         os.makedirs(out_path)
#     except OSError as exception:
#         if exception.errno != errno.EEXIST:
#             keep_logging('\nErrors in output folder path! please change the output path or analysis name\n',
#                          '\nErrors in output folder path! please change the output path or analysis name\n', logger,
#                          'info')
#             exit()

def make_sure_files_exists(vcf_file_array, Config, logger):
    """
    Function to make sure the variant call output files exists and are not empty.
    :param: vcf_file_array
    :return: null or exception
    """
    not_found_files = []
    for files in vcf_file_array:
        ori_unmapped_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                          "unmapped.bed_positions")
        ori_proximate_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                           "filter2_final.vcf_no_proximate_snp.vcf_positions_array")
        ori_variant_position_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                                  "filter2_final.vcf_no_proximate_snp.vcf")
        ori_5bp_mpileup_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                             "aln_mpileup_raw.vcf_5bp_indel_removed.vcf")
        ori_mpileup_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                         "aln_mpileup_raw.vcf")
        ori_filter_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                        "filter2_final.vcf")
        ori_indel_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                       "filter2_indel_final.vcf")

        if os.path.isfile(ori_unmapped_file) and os.path.isfile(ori_proximate_file) and os.path.isfile(
                ori_variant_position_file) and os.path.isfile(ori_5bp_mpileup_file) and os.path.isfile(
                ori_mpileup_file) and os.path.isfile(ori_filter_file) and os.path.isfile(ori_indel_file):
            continue
        else:
            not_found_files.append(files)
    if len(not_found_files) > 0:
        for i in not_found_files:
            keep_logging('Error finding variant calling output files for: %s' % os.path.basename(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                         'Error finding variant calling output files for: %s' % os.path.basename(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), logger, 'exception')
        exit()

def make_sure_label_files_exists(vcf_file_array, uniq_snp_positions, uniq_indel_positions, Config, logger):
    """
    Function to make sure the variant call output files exists and are not empty.
    :param: vcf_file_array
    :return: null or exception
    """
    not_found_files = []
    found_incomplete = []
    for files in vcf_file_array:

        snps_label_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                                  "filter2_final.vcf_no_proximate_snp.vcf_positions_label")

        indel_label_file = files.replace("filter2_final.vcf_no_proximate_snp.vcf",
                                        "filter2_indel_final.vcf_indel_positions_label")

        num_snps_label_lines = sum(1 for line in open('%s' % snps_label_file))
        num_indel_label_lines = sum(1 for line in open('%s' % indel_label_file))

        if os.path.isfile(snps_label_file) and os.path.isfile(indel_label_file):
            if num_snps_label_lines == uniq_snp_positions and num_indel_label_lines == uniq_indel_positions:
                continue
            else:
                found_incomplete.append(files)
        else:
            not_found_files.append(files)
            # keep_logging('Error finding variant calling output files: %s' % files, 'Error finding variant calling output files: %s' % files, logger, 'exception')
    if len(not_found_files) > 0:
        for i in not_found_files:
            keep_logging('Error finding core_prep output files for: %s' % os.path.basename(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                         'Error finding core_prep output files for: %s' % os.path.basename(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), logger, 'exception')
        exit()

    if len(found_incomplete) > 0:
        for i in found_incomplete:
            tmp_file = os.path.basename(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', ''))
            keep_logging('core_prep step failed for: %s. The number of unique positions doesn\'t match. Rerun %s.pbs' % (tmp_file, i),
                         'core_prep step failed for: %s. The number of unique positions doesn\'t match. Rerun %s.pbs' % (tmp_file, i), logger, 'exception')
        exit()