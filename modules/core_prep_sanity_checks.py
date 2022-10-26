from __future__ import division
import sys
import os
import errno
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
import glob


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

        if not os.path.isfile(ori_unmapped_file):
            not_found_files.append(ori_unmapped_file)

        elif not os.path.isfile(ori_proximate_file):
            not_found_files.append(ori_proximate_file)

        elif not os.path.isfile(ori_variant_position_file):
            not_found_files.append(ori_variant_position_file)

        elif not os.path.isfile(ori_5bp_mpileup_file):
            not_found_files.append(ori_5bp_mpileup_file)

        elif not os.path.isfile(ori_mpileup_file):
            not_found_files.append(ori_mpileup_file)

        elif not os.path.isfile(ori_filter_file):
            not_found_files.append(ori_filter_file)

        elif not os.path.isfile(ori_indel_file):
            not_found_files.append(ori_indel_file)
        else:
            continue
    
    if len(not_found_files) > 0:
        for i in not_found_files:
            keep_logging('Error finding variant calling output files for: %s' % os.path.basename(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')),
                         'Error finding variant calling output files for: %s' % os.path.basename(i.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')), logger, 'exception')
            keep_logging('File not found: %s' % i,
                         'File not found: %s' % i, logger, 'debug')
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
                print (num_indel_label_lines)
                print (uniq_indel_positions)
                print (files)
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