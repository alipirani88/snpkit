import sys
import os
import argparse
import errno
from config_settings import ConfigSectionMap
import csv
from modules.logging_subprocess import *
from modules.log_modules import keep_logging
from modules.gatk import *

""" Initialize the arrays """
indel_positions = []
indel_range_positions = []

""" Remove SNPs that are within 5 bp in proximity to an indel """
def remove_5_bp_snp_indel(raw_vcf_file, out_path, analysis, reference, logger, Config):
    #variant_caller = eval(ConfigSectionMap("pipeline", Config)['variant_caller'])
    if ConfigSectionMap("pipeline", Config)['variant_caller'] == "samtools":
        remove_snps_5_bp_snp_indel_file_name = raw_vcf_file + "_5bp_indel_removed.vcf"
        with open(raw_vcf_file, 'rU') as csv_file:
            for line in csv_file:
                if not line.startswith('#'):
                    line_array = line.split('\t')
                    if line_array[7].startswith('INDEL;'):
                         indel_positions.append(line_array[1])
            for i in indel_positions:
                lower_range = int(i) - 5
                upper_range = int(i) + 6
                for positions in range(lower_range,upper_range):
                    indel_range_positions.append(positions)
        f1=open(remove_snps_5_bp_snp_indel_file_name, 'w+')
        with open(raw_vcf_file, 'rU') as csv_file2:
            for line in csv_file2:
                if not line.startswith('#'):
                   line_array = line.split('\t')
                   if int(line_array[1]) not in indel_range_positions:
                       print_string = line
                       f1.write(print_string)
                else:
                    print_string = line
                    f1.write(print_string)
        return remove_snps_5_bp_snp_indel_file_name

    elif ConfigSectionMap("pipeline", Config)['variant_caller'] == "gatkhaplotypecaller":
        remove_snps_5_bp_snp_indel_file_name = raw_vcf_file + "_5bp_indel_removed.vcf"
        indel_file_name = raw_vcf_file + "_indel.vcf"
        base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
        cmd = "%s SelectVariants -R %s -V %s -select-type INDEL -O %s" % (
        base_cmd, reference, raw_vcf_file, indel_file_name)
        call(cmd, logger)
        # keep_logging('Running Command: [%s]' % cmd, 'Running Command: [%s]' % cmd, logger, 'info')
        with open(indel_file_name, 'rU') as csv_file:
            for line in csv_file:
                if not line.startswith('#'):
                    line_array = line.split('\t')
                    indel_positions.append(line_array[1])
            for i in indel_positions:
                lower_range = int(i) - 5
                upper_range = int(i) + 6
                for positions in range(lower_range,upper_range):
                    indel_range_positions.append(positions)
        f1 = open(remove_snps_5_bp_snp_indel_file_name, 'w+')
        with open(raw_vcf_file, 'rU') as csv_file2:
            for line in csv_file2:
                if not line.startswith('#'):
                    line_array = line.split('\t')
                    if int(line_array[1]) not in indel_range_positions:
                        print_string = line
                        f1.write(print_string)
                else:
                    print_string = line
                    f1.write(print_string)
        return remove_snps_5_bp_snp_indel_file_name

""" Extract Indels from raw vcf files """
def prepare_indel(raw_vcf_file, out_path, analysis, reference, logger, Config):

    #variant_caller = eval(ConfigSectionMap("pipeline", Config)['variant_caller'])
    if ConfigSectionMap("pipeline", Config)['variant_caller'] == "samtools":
        indel_file_name = raw_vcf_file + "_indel.vcf"
        with open(raw_vcf_file, 'rU') as csv_file:
            for line in csv_file:
                if not line.startswith('#'):
                    line_array = line.split('\t')
                    if line_array[7].startswith('INDEL;'):
                         indel_positions.append(int(line_array[1]))
        #print indel_positions
        f1=open(indel_file_name, 'w+')
        with open(raw_vcf_file, 'rU') as csv_file2:
            for line in csv_file2:
                if not line.startswith('#'):
                   line_array = line.split('\t')
                   if int(line_array[1]) in indel_positions:
                       print_string = line
                       f1.write(print_string)
                else:
                    print_string = line
                    f1.write(print_string)
        #print indel_file_name
        return indel_file_name
    elif ConfigSectionMap("pipeline", Config)['variant_caller'] == "gatkhaplotypecaller":
        indel_file_name = raw_vcf_file + "_indel.vcf"
        base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
        cmd = "%s SelectVariants -R %s -V %s -select-type INDEL -O %s" % (
            base_cmd, reference, raw_vcf_file, indel_file_name)
        call(cmd, logger)
        # keep_logging('Running Command: [%s]' % cmd, 'Running Command: [%s]' % cmd, logger, 'info')
        return indel_file_name

""" Extract Indels from raw vcf files """
def prepare_indel_gatk(out_finalbam, out_path, analysis, reference, logger, Config):
    reference_filename = ConfigSectionMap(reference, Config)['ref_path'] + "/" + ConfigSectionMap(reference, Config)['ref_name']
    final_raw_vcf = gatkhaplotypecaller(out_finalbam, out_path, reference, analysis, logger, Config)
    if not os.path.isfile(final_raw_vcf):
        keep_logging('- Error in GATK Haplotype Variant Calling step. Exiting.', '- Error in GATK Haplotype Variant Calling step. Exiting.', logger, 'exception')
        exit()
    else:
        indel_file_name = final_raw_vcf + "_indel.vcf"
        base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
        cmd = "%s SelectVariants -R %s -V %s -select-type INDEL -O %s" % (
            base_cmd, reference_filename, final_raw_vcf, indel_file_name)
        call(cmd, logger)
        # keep_logging('Running Command: [%s]' % cmd, 'Running Command: [%s]' % cmd, logger, 'info')
        return indel_file_name



    # #variant_caller = eval(ConfigSectionMap("pipeline", Config)['variant_caller'])
    # if ConfigSectionMap("pipeline", Config)['variant_caller'] == "samtools":
    #     print "Samtools: Extracting indels from raw vcf files"
    #     indel_file_name = raw_vcf_file + "_indel.vcf"
    #     with open(raw_vcf_file, 'rU') as csv_file:
    #         for line in csv_file:
    #             if not line.startswith('#'):
    #                 line_array = line.split('\t')
    #                 if line_array[7].startswith('INDEL;'):
    #                      indel_positions.append(int(line_array[1]))
    #     #print indel_positions
    #     f1=open(indel_file_name, 'w+')
    #     with open(raw_vcf_file, 'rU') as csv_file2:
    #         for line in csv_file2:
    #             if not line.startswith('#'):
    #                line_array = line.split('\t')
    #                if int(line_array[1]) in indel_positions:
    #                    print_string = line
    #                    f1.write(print_string)
    #             else:
    #                 print_string = line
    #                 f1.write(print_string)
    #     #print indel_file_name
    #     return indel_file_name
    # elif ConfigSectionMap("pipeline", Config)['variant_caller'] == "gatkhaplotypecaller":
    #     print "GATK Haplotype caller: Extracting indels from raw vcf files"
    #     indel_file_name = raw_vcf_file + "_indel.vcf"
    #     base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)[
    #         'gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    #     cmd = "java -jar %s -T SelectVariants -R %s -V %s -selectType INDEL -o %s" % (
    #         base_cmd, reference, raw_vcf_file, indel_file_name)
    #     call(cmd, logger)
    #     keep_logging('Running Command: [%s]' % cmd, 'Running Command: [%s]' % cmd, logger, 'info')
    #     return indel_file_name


# Remove SNPS that are within 10 bp in proximity to each other
def remove_proximate_snps(gatk_filter2_final_vcf_file, out_path, analysis, reference, logger, Config):
    #variant_caller = eval(ConfigSectionMap("pipeline", Config)['variant_caller'])
    if ConfigSectionMap("pipeline", Config)['variant_caller'] == "samtools":
        filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']
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
                if diff < int(ConfigSectionMap(filter_criteria, Config)['prox']):
                    if position not in remove_proximate_position_array and all_position[next_position_index] not in remove_proximate_position_array:
                        remove_proximate_position_array.append(int(position))
                        remove_proximate_position_array.append(int(all_position[next_position_index]))
        #print remove_proximate_position_array
        f1=open(gatk_filter2_final_vcf_file_no_proximate_snp, 'w+')
        with open(gatk_filter2_final_vcf_file, 'rU') as csv_file2:
            for line in csv_file2:
                if not line.startswith('#'):
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


    elif ConfigSectionMap("pipeline", Config)['variant_caller'] == "gatkhaplotypecaller":
        
        keep_logging('GATK Haplotype caller: Removing proximate SNPs', 'GATK Haplotype caller: Removing proximate SNPs', logger, 'info')

        filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']
























