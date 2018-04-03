__author__ = 'alipirani'

import sys
import os
import argparse
import errno


################################################### Command Line Argument Parsing ###################################################
parser = argparse.ArgumentParser(description='FQ analysis: Applying modified vcf filters to retail low FQ values for further analysis')
parser.add_argument('-final_raw_vcf', action='store', dest="final_raw_vcf", help='Path to final_raw_vcf(5bp snp indel vcf file)', required=True)
parser.add_argument('-o', action='store', dest="output_folder", help='Output Path ending with output directory name to save the results', required=True)
parser.add_argument('-index', action='store', dest="index", help='Reference Index Name. Full path to reference', required=True)
args = parser.parse_args()
############################################################### End #################################################################

analysis = args.final_raw_vcf.replace("_aln_mpileup_raw.vcf_5bp_indel_removed.vcf", "")
out_path = args.output_folder
reference = args.index
final_raw_vcf = args.final_raw_vcf

def gatk_filter2(final_raw_vcf, out_path, analysis, reference):
    gatk_filter2_parameter_expression = "MQ > 50 && QUAL > 100 && DP > 15"
    gatk_filter2_command = "java -jar ~/bin/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s --filterExpression \"%s\" --filterName PASS_filter2" % (reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
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
    #print remove_proximate_position_array
    f1=open(gatk_filter2_final_vcf_file_no_proximate_snp, 'w+')
    with open(gatk_filter2_final_vcf_file, 'rU') as csv_file2:
        for line in csv_file2:
            if line.startswith('gi') or line.startswith('MRSA_8058'):
               line_array = line.split('\t')
               if int(line_array[1]) not in remove_proximate_position_array:
                   #print line_array[1]
                   #print line_array[1]
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
#remove_proximate_snps(gatk_filter2_final_vcf_file, out_path, analysis)

gatk_filter2_final_vcf_file = gatk_filter2(final_raw_vcf, out_path, analysis, reference)

remove_proximate_snps(gatk_filter2_final_vcf_file, out_path, analysis, reference)
