__author__ = 'alipirani'

import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
from cyvcf2 import VCF
import timeit
import time
import ConfigParser
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
from memory_profiler import profile

parser = argparse.ArgumentParser(description='Creating Label files individual jobs')
parser.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-filter2_only_snp_vcf_file', action='store', dest="filter2_only_snp_vcf_file",
                    help='Names of filter2 only SNP vcf file')
parser.add_argument('-unique_position_file', action='store', dest="unique_position_file",
                    help='Names of unique positions file')
parser.add_argument('-tmp_dir', action='store', dest="tmp_dir",
                    help='Names of temporary directory')
args = parser.parse_args()


"""Set variables and set up the tmp directories"""
dir = args.filter2_only_snp_vcf_dir
unique_positions_file = args.unique_position_file
os.system("mkdir %s" % args.tmp_dir)
os.system("cp -f %s %s/%s" % (args.filter2_only_snp_vcf_file, args.tmp_dir, os.path.basename(args.filter2_only_snp_vcf_file)))

""" Generate unique positions array"""
position_array_sort = []
f = open(unique_positions_file, 'r+')
for line in f:
    line = line.strip()
    position_array_sort.append(line)
f.close()

""" Prepare output label file """
file = args.tmp_dir + "/" + os.path.basename(args.filter2_only_snp_vcf_file)
print "Processing %s" % file
out_file_name = args.filter2_only_snp_vcf_file + "_positions_label"

""" Get the prefix for all the arrays """
array_name = os.path.basename(out_file_name)

#Changed 8 March
""" Generate proximate, unmapped, variant positions array"""
ori_unmapped_file = out_file_name.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "unmapped.bed_positions")
ori_proximate_file = out_file_name.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "filter2_final.vcf_no_proximate_snp.vcf_positions_array")
ori_variant_position_file = out_file_name.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "filter2_final.vcf_no_proximate_snp.vcf")
ori_mpileup_file = out_file_name.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "aln_mpileup_raw.vcf_5bp_indel_removed.vcf")

current_unmapped_file = args.tmp_dir + "/%s" % (os.path.basename(ori_unmapped_file))
current_proximate_file = args.tmp_dir + "/%s" % (os.path.basename(ori_proximate_file))
current_variant_position_file = args.tmp_dir + "/%s" % (os.path.basename(ori_variant_position_file))
current_mpileup_file = args.tmp_dir + "/%s" % (os.path.basename(ori_mpileup_file))

os.system("cp -f %s %s/" % (ori_unmapped_file, args.tmp_dir))
os.system("cp -f %s %s/" % (ori_proximate_file, args.tmp_dir))
os.system("cp -f %s %s/" % (ori_variant_position_file, args.tmp_dir))
os.system("cp -f %s %s/" % (ori_mpileup_file, args.tmp_dir))


# Optimization changes
# #variant position array
# variant_position_array = "variant_" + str(array_name)
# variant_position_array = []
# with open(current_variant_position_file, 'rU') as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter='\t')
#     for row in csv_reader:
#         position = row[0]
#         if not position.startswith('#'):
#             variant_position_array.append(row[1])
# csv_file.close()

# Optimization changes
def generate_dicts():
    #unmapped position dict
    program_starts = time.time()
    global unmapped_array
    unmapped_array = "unmapped_" + str(array_name)
    unmapped_array = {}
    with open(current_unmapped_file, 'rU') as fp1:
        for line in fp1:
            line = line.strip()
            unmapped_array[line] = ""
    fp1.close()
    now = time.time()
    print "Time taken to load unmapped positions array - {0} seconds".format(now - program_starts)

    #proximate position dict
    program_starts = time.time()
    global proximate_array
    proximate_array = "proximate_" + str(array_name)
    proximate_array = {}
    with open(current_proximate_file, 'rU') as fp2:
        for liness in fp2:
            liness = liness.strip()
            proximate_array[liness] = ""
    fp2.close()
    now = time.time()
    print "Time taken to load proximate positions array - {0} seconds".format(now - program_starts)

    """ Prepare cyvcf vcf files - Load Cyvcf objects """
    # Optimization changes
    program_starts = time.time()
    global positions_final_vcf
    global positions_mpileup_vcf
    positions_final_vcf = defaultdict(list)
    positions_mpileup_vcf = defaultdict(list)

    for variants in VCF(args.filter2_only_snp_vcf_file + ".gz"):
        positions_final_vcf[int(variants.POS)].append(variants.INFO.get('DP'))
    now = time.time()
    print "Time taken to load filtered positions array - {0} seconds".format(now - program_starts)

    program_starts = time.time()
    for variants in VCF(ori_mpileup_file + ".gz"):
        positions_mpileup_vcf[int(variants.POS)].append(variants.INFO.get('DP'))
        positions_mpileup_vcf[int(variants.POS)].append(variants.INFO.get('FQ'))
        positions_mpileup_vcf[int(variants.POS)].append(variants.QUAL)
        positions_mpileup_vcf[int(variants.POS)].append(variants.INFO.get('MQ'))
        positions_mpileup_vcf[int(variants.POS)].append(variants.INFO.get('AF1'))
    now = time.time()
    print "Time taken to load raw vcf data array - {0} seconds".format(now - program_starts)

# @profile
def get_reason():
    generate_dicts()
    #print "Time taken to generate dictionaries: %s" % (timeit.timeit(generate_dicts, number=1))
    #program_starts = time.time()
    f1=open(out_file_name, 'w+')

    # Older chunk of code -slower
    # for j in position_array_sort:
    #     """ Check if the unique position is present in the final no_proximate_snp.vcf file """
    #     if int(j) not in positions_final_vcf.keys():
    #         if int(j) not in positions_mpileup_vcf.keys():
    #             if j in unmapped_array.keys():
    #                 st = "reference_unmapped_position\n"
    #                 f1.write(st)
    #             else:
    #                 st = "reference_allele\n"
    #                 f1.write(st)
    #         else:
    #             if j in proximate_array.keys():
    #                 pst = "_proximate_SNP"
    #             else:
    #                 pst = ""
    #             if positions_mpileup_vcf[int(j)][1] < -40:
    #                 st = "HighFQ"
    #                 if positions_mpileup_vcf[int(j)][2] < 100.00:
    #                     st = st + "_QUAL"
    #                 if positions_mpileup_vcf[int(j)][0] < 15:
    #                     st = st + "_DP"
    #             else:
    #                 st = "LowFQ"
    #                 if positions_mpileup_vcf[int(j)][2] < 100.00:
    #                     st = st + "_QUAL"
    #                 if positions_mpileup_vcf[int(j)][0] < 15:
    #                     st = st + "_DP"
    #             st = st + pst + "\n"
    #             f1.write(st)
    #     else:
    #         st = "VARIANT" + "\n"
    #         f1.write(st)
    #     now = time.time()
    #     print "Time taken to iterate the loop once - {0} seconds".format(now - program_starts)

    # Newer chunk of code -faster
    for j in position_array_sort:
        """ Check if the unique position is present in the final no_proximate_snp.vcf file """
        if not positions_final_vcf.has_key(int(j)):
            if not positions_mpileup_vcf.has_key(int(j)):
                if unmapped_array.has_key(j):
                    st = "reference_unmapped_position\n"
                    f1.write(st)
                else:
                    st = "reference_allele\n"
                    f1.write(st)
            else:
                if proximate_array.has_key(j):
                    pst = "_proximate_SNP"
                else:
                    pst = ""
                if positions_mpileup_vcf[int(j)][1] < -40:
                    st = "HighFQ"
                    if positions_mpileup_vcf[int(j)][2] < 100.00:
                        st = st + "_QUAL"
                    if positions_mpileup_vcf[int(j)][0] < 15:
                        st = st + "_DP"
                else:
                    st = "LowFQ"
                    if positions_mpileup_vcf[int(j)][2] < 100.00:
                        st = st + "_QUAL"
                    if positions_mpileup_vcf[int(j)][0] < 15:
                        st = st + "_DP"
                st = st + pst + "\n"
                f1.write(st)
        else:
            st = "VARIANT" + "\n"
            f1.write(st)
        now = time.time()
        #print "Time taken to iterate the loop once - {0} seconds".format(now - program_starts)
    f1.close()

print "Time taken to execute this code block: %s" % (timeit.timeit(get_reason, number=1))

# """ Generate label files and check why the positions were filtered out from the final vcf file """
# def get_reason():
#     f1=open(out_file_name, 'w+')
#     for j in position_array_sort:
#         """ Check if the unique position is present in the final no_proximate_snp.vcf file """
#         if int(j) not in positions_final_vcf.keys():
#             if int(j) not in positions_mpileup_vcf.keys():
#                 if j in unmapped_array.keys():
#                     st = "reference_unmapped_position\n"
#                     f1.write(st)
#                 else:
#                     st = "reference_allele\n"
#                     f1.write(st)
#             else:
#                 if j in proximate_array.keys():
#                     pst = "_proximate_SNP"
#                 else:
#                     pst = ""
#                 if positions_mpileup_vcf[int(j)][1] < -40:
#                     st = "HighFQ"
#                     if positions_mpileup_vcf[int(j)][2] < 100.00:
#                         st = st + "_QUAL"
#                     if positions_mpileup_vcf[int(j)][0] < 15:
#                         st = st + "_DP"
#                 else:
#                     st = "LowFQ"
#                     if positions_mpileup_vcf[int(j)][2] < 100.00:
#                         st = st + "_QUAL"
#                     if positions_mpileup_vcf[int(j)][0] < 15:
#                         st = st + "_DP"
#                 st = st + pst + "\n"
#                 f1.write(st)
#         else:
#             st = "VARIANT" + "\n"
#             f1.write(st)
#     f1.close()