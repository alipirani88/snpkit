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

parser = argparse.ArgumentParser(description='Creating Label files individual jobs')
parser.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-filter2_only_snp_vcf_file', action='store', dest="filter2_only_snp_vcf_file",
                    help='Names of filter2 only SNP vcf file')
args = parser.parse_args()

"""Set variables"""
dir = args.filter2_only_snp_vcf_dir
unique_positions_file = dir + "unique_positions_file"

""" Generate unique positions array"""
position_array_sort = []
f = open(unique_positions_file, 'r+')
for line in f:
    line = line.strip()
    position_array_sort.append(line)
f.close()
""" End: Generate unique positions array"""

file = args.filter2_only_snp_vcf_file
print "Processing %s" % file
out_file_name = file + "_positions_label"

#Changed 21 sept
""" Generate proximate, unmapped, variant positions array"""
current_unmapped_file = out_file_name.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "unmapped.bed_positions")
current_proximate_file = out_file_name.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "filter2_final.vcf_no_proximate_snp.vcf_positions_array")
current_variant_position_file = out_file_name.replace("filter2_final.vcf_no_proximate_snp.vcf_positions_label", "filter2_final.vcf_no_proximate_snp.vcf")

#initialize array
array_name = os.path.basename(out_file_name)

#variant position array
variant_position_array = "variant_" + str(array_name)
variant_position_array = []
with open(current_variant_position_file, 'rU') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        position = row[0]
        if not position.startswith('#'):
            variant_position_array.append(row[1])
csv_file.close()

#unmapped position array
unmapped_array = "unmapped_" + str(array_name)
unmapped_array = []
with open(current_unmapped_file, 'rU') as fp1:
    for line in fp1:
        line = line.strip()
        unmapped_array.append(line)
fp1.close()

#proximate position array
proximate_array = "proximate_" + str(array_name)
proximate_array = []
with open(current_proximate_file, 'rU') as fp2:
    for liness in fp2:
        liness = liness.strip()
        proximate_array.append(liness)
fp2.close()

################# Generate label files and check why the positions were filtered out from the final vcf file #########
f1=open(out_file_name, 'w+')
for j in position_array_sort:
    cmd = "grep -v \'^#\' %s | awk -F\'\t\' \'{print $2}\' | grep -w \'%s\'" % (file, j)
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    if not out:
        mpileup_file = file.replace("filter2_final.vcf_no_proximate_snp.vcf", "aln_mpileup_raw.vcf_5bp_indel_removed.vcf")
        final_file = mpileup_file
        cmd2 = "grep -P \'\s+" + j + "\s+\' " + final_file
        proc = subprocess.Popen([cmd2], stdout=subprocess.PIPE, shell=True)
        (out2, err2) = proc.communicate()
        if not out2:
            if j in unmapped_array:
                st = "reference_unmapped_position\n"
                f1.write(st)
            else:
                st = "reference_allele\n"
                f1.write(st)
        else:
            line_string_array = out2.split('\t')
            if line_string_array[1] in proximate_array:
                pst = "_proximate_SNP"
            else:
                pst = ""
            format_string_array = line_string_array[7].split(';')
            for i in format_string_array:
                if i.startswith('FQ='):
                    m = re.match("FQ=(\W*\d+)", i)
                    FQ_value = m.group(1)
                    if FQ_value.startswith('-'):
                        st = "HighFQ"
                        #st = "high_qualty_filtered_out\n"
                        qual_string = line_string_array[5]
                        if float(qual_string) < 100.00:
                            st = st + "_QUAL"
                        for i in format_string_array:
                            if i.startswith('DP='):
                                m = re.match("DP=(\W*\d+)", i)
                                DP_value = m.group(1)
                                if int(DP_value) <= 15:
                                    st = st + "_DP"
                    else:
                        #st = "Uncertain_low_quality\n"
                        st = "LowFQ"
                        qual_string = line_string_array[5]
                        if float(qual_string) < 100.00:
                            st = st + "_QUAL"
                        for i in format_string_array:
                            if i.startswith('DP='):
                                m = re.match("DP=(\W*\d+)", i)
                                DP_value = m.group(1)
                                if int(DP_value) < 15:
                                    st = st + "_DP"
            st = st + pst + "\n"
            f1.write(st)
    else:
        st = "VARIANT" + "\n"
        f1.write(st)
f1.close()
############# End: Generate label files and check why the positions were filtered out from the final vcf file #########
