__author__ = 'alipirani'

import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict


parser = argparse.ArgumentParser(description='Parsing All position with label file and investigating positions to determine the reason why it was filtered out from the final list')
#All raw only snp pileup files should be store in the same directory where filter2 only snp vcf files are.
parser.add_argument('-positions_file_dir', action='store', dest="positions_file_dir", help='Directory where all the filter2 only SNP vcf files are saved.')
parser.add_argument('-label_filename', action='store', dest="label_filename", help='Names of All_label_final_raw file created after running paste.sh script.')
parser.add_argument('-unique_positions', action='store', dest="unique_positions", help='Names of unique_positions_file')
args = parser.parse_args()

# def generate_label_report():
#     MyValues = []
#     cmd = "ls %s" % args.positions_file_dir
#     proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
#     (out, err) = proc.communicate()
#     with open(All_position_file, 'rU') as csv_file:
#         csv_reader = csv.reader(csv_file, delimiter='\t')
#         for row in csv_reader:
#             for i in range(1, len(row[1:])):
#                 total_variant_positions = row[1].count("1")
#             print total_variant_positions
#             #MyValues.append(row[1])
#         #print MyValues
# generate_label_report()

def generate_sed_command():
    sed_file = args.positions_file_dir + "/sed_reason.sh"
    f4=open(sed_file, 'w')
    sed_command = "sed -i 's/reference_unmapped_position/0/g' All_label_final_raw\nsed -i 's/reference_allele/1/g' All_label_final_raw\nsed -i 's/VARIANT/1/g' All_label_final_raw\nsed -i 's/LowFQ_QUAL_DP_proximate_SNP/2/g' All_label_final_raw\nsed -i 's/LowFQ_DP_QUAL_proximate_SNP/2/g' All_label_final_raw\nsed -i 's/LowFQ_QUAL_proximate_SNP/2/g' All_label_final_raw\nsed -i 's/LowFQ_DP_proximate_SNP/2/g' All_label_final_raw\nsed -i 's/LowFQ_proximate_SNP/2/g' All_label_final_raw\nsed -i 's/LowFQ_QUAL_DP/2/g' All_label_final_raw\nsed -i 's/LowFQ_DP_QUAL/2/g' All_label_final_raw\nsed -i 's/LowFQ_QUAL/2/g' All_label_final_raw\nsed -i 's/LowFQ_DP/2/g' All_label_final_raw\nsed -i 's/HighFQ_QUAL_DP_proximate_SNP/4/g' All_label_final_raw\nsed -i 's/HighFQ_DP_QUAL_proximate_SNP/4/g' All_label_final_raw\nsed -i 's/HighFQ_QUAL_proximate_SNP/4/g' All_label_final_raw\nsed -i 's/HighFQ_DP_proximate_SNP/4/g' All_label_final_raw\nsed -i 's/HighFQ_proximate_SNP/7/g' All_label_final_raw\nsed -i 's/HighFQ_QUAL_DP/3/g' All_label_final_raw\nsed -i 's/HighFQ_DP_QUAL/3/g' All_label_final_raw\nsed -i 's/HighFQ_QUAL/3/g' All_label_final_raw\nsed -i 's/HighFQ_DP/3/g' All_label_final_raw\nsed -i 's/LowFQ/5/g' All_label_final_raw\nsed -i 's/HighFQ/6/g' All_label_final_raw"
    print sed_command
    f4.write(sed_command)
    os.system(sed_command)
generate_sed_command()


All_position_file = args.label_filename
position_label = OrderedDict()
with open(All_position_file, 'rU') as csv_file:
    print "reading position file"
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        position_label[row[0]] = row[1:]

##### Filter out those position array that only contain Reference allele and True Variant
##### This is for the sake of generating heatmap so that we can reduce nonrelevant data from heatmap
def generate_heatmap_position():
    print "generate heatmap matrix"
    f1=open("Only_ref_variant_positions_for_closely", 'w+')
    f2=open("Only_ref_variant_positions_for_closely_matrix", 'w+')
    f3=open("Only_filtered_positions_for_closely_matrix", 'w+')
    for value in position_label:
        lll = ['0', '2', '3', '4', '5', '6', '7']
        ref_var = ['1', '1']
        if set(ref_var) & set(position_label[value]):
            if set(lll) & set(position_label[value]):
                #print "bakwaas"
                STRR3 = value + "\t" + str(position_label[value]) + "\n"
                f3.write(STRR3)
            else:
                strr = value + "\n"
                f1.write(strr)
                STRR2 =	value +	"\t" + str(position_label[value]) + "\n"
                f2.write(STRR2)
generate_heatmap_position()







