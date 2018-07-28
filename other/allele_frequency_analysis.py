from __future__ import division
__author__ = 'alipirani'

##Pending Changes

import sys
import os
import argparse
import errno
import subprocess

################################################### Command Line Argument Parsing ###################################################
parser = argparse.ArgumentParser(description='Allele Frequency analysis: For the given variant positions; extract DP4 information from raw vcf file and calculate reference allele fraction and variant allele fraction')
parser.add_argument('-raw_mpileup_file', action='store', dest="raw_mpileup_file", help='Path to raw_mpileup_file', required=True)
parser.add_argument('-o', action='store', dest="output_folder", help='Output Path ending with output directory name to save the results', required=True)
parser.add_argument('-position', action='store', dest="position", help='File containing positions', required=True)
args = parser.parse_args()
############################################################### End #################################################################
print "Processing: %s..." % args.raw_mpileup_file
position_array = []

with open(args.position) as fp:
        for line in fp:
            line = line.strip()
            position_array.append(int(line))
out_file = args.output_folder + args.raw_mpileup_file + "_ref_fraction"
#print position_array
f1=open(out_file, 'w+')
for i in position_array:
    cmd = "grep -w \'" + str(i) + "\' " + args.raw_mpileup_file
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    if not out:
        ref_fraction = "REF"
	    f1.write(str(ref_fraction) + "\n")
    else:
        line_array = out.split('\t')
        #print line_array[1]
        format_string_array = line_array[7].split(';')
        #print format_string_array
        for i in format_string_array:
            if i.startswith('DP4='):
                allele_freq = i.split(',')
                forward_ref = allele_freq[0].replace('DP4=', '')
                reverse_ref = allele_freq[1]
                forward_nonref = allele_freq[2]
                reverse_nonref = allele_freq[3]
                total = int(forward_ref) + int(reverse_ref) + int(forward_nonref) + int(reverse_nonref)
                #print total
                ref_alleles = (int(forward_ref) + int(reverse_ref))
		        ref_fraction = ref_alleles / total
                #print ref_alleles
                #print total
                f1.write(str(ref_fraction) + "\n")
                #print ref_fraction













#with open(args.raw_mpileup_file, 'rU') as csv_file2:
        # f1=open(out_file, 'w+')
        # for line in csv_file2:
         #    if line.startswith('gi') or line.startswith('MRSA_8058'):
         #        #print line
         #        line_array = line.split('\t')
		# #print line_array[1]
         #        if int(line_array[1]) in position_array:
         #            print line_array[1]
         #            format_string_array = line_array[7].split(';')
         #            print format_string_array
         #            for i in format_string_array:
         #                if i.startswith('DP4='):
         #                    allele_freq = i.split(',')
         #                    forward_ref = allele_freq[0].replace('DP4=', '')
         #                    reverse_ref = allele_freq[1]
         #                    forward_nonref = allele_freq[2]
         #                    reverse_nonref = allele_freq[3]
         #                    total = int(forward_ref) + int(reverse_ref) + int(forward_nonref) + int(reverse_nonref)
		# 	                #print total
         #                    ref_alleles = (int(forward_ref) + int(reverse_ref))
         #                    print ref_alleles
         #                    print total
         #                    f1.write(str(ref_fraction) + "\n")
         #                    #print ref_fraction
