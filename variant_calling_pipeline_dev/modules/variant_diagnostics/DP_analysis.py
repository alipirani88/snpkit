from __future__ import division
import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
import thread
import glob
import readline
import pandas as pd
import errno
from pyfasta import Fasta
from datetime import datetime
import threading


parser = argparse.ArgumentParser(description='Extract Only reference and variant positions and generate a fasta file out of it.')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
required.add_argument('-filter2_only_snp_vcf_filename', action='store', dest="filter2_only_snp_vcf_filename",
                    help='Name of filter2 only SNP vcf file')
args = parser.parse_args()







def DP_analysis():
    print "Analyzing positions that were filtered out due to Depth..."
    extract_DP_positions = "awk -F\'\\t\' \'{print $1}\' temp_Only_filtered_positions_for_closely_matrix_DP.txt | sed \'/^$/d\' > extract_DP_positions.txt"
    os.system(extract_DP_positions)


    filename_base = os.path.basename(args.filter2_only_snp_vcf_filename)
    aln_mpileup_vcf_file = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf_5bp_indel_removed.vcf')
    analysis = filename_base.replace('_filter2_final.vcf_no_proximate_snp.vcf', '')
    grep_reference_file = "grep \'^##reference\' %s" % aln_mpileup_vcf_file
    proc = subprocess.Popen([grep_reference_file], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    reference_file = out.split(':')
    #gatk_filter2_final_vcf_file = gatk_filter2(aln_mpileup_vcf_file, temp_dir, analysis, reference_file[1])
    #gatk_filter2_final_vcf_file_no_proximate_snp = remove_proximate_snps(gatk_filter2_final_vcf_file, temp_dir, analysis, reference_file[1])
    DP_values_file = "%s/%s_DP_values" % (args.filter2_only_snp_vcf_dir, analysis)
    f2=open(DP_values_file, 'w+')


    with open("%s/temp_Only_filtered_positions_for_closely_matrix_DP.txt" % args.filter2_only_snp_vcf_dir, 'rU') as csv_filess:
        csv_readerr = csv.reader(csv_filess, delimiter='\t')
        next(csv_readerr, None)
        for rows in csv_readerr:
            #print rows
            #grep_dp_field = "grep -wP \'^\S+\s+%s\s+\b\' %s | awk -F\'\\t\' \'{print $8}\' | grep -o \'DP=.*\' | sed \'s/DP=//g\' | awk -F\';\' \'{print $1}\'" % (rows[0], aln_mpileup_vcf_file)
            grep_dp_field = "grep -w \'%s\' %s" % (rows[0], aln_mpileup_vcf_file)
            awk_dp_field = "awk -F\'\t\' \'$2 == %s\' %s | awk -F\'\t\' \'{print $8}\' | awk -F\';\' \'{print $1}\' | sed \'s/DP=//g\'" % (rows[0], aln_mpileup_vcf_file)
            #print awk_dp_field
            #proc = subprocess.Popen([grep_dp_field], stdout=subprocess.PIPE, shell=True)
            #(out2, err2) = proc.communicate()
            #out_split = out.split('\n')
            #out = out.strip()
            proc = subprocess.Popen([awk_dp_field], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out2, err2) = proc.communicate()
            #print out2.strip()
            if out2:
                #print out2.strip()
                if "INDEL" in out2:
                    #print awk_dp_field
                    out2 == "NA"
                f2.write(out2.strip() + '\n')
                # if len(out_split) > 1:
                #     print out_split[0]
                # # for i in out:
                # #     print i
                # line_split = out.split('\t')
                # #print line_split
                # if line_split[1] == rows[0]:
                #     DP_field = line_split[7].split(';')
                #     DP_value = DP_field[0].replace('DP=', '')
                #print out
            else:
                f2.write("NA\n")
                #print "NA"



    #os.system(cmd) change
    #subprocess.call(["%s" % cmd], shell=True)
    #subprocess.check_call('%s' % cmd)


DP_analysis()