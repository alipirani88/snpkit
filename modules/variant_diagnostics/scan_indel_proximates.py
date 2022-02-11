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
#import pandas as pd
import errno
from pyfasta import Fasta
from datetime import datetime
import time
import threading
import json
from cyvcf2 import VCF
import ConfigParser
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
from tabix import *
from Bio import SeqIO
from phage_detection import *
from find_repeats import *
from mask_regions import *
from fasttree import fasttree
from gubbins import *
from raxml import raxml
from pyfasta import Fasta
from core_prep_sanity_checks import *
from core_prep_functions import *
from iqtree import iqtree
from memory_profiler import profile
# import numpy as np


# Parse Command line Arguments
parser = argparse.ArgumentParser(description='Scan vcf files and determine SNPs filtered with 5bp Indel proximate filter. ')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
parser.add_argument('-matrix', action='store', dest="matrix",
                    help='SNP Code Matrix to perform scan')
parser.add_argument('-vcffiles', action='store', dest="vcffiles",
                    help='VCF filenames for parsing')
args = parser.parse_args()


indel_proximate_variants = []
with open(args.vcffiles) as fp:
    for line in fp:
        line = line.strip()
        print "Reading vcf file - %s" % line
        before_indel_proximate_variants = []
        after_indel_proximate_variants = []

        for variants in VCF(line.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf')):
            if variants.POS not in before_indel_proximate_variants and variants.INFO.get('INDEL') != True:
                # print variants.INFO.get('INDEL')
                # print variants.POS
                before_indel_proximate_variants.append(variants.POS)
            # else:
            #     print variants.POS
            #     print variants.INFO.get('INDEL')
        for variants in VCF(line.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf_5bp_indel_removed.vcf')):
            if variants.POS not in after_indel_proximate_variants and variants.INFO.get('INDEL') != True:
                # print variants.INFO.get('INDEL')
                # print variants.POS
                after_indel_proximate_variants.append(variants.POS)
            # else:
            #     print variants.POS
            #     print variants.INFO.get('INDEL')
        print "Raw pre-filtered variants - %s" % len(before_indel_proximate_variants)
        print "After Indel proximate filter - %s" % len(after_indel_proximate_variants)

        set_difference = set(before_indel_proximate_variants) - set(after_indel_proximate_variants)
        #indel_proximate_variants = list(set_difference)
        indel_proximate_variants.extend(list(set_difference))
        print "Number of positions filtered with 5 bp Indel proximate filter - %s" % len(indel_proximate_variants)
        with open(line.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_proximate_indel_filtered_positions.txt'), 'w+') as fopen:
            for i in list(set_difference):
                fopen.write(str(i) + "\n")
        fopen.close()
print len(indel_proximate_variants)

# with open("before_indel_proximate_variants_sorted_uniq_test.txt") as fp:
#     for line in fp:
#         before_indel_proximate_variants.append(line.strip())
#     fp.close()
#
# with open("after_indel_proximate_variants_sorted_uniq_test.txt") as fp:
#     for line in fp:
#         after_indel_proximate_variants.append(line.strip())
#     fp.close()


# print "Raw pre-filtered variants - %s" % len(before_indel_proximate_variants)
# print "After Indel proximate filter - %s" % len(after_indel_proximate_variants)
#
# set_difference = set(before_indel_proximate_variants) - set(after_indel_proximate_variants)
#
# indel_proximate_variants = list(set_difference)
# print "Number of positions filtered with 5 bp Indel proximate filter - %s" % len(indel_proximate_variants)

print indel_proximate_variants

print "Parsing Code Matrix..."
scan_log = open("Scan_log.txt", 'w+')

with open("%s" % args.matrix, 'rU') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    next(csv_reader, None)
    for row in csv_reader:
        position = (row[0]).split(' ')
        #print position[3]
        if int(position[3]) in indel_proximate_variants:
            if "1" in set(row[1:]) or "3" in set(row[1:]):
                if "-1" in set(row[1:]) or "-2" in set(row[1:]) or "-3" in set(row[1:]) or "2" in set(row[1:]) or "-4" in set(row[1:]):
                    print "WARNING variant allele with other filter codes: This positions was filtered out by Indel proximate filter and has a variant allele in another sample - %s" % position[3]
                    scan_log.write("WARNING variant allele with other filter codes: This positions was filtered out by Indel proximate filter and has a variant allele in another sample - %s\n" % position[3])
                else:
                    print "WARNING variant allele with no other filter codes: This positions was filtered out by Indel proximate filter and has a variant allele in another sample - %s" % position[3]
                    scan_log.write("WARNING variant allele with no other filter codes: This positions was filtered out by Indel proximate filter and has a variant allele in another sample - %s\n" % position[3])
            else:
                print "WARNING no variant allele: This positions was filtered out by Indel proximate filter and doesn't have a variant allele in another sample - %s" % position[3]
                scan_log.write("WARNING no variant allele: This positions was filtered out by Indel proximate filter and doesn't have a variant allele in another sample - %s\n" % position[3])
        # else:
        #     print position[3]
    csv_file.close()
