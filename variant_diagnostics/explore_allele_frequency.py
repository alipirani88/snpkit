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
from datetime import datetime
import threading
import json
from cyvcf2 import VCF
import ConfigParser
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
from tabix import *
from memory_profiler import profile


parser = argparse.ArgumentParser(description='Creating Label files individual jobs')
parser.add_argument('-vcf', action='store', dest="vcf",
                    help='VCF file to extract allele frequency from.')
args = parser.parse_args()

def extract_allele_frequency():
    filter_passed_final_vcf = VCF((args.vcf).replace('_aln_mpileup_raw.vcf', '_filter2_final.vcf_no_proximate_snp.vcf'))
    print filter_passed_final_vcf.POSITION
    for variants in VCF("%s" % args.vcf):
        #grep -w "2872079" MRSA_CO_HA_426__aln_mpileup_raw.vcf | cut -f8 | cut -d';' -f11 | sed 's/DP4=//g' | awk -F',' '{print ($3+$4)/(($1+$2) + ($3+$4))}'
        DP4_value_list = str(variants.INFO.get('DP4')).replace('(', '').replace(')', '').split(',')
        #print DP4_value_list
        DP4_value_list = map(int, DP4_value_list)
        numerator =  DP4_value_list[2] + DP4_value_list[3]
        deno = DP4_value_list[0] + DP4_value_list[1] + DP4_value_list[2] + DP4_value_list[3]
        allele_frequency = float(numerator / deno)
        print "%s, %s, %s" % (variants.POS, variants.INFO.get('DP'), allele_frequency)
        filter_passed_final_vcf = (args.vcf).replace('_aln_mpileup_raw.vcf', '_filter2_final.vcf_no_proximate_snp.vcf')


extract_allele_frequency()