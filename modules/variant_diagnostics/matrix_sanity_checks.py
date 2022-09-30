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

parser = argparse.ArgumentParser(description='Sanity Check SNP matrix file')
parser.add_argument('-matrix', action='store', dest="matrix",
                    help='SNP allele Matrix to perform sanity checks.')
parser.add_argument('-functional_annotation', action='store', dest="functional_annotation",
                    help='Functional Annotation Positions file.')
args = parser.parse_args()

functional_annotation_positions = []
with open(args.functional_annotation) as fp:
    for line in fp:
        line = line.strip()
        functional_annotation_positions.append(int(line))
    fp.close()
f_handle=open("matrix_sanity_check.log.txt", 'w+')

print "- Parsing Matrix..."
N_string = ["N"]
count = 0
with open("%s" % args.matrix, 'rU') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    next(csv_reader, None)
    for row in csv_reader:
        position = (row[0]).split(' ')
        if set(row[1:]) == set(N_string):
            #print "%s\t%s" % (set(row[1:]), set(N_string))
            if int(position[3]) in functional_annotation_positions:
                f_handle.write("Functional Position %s Masked in all samples\n" % int(position[3]))
            else:
                count = count + 1
                print "- Error - Wrong position masked - %s" % int(position[3])
                f_handle.write("- Error - Wrong position masked - %s" % int(position[3]))

print "- No. of wrongly masked variants %s" % count
exit()