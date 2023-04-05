from __future__ import division
__author__ = 'alipirani'
import os
import readline
import argparse
from itertools import islice
import subprocess
import numpy as np
import pandas as pd
import glob


# Parse Command line Arguments
parser = argparse.ArgumentParser(description='Extract Read Depth signal across the genome from GATK Depth of Coverage file')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-depth', action='store', dest="depth",
                      help='Read Depth coverage file generated by GATK Depth of Coverage.')
args = parser.parse_args()



def generate_depth_dataframe(list_of_files):
    for file in list_of_files:
        print "Analyzing file - %s" % file
        df = pd.read_csv("%s" % file, sep='\t', header=0)
        n = 100
        count = 0
        fp = open(file.replace('_depth_of_coverage', '_depth_of_coverage_window_mean.txt'), 'w+')
        fp.write('%s_Mean_1000bp\n' % os.path.basename(file))
        fp1 = open('window.txt', 'w+')
        fp1.write('Window\n')
        df = df.groupby(np.arange(len(df['Total_Depth']))//n).mean()
        for index, row in df.iterrows():
            count = count + 1
            fp1.write(str(count) + '\n')
            fp.write(str(int(row[0])) + '\n')
        fp.close()
        fp1.close()
    paste_command = "paste window.txt"
    for file in list_of_files:
        paste_command = paste_command + " " + file.replace('_depth_of_coverage', '_depth_of_coverage_window_mean.txt')
    paste_command = paste_command + " > Mean_Depth_of_coverage_matrix.tsv"
    print (paste_command)
    os.system(paste_command)

def generate_depth_dataframe_Plasmid_SUR(list_of_files):

    gene_lis

    for file in list_of_files:
        print "Analyzing file - %s" % file
        df = pd.read_csv("%s" % file, sep='\t', header=0)
        n = 100
        count = 0
        fp = open(file.replace('_depth_of_coverage', '_depth_of_coverage_window_mean.txt'), 'w+')
        fp.write('%s_Mean_1000bp\n' % os.path.basename(file))
        fp1 = open('window.txt', 'w+')
        fp1.write('Window\n')
        df = df.groupby(np.arange(len(df['Total_Depth']))//n).mean()
        for index, row in df.iterrows():
            count = count + 1
            fp1.write(str(count) + '\n')
            fp.write(str(int(row[0])) + '\n')
        fp.close()
        fp1.close()
    paste_command = "paste window.txt"
    for file in list_of_files:
        paste_command = paste_command + " " + file.replace('_depth_of_coverage', '_depth_of_coverage_window_mean.txt')
    paste_command = paste_command + " > Mean_Depth_of_coverage_matrix.tsv"
    print (paste_command)
    os.system(paste_command)

list_of_files = sorted(glob.glob(args.depth + "/*/*/*_depth_of_coverage"))
print (list_of_files)
generate_depth_dataframe(list_of_files)
