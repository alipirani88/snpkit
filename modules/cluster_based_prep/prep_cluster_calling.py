import sys
import os
from os import path
sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
import argparse
import errno
import glob
import re
import multiprocessing
import subprocess as subprocess
from datetime import datetime
from joblib import Parallel, delayed
from modules.logging_subprocess import *
from modules.log_modules import *
from argparse import RawTextHelpFormatter

# Command Line Argument Parsing
def parser():
    parser = argparse.ArgumentParser(description='\nPrepare Cluster based Reference genome and config files\n', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-cluster_names', action='store', dest="cluster_names", help='Provide a file describing cluster group. Each line represents a cluster and each sample name in a cluster is seperated by a tab.', required=True)
    required.add_argument('-report', action='store', dest="report", help='comma-seperated Assembly report', required=True)
    required.add_argument('-assemblies', action='store', dest="assemblies", help='Path to final assembly fasta files', required=True)
    required.add_argument('-config', action='store', dest="index", help='path to config file', required=True)
    required.add_argument('-reads_dir', action='store', dest="reads_dir", help='path to reads directory', required=True)
    return parser




def prep_cluster_names():
    report_dict_contigs = {}
    report_dict_length = {}
    os.system("mkdir filename_meta")
    os.system("mkdir filename_ls_commands")
    with open(args.report) as fp:
        for line in fp:
            line = line.strip()
            line_split = line.split(',')
            first_part = re.sub("_S.*", "", line_split[0])
            report_dict_contigs[first_part] = line_split[1]
            report_dict_length[first_part] = line_split[2]
    print "preparing cluster group reference genomes"
    cluster_groups = {}
    count = 0
    with open(args.cluster_names) as fp:
        for line in fp:
            count = count + 1
            line = line.strip()
            line_split = line.split('\t')
            cluster_groups[count] = line_split
        fp.close()
    f = open("cluster_reference_genome.sh", 'w+')
    f_temp = open("cluster_temp.txt", 'w+')
    for key in cluster_groups.keys():
        cluster_no = key
        cluster_dir = "cluster_" + str(key) + "_size_" + str(len(cluster_groups[key]))
        os.system("mkdir %s" % cluster_dir)
        cluster_file = "cluster_" + str(key) + "_size_" + str(len(cluster_groups[key])) + "_meta.txt"
        cluster_file_names = "cluster_" + str(key) + "_size_" + str(len(cluster_groups[key])) + "_filenames.txt"
        f1 = open("%s/%s" % (cluster_dir, cluster_file), 'w+')
        f2 = open("%s/%s" % (cluster_dir, cluster_file_names), 'w+')
        for i in cluster_groups[key]:
            if i in report_dict_length.keys():
                print_string = "cluster_" + str(key) + "_size_" + str(len(cluster_groups[key])) + ",%s,%s,%s" % (i, report_dict_contigs[i], report_dict_length[i])
                f1.write(print_string + '\n')
                cmd = "basename %s/%s_*R1.fastq.gz\n" % (args.reads_dir, i)
                proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                f2.write("%s" % (out))
            else:
                print "%s: Not found" % i
        get_best_ref_cmd = "sort -t',' -k3n %s/%s | head -n1 | cut -d',' -f2\n" % (cluster_dir, cluster_file)
        f.write(get_best_ref_cmd)
        f_temp.write("cluster_" + str(key) + '\n')










# Start of Main Method/Pipeline
if __name__ == '__main__':
    # Set up logging modules and config file
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    args = parser().parse_args()
    global log_unique_time
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    global logger
    if args.cluster_names:
        prep_cluster_names()