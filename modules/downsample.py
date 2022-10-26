from __future__ import division
import sys
import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict

# Parse Command line Arguments
parser = argparse.ArgumentParser(description='Downsample Raw reads data to specified Depth or default')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-coverage_depth', action='store', dest="coverage_depth",
                    help='Downsample raw data to this coverage')
required.add_argument('-read1', action='store', dest="read1",
                    help='Forward Read')
required.add_argument('-read2', action='store', dest="read2",
                    help='Reverse Read')
optional.add_argument('-genome_size', action='store', dest="genome_size",
                    help='Genome Size. If not provided, will be estimated from Mash')
optional.add_argument('-analysis', action='store', dest="analysis",
                    help='Analysis name.')


args = parser.parse_args()

#def downsample(forward_read, reverse_read, coverage_depth, genome_size):
    # Use forward read and run seqtk to extract basic information: min_len: 251; max_len: 251; avg_len: 251.00; 31 distinct quality values
    # Use total number of bases from the output file
    # calculate genome size if not provided using Mash: mash sketch -o sketch_out -k 32 -m 3 -r PCMP_H159_R1.fastq.gz
    # Calculate depth with total bases and estimated genome size if gsize is not provided

if __name__ == '__main__':

    # Set Default 
    if args.coverage_depth:
        depth = args.coverage_depth
    else:
        depth = 100

    print "Downsampling Coverage Depth to: %s" % depth

    seqtk_check = "seqtk fqchk -q3 %s > /tmp/%s_fastqchk.txt" % (args.read1, args.read1)
    print "Running: %s" % seqtk_check
    os.system(seqtk_check)
    with open("/tmp/%s_fastqchk.txt" % os.path.basename(args.read1), 'rU') as file_open:
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()
    print "Average Read Length: %s" % avg_len
    print "Total number of bases in fastq: %s" % total_bases

    if not args.genome_size:
        print "Running: mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s" % args.read1
        #os.system("mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s" % args.read1)
        mash_cmd = "mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout" % args.read1
        os.system(mash_cmd)
        with open("/tmp/sketch_stdout", 'rU') as file_open:
            for line in file_open:
                if line.startswith('Estimated genome size:'):
                    gsize = float(line.split(': ')[1].strip())
                if line.startswith('Estimated coverage:'):
                    est_cov = float(line.split(': ')[1].strip())
        file_open.close()
        print "Estimated Genome Size: %s" % gsize
        #print "Estimated Coverage: %s" % est_cov
    else:
        gsize = int(args.genome_size)

    ori_coverage_depth = int(total_bases/gsize)
    print "Original Covarage Depth: %s x" % ori_coverage_depth

    if not args.coverage_depth and ori_coverage_depth > 100:
        # Downsample to 100
        factor = float(100/ori_coverage_depth)
        print round(factor, 3)
    else:
        # Bug found.
        # factor = float(args.coverage_depth / ori_coverage_depth)
        # print round(factor, 3)
        factor = 1


    proc = subprocess.Popen(["nproc"], stdout=subprocess.PIPE, shell=True)
    (nproc, err) = proc.communicate()
    nproc = nproc.strip()


    r1_sub = "/tmp/%s" % os.path.basename(args.read1)
    r2_sub = "/tmp/%s" % os.path.basename(args.read2)

    
    os.system("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (args.read1, factor, nproc, os.path.basename(args.read1)))
    os.system("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (args.read2, factor, nproc, os.path.basename(args.read2)))
    ## downsample(forward_read, reverse_read, coverage_depth, genome_size)