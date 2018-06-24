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
from cyvcf2 import VCF
import ConfigParser
from config_settings import ConfigSectionMap
# from logging_subprocess import *
# from log_modules import *


parser = argparse.ArgumentParser(description='Extract Only reference and variant positions and generate a fasta file out of it.')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-filter2_only_snp_vcf_dir', action='store', dest="filter2_only_snp_vcf_dir",
                    help='Directory where all the filter2 only SNP vcf files are saved.')
required.add_argument('-filter2_only_snp_vcf_filename', action='store', dest="filter2_only_snp_vcf_filename",
                    help='Name of filter2 only SNP vcf file')
required.add_argument('-reference', action='store', dest="reference",
                    help='Path to Reference Fasta File')
required.add_argument('-out_core', action='store', dest="out_core",
                    help='Path to core results directory')
required.add_argument('-config', action='store', dest="config",
                    help='Path to core results directory')


args = parser.parse_args()

if args.config:
    config_file = args.config
else:
    config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
global Config
Config = ConfigParser.ConfigParser()
Config.read(config_file)


def extract_only_ref_variant_fasta_unique_positions():
    #print "here"

    # Get reference genome ID
    get_reference = Fasta(args.reference)
    if len(get_reference.keys()) == 1:
        ref_id = get_reference.keys()


    c_reader = csv.reader(open('%s/SNP_matrix_allele_new.csv' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
    c_reader_2 = csv.reader(open('%s/SNP_matrix_allele_new.csv' % args.filter2_only_snp_vcf_dir, 'r'), delimiter='\t')
    columns = list(zip(*c_reader))
    ncol = len(next(c_reader_2))


    unique_position_array = []
    for i in columns[0][1:]:
        replace_string = i.split(' ')
        if replace_string[0] != "None":
            unique_position_array.append(int(replace_string[3]))
        else:
            unique_position_array.append(int(replace_string[2]))
    #print unique_position_array

    counts = 1
    end = ncol
    for i in xrange(1, end, 1):
        print_string = ""
        ref_print_string = ""
        sample_name = str(columns[i][0])
        sample_name_re = re.sub('_R1.fastq.gz', '', sample_name)
        sample_name_re = re.sub('_*1*.fastq.gz', '', sample_name_re)

        if sample_name_re == os.path.basename(args.filter2_only_snp_vcf_filename).replace('_filter2_final.vcf_no_proximate_snp.vcf', ''):
            vcf_header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name_re
            print_string = print_string + ">%s\n" % sample_name_re
            ref_print_string = ref_print_string + ">%s\n" % sample_name_re
            #variant_allele = ''.join(columns[i][1:])
            variant_allele = ""
            for ntd in columns[i][1:]:
                if "/" in ntd:
                    variant_allele = variant_allele + ntd[0]
                else:
                    variant_allele = variant_allele + ntd
            #print variant_allele
            print_string = print_string + str(variant_allele) + "\n"
            allele_variant_fasta = open("%s/%s_allele_variants.fa" % (args.filter2_only_snp_vcf_dir, sample_name_re), 'w+')
            allele_ref_variant_fasta = open("%s/%s_ref_allele_variants.fa" % (args.filter2_only_snp_vcf_dir, sample_name_re), 'w+')
            allele_ref_variant_vcf = open("%s/%s_ref_allele_variants.vcf" % (args.filter2_only_snp_vcf_dir, sample_name_re), 'w+')
            allele_ref_variant_vcf.write(vcf_header)
            allele_variant_fasta.write(print_string)
            allele_variant_fasta.close()
            variant_allele_array = []
            variant_allele_array.append(columns[i][1:])
            get_sample_reference = Fasta("%s/%s_allele_variants.fa" % (args.filter2_only_snp_vcf_dir, sample_name_re))
            if len(get_sample_reference.keys()) == 1:
                sample_ref_id = get_sample_reference.keys()
            for positions in unique_position_array:
                pos_index = unique_position_array.index(positions)

                if "/" in str(variant_allele_array[0][pos_index]):
                    allele_var = str(variant_allele_array[0][pos_index][0])
                    #print allele_var
                else:
                    allele_var = str(variant_allele_array[0][pos_index])
                ref_allele = str(get_reference.sequence({'chr': str(get_reference.keys()[0]), 'start': int(positions), 'stop': int(positions)}))
                generate_vcf_string = "%s\t%s\t.\t%s\t%s\t221.999\t.\t.\t.\n" % (ref_id[0].split(' ')[0], positions, ref_allele, allele_var)
                allele_ref_variant_vcf.write(generate_vcf_string)
            allele_ref_variant_vcf.close()
            filename = "%s/consensus_ref_allele_variant.sh" % args.filter2_only_snp_vcf_dir

            vcf_filename = "%s/%s_ref_allele_variants.vcf" % (args.filter2_only_snp_vcf_dir, sample_name_re)
            f1 = open(filename, 'a+')
            bgzip_cmd = "%s/%s/bgzip -f %s\n" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], vcf_filename)
            f1.write(bgzip_cmd)
            subprocess.call([bgzip_cmd], shell=True)
            tabix_cmd = "%s/%s/tabix -f -p vcf %s.gz\n" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], vcf_filename)
            f1.write(tabix_cmd)
            subprocess.call([tabix_cmd], shell=True)
            base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
            fasta_cmd = "cat %s | %s/vcf-consensus %s.gz > %s_ref_allele_variants.fa\n" % (args.reference, base_vcftools_bin, vcf_filename, sample_name_re)
            f1.write(fasta_cmd)
            subprocess.call([fasta_cmd], shell=True)

            sed_command = "sed -i 's/>.*/>%s/g' %s_ref_allele_variants.fa\n" % (sample_name_re, sample_name_re)
            subprocess.call([sed_command], shell=True)
            f1.write(sed_command)

            #os.system("bash %s" % filename)
            #sequence_lgth_cmd = "for i in %s/*.fa; do %s/%s/bioawk -c fastx \'{ print $name, length($seq) }\' < $i; done" % (args.filter2_only_snp_vcf_dir, ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bioawk", Config)['bioawk_bin'])
            #os.system(sequence_lgth_cmd)
            #call("%s" % sequence_lgth_cmd, logger)


        else:
            print "Sample name %s does not match with column name %s" % (os.path.basename(args.filter2_only_snp_vcf_filename).replace('_filter2_final.vcf_no_proximate_snp.vcf', ''), sample_name_re)



extract_only_ref_variant_fasta_unique_positions()