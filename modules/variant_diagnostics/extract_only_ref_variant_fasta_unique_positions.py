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
import itertools
from pyfasta import Fasta
from datetime import datetime
import threading
from cyvcf2 import VCF
import ConfigParser
from config_settings import ConfigSectionMap

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

def Generate_core_plus_noncore_alignment():
    # Get reference genome ID from reference fasta file
    get_reference = Fasta(args.reference)
    if len(get_reference.keys()) == 1:
        ref_id = get_reference.keys()

    # Read in the SNP Matrix file and seperate the columns.
    core_temp = os.path.dirname(os.path.dirname(os.path.dirname(args.filter2_only_snp_vcf_dir))) + "/core_temp_dir/"
    c_reader_2 = csv.reader(open('%s/SNP_matrix_allele_new.tsv' % core_temp, 'r'), delimiter='\t')
    ncol = len(next(c_reader_2))
    del c_reader_2
    c_reader = csv.reader(open('%s/SNP_matrix_allele_new.tsv' % core_temp, 'r'), delimiter='\t')
    columns = list(zip(*c_reader))
    #columns = list(itertools.izip(*c_reader))
    del c_reader
    
    # Generate an array of all the unique variant positions that were called in all the samples
    unique_position_array = []
    for i in columns[0][1:]:
        replace_string = i.split(' ')
        if replace_string[0] != "None":
            unique_position_array.append(int(replace_string[3]))
        else:
            unique_position_array.append(int(replace_string[2]))
    counts = 1
    end = ncol
    
    # Loop over each column, check if the column name matches the sample name provided with argument args.filter2_only_snp_vcf_filename
    for i in xrange(1, end, 1):
        print_string = ""
        ref_print_string = ""
        grab_vcf_filename = len(os.path.basename(args.filter2_only_snp_vcf_filename).replace('_filter2_final.vcf_no_proximate_snp.vcf', ''))
        sample_name_re = columns[i][0][:grab_vcf_filename]
        if sample_name_re == os.path.basename(args.filter2_only_snp_vcf_filename).replace('_filter2_final.vcf_no_proximate_snp.vcf', '') or sample_name_re in os.path.basename(args.filter2_only_snp_vcf_filename).replace('_filter2_final.vcf_no_proximate_snp.vcf', ''):
            vcf_header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name_re
            print_string = print_string + ">%s\n" % sample_name_re
            ref_print_string = ref_print_string + ">%s\n" % sample_name_re
            #variant_allele = ''.join(columns[i][1:])
            variant_allele = ""
            for ntd in columns[i][1:]:
                if "/" in ntd or len(ntd) > 1:
                    variant_allele = variant_allele + ntd[0]
                else:
                    variant_allele = variant_allele + ntd
            print_string = print_string + str(variant_allele) + "\n"
            allele_variant_fasta = open("%s/%s_allele_variants.fa" % (args.filter2_only_snp_vcf_dir, sample_name_re), 'w+')
            allele_ref_variant_fasta = open("%s/%s_ref_allele_variants.fa" % (args.filter2_only_snp_vcf_dir, sample_name_re), 'w+')
            allele_ref_variant_vcf = open("%s/%s_ref_allele_variants.vcf" % (args.filter2_only_snp_vcf_dir, sample_name_re), 'w+')
            allele_ref_variant_vcf.write(vcf_header)
            allele_variant_fasta.write(print_string)
            allele_variant_fasta.close()
            variant_allele_array = []
            variant_allele_array_dict = {}

            count_index = 0
            end_index = len(unique_position_array) + 1
            for start_count in xrange(1, end_index, 1):
                pos = columns[0][start_count]
                get_positions_string = pos.split(' ')
                if get_positions_string[0] != "None":
                    get_positions = int(get_positions_string[3])
                else:
                    get_positions = int(get_positions_string[2])

                variant_allele_array_dict[get_positions] = columns[i][start_count]
            get_sample_reference = Fasta("%s/%s_allele_variants.fa" % (args.filter2_only_snp_vcf_dir, sample_name_re))
            if len(get_sample_reference.keys()) == 1:
                sample_ref_id = get_sample_reference.keys()
            for positions in unique_position_array:
                if "/" in str(variant_allele_array_dict[positions]) or len(variant_allele_array_dict[positions]) > 1:
                    allele_var = str(variant_allele_array_dict[positions][0])
                else:
                    allele_var = str(variant_allele_array_dict[positions])
                ref_allele = str(get_reference.sequence({'chr': str(get_reference.keys()[0]), 'start': int(positions), 'stop': int(positions)}))
                generate_vcf_string = "%s\t%s\t.\t%s\t%s\t221.999\t.\t.\t.\t.\n" % (ref_id[0].split(' ')[0], positions, ref_allele, allele_var)
                allele_ref_variant_vcf.write(generate_vcf_string)
            allele_ref_variant_vcf.close()
            
            filename = "%s/consensus_ref_allele_variant.sh" % args.filter2_only_snp_vcf_dir
            vcf_filename = "%s/%s_ref_allele_variants.vcf" % ((args.filter2_only_snp_vcf_dir).replace('core_temp_dir', '%s/%s_vcf_results' % (sample_name_re, sample_name_re)), sample_name_re)
            
            f1 = open(filename, 'a+')
            bgzip_cmd = "bgzip -f %s\n" % (vcf_filename)
            f1.write(bgzip_cmd)
            subprocess.call([bgzip_cmd], shell=True)
            tabix_cmd = "tabix -f -p vcf %s.gz\n" % (vcf_filename)
            f1.write(tabix_cmd)
            subprocess.call([tabix_cmd], shell=True)
            base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
            fasta_cmd = "cat %s | vcf-consensus %s.gz > %s_ref_allele_variants.fa\n" % (args.reference, vcf_filename, sample_name_re)
            f1.write(fasta_cmd)
            subprocess.call([fasta_cmd], shell=True)

            sed_command = "sed -i 's/>.*/>%s/g' %s_ref_allele_variants.fa\n" % (sample_name_re, sample_name_re)
            subprocess.call([sed_command], shell=True)
            f1.write(sed_command)

            unmapped_positions_file = "%s/%s_unmapped.bed_positions" % ((args.filter2_only_snp_vcf_dir).replace('core_temp_dir', '%s/%s_vcf_results' % (sample_name_re, sample_name_re)), sample_name_re)
            unmapped_vcf_file = "%s/%s_unmapped.vcf" % ((args.filter2_only_snp_vcf_dir).replace('core_temp_dir', '%s/%s_vcf_results' % (sample_name_re, sample_name_re)), sample_name_re)
            unmapped_vcf = open(unmapped_vcf_file, 'w+')
            unmapped_vcf.write(vcf_header)
            
            with open(unmapped_positions_file, 'r') as fpp:
                for lines in fpp:
                    lines = lines.strip()
                    ref_allele = str(get_reference.sequence(
                        {'chr': str(get_reference.keys()[0]), 'start': int(lines), 'stop': int(lines)}))
                    generate_vcf_string_unmapped = "%s\t%s\t.\t%s\t-\t221.999\t.\t.\t.\t.\n" % (
                    ref_id[0].split(' ')[0], lines, ref_allele)
                    unmapped_vcf.write(generate_vcf_string_unmapped)
            unmapped_vcf.close()
            
            bgzip_cmd = "bgzip -f %s\n" % (unmapped_vcf_file)
            
            tabix_cmd = "tabix -f -p vcf %s.gz\n" % (unmapped_vcf_file)
            
            subprocess.call([bgzip_cmd], shell=True)
            subprocess.call([tabix_cmd], shell=True)

            vcf_filename_unmapped = "%s/%s_ref_allele_unmapped.vcf" % ((args.filter2_only_snp_vcf_dir).replace('core_temp_dir', '%s/%s_vcf_results' % (sample_name_re, sample_name_re)), sample_name_re)
            
            bcftools_merge_cmd =  "bcftools merge --merge snps --force-samples %s.gz %s.gz -O v -o %s" % (unmapped_vcf_file, vcf_filename, vcf_filename_unmapped)

            bgzip_cmd = "bgzip -f %s\n" % (vcf_filename_unmapped)

            subprocess.call([bcftools_merge_cmd], shell=True)

            tabix_cmd = "tabix -f -p vcf %s.gz\n" % (vcf_filename_unmapped)

            fasta_cmd = "cat %s | vcf-consensus %s.gz > %s_ref_allele_unmapped_variants.fa\n" % (
                args.reference, vcf_filename_unmapped, sample_name_re)

            filename = "%s/%s_consensus_ref_allele_unmapped_variant.sh" % (args.filter2_only_snp_vcf_dir, sample_name_re)
            f1 = open(filename, 'w+')
            f1.write(bgzip_cmd)
            f1.write(tabix_cmd)
            f1.write(fasta_cmd)
            
            subprocess.call(['pwd'], shell=True)
            subprocess.call(bgzip_cmd, shell=True)
            subprocess.call(tabix_cmd, shell=True)
            subprocess.call(fasta_cmd, shell=True)
            
            sed_command = "sed -i 's/>.*/>%s/g' %s_ref_allele_unmapped_variants.fa\n" % (sample_name_re, sample_name_re)
            subprocess.call([sed_command], shell=True)
            f1.write(sed_command)
            f1.close()

        else:
            print "- Sample name %s does not match with column name %s" % (os.path.basename(args.filter2_only_snp_vcf_filename).replace('_filter2_final.vcf_no_proximate_snp.vcf', ''), sample_name_re)

Generate_core_plus_noncore_alignment()
