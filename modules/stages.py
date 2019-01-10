__author__ = 'alipirani'

from modules.trim import *
import os
import gzip
import re
from modules.log_modules import keep_logging
from modules.bwa import align_bwa
from modules.samtools import *
from modules.picard import *
from modules.gatk import *
from modules.vcftools import *
from modules.remove_5_bp_snp_indel import *
from modules.qualimap import *
from modules.bowtie import *

""" Prepare ReadGroup option for BWA alignment """
def prepare_readgroup(forward_read, aligner, logger):
    keep_logging('Preparing ReadGroup Info', 'Preparing ReadGroup Info', logger, 'info')
    samplename = os.path.basename(forward_read)
    if forward_read.endswith(".gz"):
        #output = gzip.open(forward_read, 'rb')
        #firstLine = output.readline()
        #split_field = re.split(r":",firstLine)
        #id_name = split_field[1]
        #id_name = id_name.strip()
        #split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        #return split_field
        output = gzip.open(forward_read, 'rb')
        firstLine = output.readline()
        if ":" in firstLine:
            split_field = re.split(r":",firstLine)
            id_name = split_field[1].rstrip()
            id_name = id_name.rstrip()
        # if aligner == "bowtie":
        #     split_field = "--rg-id %s --rg SM:%s --rg LB:1 --rg PL:Illumina" % (split_field[1], samplename)
        # elif aligner == "bwa":
        #     split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""

        ###Pending
        elif "/" in firstLine:
            split_field = re.split(r"/",firstLine)
            id_name = split_field[1].rstrip()
            id_name = id_name.rstrip()
            #id_name = split_field[1].rstrip()
            #id_name = id_name.rstrip()
            split_field = "\"" + "@RG" + "\\tID:" + id_name + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        else:
            id_name = "1"
            split_field = "\"" + "@RG" + "\\tID:" + id_name + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        if aligner == "bowtie":
            split_field = "--rg-id %s --rg SM:%s --rg LB:1 --rg PL:Illumina" % (split_field[1], samplename)
        elif aligner == "bwa":
            split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field

    elif forward_read.endswith(".fastq"):
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field

    elif forward_read.endswith(".fq"):
        ###
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field

""" Raw data Pre-processing using Trimmomatic """
def trimmomatic(input1, input2, out_path, crop, logger, Config):
    trim(input1, input2, out_path, crop, logger, Config)

""" bwa, smalt, bowtie: Alignment """
def align(out_path, ref_index, split_field, analysis, files_to_delete, logger, Config, type):
    reference = ConfigSectionMap(ref_index, Config)['ref_path'] + "/" + ConfigSectionMap(ref_index, Config)['ref_name']
    forward_clean = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['f_p']
    reverse_clean = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['r_p']
    forward_unpaired = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['f_up']
    reverse_unpaired = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['r_up']

    aligner = ConfigSectionMap("pipeline", Config)['aligner']
    if aligner == "bwa":
        base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bwa", Config)['bwa_bin'] + "/" + ConfigSectionMap("bwa", Config)['base_cmd']
        out_file = align_bwa(base_cmd,forward_clean, reverse_clean, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type)
        return out_file
    elif aligner == "smalt":
        print "Smalt addition pending"
        exit()
        usage()
    elif aligner == "bowtie":
        parameters = ConfigSectionMap("bowtie", Config)['parameters']
        base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bowtie", Config)['bowtie_bin'] + "/" + ConfigSectionMap("bowtie", Config)['align_cmd']
        out_file = align_bowtie(base_cmd,forward_clean, reverse_clean, forward_unpaired, reverse_unpaired, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type, parameters)
        return out_file
        exit()
        usage()

""" samtools: Post-Alignment SAM/BAM conversion, sort, index """
def prepare_bam(out_sam, out_path, analysis, files_to_delete, logger, Config):
    out_bam = samtobam(out_sam, out_path, analysis, files_to_delete, logger, Config)
    out_sort_bam = sort_bam(out_bam, out_path, analysis, logger, Config)
    files_to_delete.append(out_sort_bam)
    out_marked_bam = markduplicates(out_sort_bam, out_path, analysis, files_to_delete, logger, Config)
    out_sort_bam = sort_bam(out_marked_bam, out_path, analysis, logger, Config)
    index_bam(out_sort_bam, out_path, logger, Config)
    if not os.path.isfile(out_sort_bam):
        keep_logging('Error in SAM/BAM conversion, sort, index. Exiting.', 'Error in SAM/BAM conversion, sort, index. Exiting.', logger, 'exception')
        exit()
    else:
        return out_sort_bam

""" samtools, gatk: Variant calling """
def variant_calling(out_finalbam, out_path, index, analysis, logger, Config):
    variant_caller = eval(ConfigSectionMap("pipeline", Config)['variant_caller'])
    final_raw_vcf = variant_caller(out_finalbam, out_path, index, analysis, logger, Config)
    if not os.path.isfile(final_raw_vcf):
        keep_logging('Error in Samtools Variant Calling step. Exiting.', 'Error in Samtools Variant Calling step. Exiting.', logger, 'exception')
        exit()
    else:
        return final_raw_vcf

""" Statistics Report """
def alignment_stats(out_sorted_bam, out_path, analysis, logger, Config):
    alignment_stats_file = flagstat(out_sorted_bam, out_path, analysis, logger, Config)
    keep_logging('The Alignments Stats file from Samtools: {}'.format(alignment_stats_file), 'The Alignments Stats file from Samtools: {}'.format(alignment_stats_file), logger, 'debug')
    return alignment_stats_file

def vcf_stats(final_raw_vcf, out_path, analysis, logger, Config):
    vcf_stats_file = vcfstats(final_raw_vcf, out_path, analysis, logger, Config)
    return vcf_stats_file

def qualimap(out_sorted_bam, out_path, analysis, logger, Config):
    qualimap_report = bamqc(out_sorted_bam, out_path, analysis, logger, Config)
    return qualimap_report

def picard_stats(out_sorted_bam, out_path, analysis, index, logger, Config):
    picard_stats_file = picardstats(out_sorted_bam, out_path, analysis, index, logger, Config)
    return picard_stats_file

""" Variant Filteration """
def filter_variants(final_raw_vcf, out_path, analysis, ref_index, logger, Config, Avg_dp):
    reference = ConfigSectionMap(ref_index, Config)['ref_path'] + "/" + ConfigSectionMap(ref_index, Config)['ref_name']
    gatk_filter_final_vcf_file = gatk_filter(final_raw_vcf, out_path, analysis, reference, logger, Config, Avg_dp)
    #gatk_filter_final_vcf_contamination_file = gatk_filter_contamination(final_raw_vcf, out_path, analysis, reference, logger, Config, Avg_dp)
    gatk_filter_final_vcf_file_no_proximate_snp = remove_proximate_snps(gatk_filter_final_vcf_file, out_path, analysis, reference, logger, Config)
    keep_logging('The vcf file with no proximate snp: {}'.format(gatk_filter_final_vcf_file_no_proximate_snp), 'The vcf file with no proximate snp: {}'.format(gatk_filter_final_vcf_file_no_proximate_snp), logger, 'debug')


def filter_indels(final_raw_indel_vcf, out_path, analysis, ref_index, logger, Config, Avg_dp):
    reference = ConfigSectionMap(ref_index, Config)['ref_path'] + "/" + ConfigSectionMap(ref_index, Config)['ref_name']
    gatk_filter_final_vcf_file = gatk_filter_indel(final_raw_indel_vcf, out_path, analysis, reference, logger, Config, Avg_dp)


""" Unused Methods """

""" Generate different VCF's """
def raw_only_snp_vcf(final_raw_vcf, out_path, analysis, ref_index):
    print "\n################## Generating different VCF ##################\n"
    reference = ConfigSectionMap(ref_index)['ref_path'] + "/" + ConfigSectionMap(ref_index)['ref_name']
    only_snp_raw_vcf_file = only_snp_raw_vcf(final_raw_vcf, out_path, analysis, reference)
    print "\nThe final raw vcf file(only SNP): %s" % only_snp_raw_vcf_file

## Variant Filteration ##
def filter1_variants(final_raw_vcf, out_path, analysis, ref_index):
    reference = ConfigSectionMap(ref_index)['ref_path'] + "/" + ConfigSectionMap(ref_index)['ref_name']
    gatk_filter1_final_vcf_file = gatk_filter1(final_raw_vcf, out_path, analysis, reference)
    gatk_filter1_final_vcf_file_no_proximate_snp = remove_proximate_snps(gatk_filter1_final_vcf_file, out_path, analysis, reference)
    #only_snp_filter1_vcf_file = only_snp_filter1_vcf(gatk_filter1_final_vcf_file, out_path, analysis, reference)
    #print only_snp_filter1_vcf_file
    gatk_vcf2fasta_filter1_file = gatk_vcf2fasta_filter1(gatk_filter1_final_vcf_file, out_path, analysis, reference)
    gatk_vcf2fasta_filter1_file_no_proximate = gatk_vcf2fasta_filter1(gatk_filter1_final_vcf_file_no_proximate_snp, out_path, analysis, reference)
    vcftools_vcf2fasta_filter1_file = vcftools_vcf2fasta_filter1(gatk_vcf2fasta_filter1_file, out_path, analysis, reference)
    vcftools_vcf2fasta_filter1_file_no_proximate = vcftools_vcf2fasta_filter1(gatk_filter1_final_vcf_file_no_proximate_snp, out_path, analysis, reference)
    print "\nThe final Consensus Fasta file from GATK: %s" % gatk_vcf2fasta_filter1_file
    print "\nThe final Consensus Fasta file from GATK with no proximate snps: %s" % gatk_vcf2fasta_filter1_file_no_proximate
    print "\nThe final Consensus Fasta file from VCF-consensus: %s" % vcftools_vcf2fasta_filter1_file
    print "\nThe final Consensus Fasta file from VCF-consensus with no proximate snps: %s" % vcftools_vcf2fasta_filter1_file_no_proximate

## Remove SAM files ##
def remove_files(analysis, out_path, out_sam, out_sorted_bam):
    os.remove(out_sam)
    if os.path.isfile(out_sorted_bam):
        raw_bam_file = "%s/%s_aln.bam" % (out_path, analysis)
        os.remove(raw_bam_file)

## picard, gatk: Mark Duplicates; Indel Realignment ##
def post_align_bam(out_sorted_bam, out_path, reference, analysis):
    print "\n################## Picard, GATK: Mark Duplicates; Indel Realignment. ##################\n"
    out_marked_bam = markduplicates(out_sorted_bam, out_path, analysis)
    if not os.path.isfile(out_marked_bam):
        print "Problem in Picard mark Duplicate\n"
        exit()
        usage()
    out_marked_sort_bam = sort_bam(out_marked_bam, out_path, analysis)
    out_marked_sort_bam_rename = "%s/%s_aln_sort_marked.bam" % (out_path, analysis)
    cp_cmd = "cp %s %s" % (out_marked_sort_bam, out_marked_sort_bam_rename)
    os.system(cp_cmd)
    index_bam(out_marked_sort_bam_rename, out_path)
    out_indel_realign_bam = indel_realign(out_marked_sort_bam_rename, reference, out_path, analysis)
    if not os.path.isfile(out_indel_realign_bam):
        print "Problem in Indel Realignment\n"
        exit()
        usage()
    else:
        print "\n################## END: Picard, GATK: Mark Duplicates; Indel Realignment. #############\n"
        return out_indel_realign_bam
