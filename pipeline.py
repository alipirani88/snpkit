__author__ = 'alipirani'

import sys
import os
import argparse
import errno
import glob
from datetime import datetime
import ConfigParser
from config_settings import ConfigSectionMap
if sys.version_info < (3, 2):
    import subprocess32 as sp
else:
    import subprocess as sp
from modules.stages import *
from modules.remove_5_bp_snp_indel import *
from modules.bedtools import *
from modules.gatk import gatk_DepthOfCoverage
from modules.logging_subprocess import *
from modules.log_modules import *
from argparse import RawTextHelpFormatter
from memory_profiler import profile
from modules.samclip import samclip
from modules.snpeff import variant_annotation
from modules.snpeff import indel_annotation
from modules.tabix import *

# Command Line Argument Parsing
def parser():
    parser = argparse.ArgumentParser(description='\nSNPKIT - A workflow for Microbial Variant Calling, Recombination detection and Phylogenetic tree reconstruction.\n', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-type', action='store', dest="type", help='Type of analysis: SE or PE', required=True)
    required.add_argument('-config', action='store', dest="config", help='Path to Config file', required=True)
    required.add_argument('-PE1', action='store', dest="forward_raw", help='Path to Paired End file 1', required=True)
    optional.add_argument('-PE2', action='store', dest="reverse_raw", help='Path to Paired End file 2', required=False)
    required.add_argument('-o', action='store', dest="output_folder", help='Output Path ending with output directory name to save the results', required=True)
    required.add_argument('-analysis', action='store', dest="analysis_name", help='Unique analysis name to save the results', required=True)
    required.add_argument('-index', action='store', dest="index", help='Reference Index Name. Change this argument in config file and mention the reference header name such as KP_NTUH_chr/KPNIH1/KPNIH32.', required=True)
    #optional.add_argument('-coverage_depth_stats', action='store', dest="coverage_depth_stats", help='Run Only Depth of Coverage Stats module after read mapping')
    optional.add_argument('-c', action='store', dest="croplength", help='Crop Length in case needed')
    required.add_argument('-steps', action='store', dest="steps", help='Variant Calling Steps in sequential order.\n'
                                                                     '1.   All : This will run all the steps starting from cleaning the reads to variant calling;\n'
                                                                     '2.   clean,align,post-align,varcall,filter,stats : This will also run all steps starting from cleaning to variant calling. \nYou can also run part of the pipeline by giving "align,post-align,varcall,filter,stats" which will skip the cleaning part.\nThe order is required to be sequential. Also, while skipping any of the step make sure you have results already present in your output folder.\n'
                                                                     '3.   coverage_depth_stats: Run Only Depth of Coverage Stats module after cleaning and read mapping steps')
    optional.add_argument('-cluster', action='store', dest='cluster', help='Run pipeline on cluster/parallel-local/local. Make Sure to check if the [CLUSTER] section in config file is set up correctly.')
    optional.add_argument('-clean', action="store_true", help='clean up intermediate files. Default: OFF')
    optional.add_argument('-coverage_depth', action='store', dest="coverage_depth",
                          help='Downsample Reads to this user specified depth')
    optional.add_argument('-genome_size', action='store', dest="genome_size",
                          help='Genome Size. If not provided, will be estimated from Mash')
    optional.add_argument('-downsample', action='store', dest="downsample",
                          help='yes/no: Downsample Reads data to default depth of 100X or user specified depth')
    optional.add_argument('-clip', action="store_true", help='Filter SAM file for soft and hard clipped alignments. Default: OFF')
    return parser

# Main Pipeline method
def pipeline(args, logger):
    keep_logging('- START: Pipeline', 'START: Pipeline', logger, 'info')

    """ SANITATION CHECKS """

    # Check Subroutines: Arguments, Input FASTQ files, Reference Index
    keep_logging('- Checking Dependencies...', 'Checking Dependencies', logger, 'info')

    # Reference Genome file name
    reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
    keep_logging('- Getting Reference Genome name from config file: {}'.format(reference), 'Getting Reference Genome name from config file: {}'.format(reference), logger, 'info')

    # Check if FASTQ files exists
    if args.type != "PE" and args.type != "BAM":
        reverse_raw = "None"
        file_exists(args.forward_raw, args.forward_raw, reference)
    elif args.type != "PE" and args.type != "SE":
        print "BAM type... Not Integrated... continue"
    else:
        file_exists(args.forward_raw, args.reverse_raw, reference)

    # Check Java Version
    java_check()

    """ Start the pipeline: """
    # split values provided with -steps argument and decide the starting point of pipeline
    steps_list = args.steps.split(',')

    # Check cluster parameter and set cluster variable, used for running pipeline locally or parallelly on local or on cluster
    if args.cluster:
        cluster = args.cluster
    else:
        cluster = "local"

    """ INDIVIDUAL SUBPROCESS FOR EACH PIPELINE STEPS"""
    ## 1. Pre-Processing Raw reads using Trimmomatic
    def clean():
        method_start_time = datetime.now()
        keep_logging('- Cleaning reads with Trimmomatic', 'Cleaning reads with Trimmomatic', logger, 'info')
        if args.type == "PE":
            trimmomatic(args.forward_raw, args.reverse_raw, args.output_folder, args.croplength, logger, Config)
        else:
            reverse_raw = "None"
            trimmomatic(args.forward_raw, reverse_raw, args.output_folder, args.croplength, logger, Config)
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Completed cleaning reads in {}'.format(method_time_taken), 'Time taken to complete the method - clean: {}'.format(method_time_taken), logger, 'info')

    ## 2. Stages: Alignment using BWA
    def align_reads():
        method_start_time = datetime.now()
        keep_logging('- Mapping Reads using BWA', 'Mapping Reads using BWA', logger, 'info')
        split_field = prepare_readgroup(args.forward_raw, ConfigSectionMap("pipeline", Config)['aligner'], logger)
        out_sam = align(args.output_folder, args.index, split_field, args.analysis_name, files_to_delete, logger, Config, args.type)
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Completed read mapping in {}'.format(method_time_taken), 'Time taken to complete the method - align_reads: {}'.format(method_time_taken), logger, 'info')
        return out_sam
    
    def coverage_depth_stats():
        method_start_time = datetime.now()
        gatk_DepthOfCoverage_file = gatk_DepthOfCoverage(out_sorted_bam, args.output_folder, args.analysis_name, reference, logger, Config)
        alignment_stats_file = alignment_stats(out_sorted_bam, args.output_folder, args.analysis_name, logger, Config)
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Time taken to complete the method - coverage_depth_stats: {}'.format(method_time_taken), 'Time taken to complete the method - coverage_depth_stats: {}'.format(method_time_taken), logger, 'info')
        return gatk_DepthOfCoverage_file

    ## 3. Stages: Post-Alignment using SAMTOOLS, PICARD etc
    def post_align(out_sam):
        method_start_time = datetime.now()
        keep_logging('- Processing and sorting alignments for variant calling', 'Processing and sorting alignments for variant calling', logger, 'info')
        # Add Samclip - 28 July 2020
        if args.clip:
            out_sam = samclip(out_sam, args.output_folder, args.analysis_name, reference, logger, Config)
        out_sorted_bam = prepare_bam(out_sam, args.output_folder, args.analysis_name, files_to_delete, logger, Config)
        #out_sorted_bam = "%s/%s_aln_sort.bam" % (args.output_folder, args.analysis_name)
        bedgraph_coverage(out_sorted_bam, args.output_folder, args.analysis_name, reference, logger, Config)
        only_unmapped_positions_file = bedtools(out_sorted_bam, args.output_folder, args.analysis_name, logger, Config)
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Completed processing, sorting and cleaning alignments in  {}'.format(method_time_taken), 'Time taken to complete the method - post_align: {}'.format(method_time_taken), logger, 'info')
        return out_sorted_bam

    ## 4. Stages: Variant Calling
    def varcall():
        method_start_time = datetime.now()
        keep_logging('- START: Variant Calling', 'START: Variant Calling', logger, 'info')
        caller = ConfigSectionMap("pipeline", Config)['variant_caller']
        if caller == "gatkhaplotypecaller":
            keep_logging('- START: Variant Calling using GATK haplotyper.', 'START: Variant Calling using GATK haplotyper.', logger, 'info')
            final_raw_vcf_mpileup = variant_calling(out_sorted_bam, args.output_folder, args.index, args.analysis_name, logger, Config)
            #final_raw_vcf_mpileup = "%s/%s_aln_mpileup_raw.vcf" % (args.output_folder, args.analysis_name)
            final_raw_vcf = remove_5_bp_snp_indel(final_raw_vcf_mpileup, args.output_folder, args.analysis_name, reference, logger, Config)
            final_raw_indel_vcf = prepare_indel(final_raw_vcf_mpileup, args.output_folder, args.analysis_name, reference, logger, Config)
            keep_logging('- The final raw VCF file: {}'.format(final_raw_vcf), 'The final raw VCF file: {}'.format(final_raw_vcf), logger, 'debug')
            keep_logging('- The final raw Indel VCF file: {}'.format(final_raw_indel_vcf),
                         'The final raw Indel VCF file: {}'.format(final_raw_indel_vcf), logger, 'debug')
            keep_logging('- END: Variant Calling using Samtools without post-align bam input files.', 'END: Variant Calling using Samtools without post-align bam input files.', logger, 'info')
            return final_raw_vcf, final_raw_indel_vcf

        elif caller == "samtools":
            keep_logging('- START: Variant Calling using Samtools without post-align bam input files.', 'START: Variant Calling using Samtools without post-align bam input files.', logger, 'info')
            final_raw_indel_vcf = prepare_indel_gatk(out_sorted_bam, args.output_folder, args.analysis_name, args.index, logger, Config)
            final_raw_vcf_mpileup = variant_calling(out_sorted_bam, args.output_folder, args.index, args.analysis_name, logger, Config)
            #final_raw_vcf_mpileup = "%s/%s_aln_mpileup_raw.vcf" % (args.output_folder, args.analysis_name)
            final_raw_vcf = remove_5_bp_snp_indel(final_raw_vcf_mpileup, args.output_folder, args.analysis_name, reference, logger, Config)
            # GATK indel calling integration
            #final_raw_indel_vcf = prepare_indel(final_raw_vcf_mpileup, args.output_folder, args.analysis_name, reference, logger, Config)
            keep_logging('- The final raw VCF file: {}'.format(final_raw_vcf), 'The final raw VCF file: {}'.format(final_raw_vcf), logger, 'debug')
            keep_logging('- END: Variant Calling using Samtools without post-align bam input files.', 'END: Variant Calling using Samtools without post-align bam input files.', logger, 'info')
            return final_raw_vcf, final_raw_indel_vcf
        else:
            keep_logging('- Please provide Variant Caller name in config file under the section [pipeline]. Options for Variant caller: 1. samtools 2. gatkhaplotypecaller', 'Please provide Variant Caller name in config file under the section [pipeline]. Options for Variant caller: 1. samtools 2. gatkhaplotypecaller', logger, 'info')
            exit()
        keep_logging('- END: Variant Calling', 'END: Variant Calling', logger, 'info')
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Time taken to complete the method - varcall: {}'.format(method_time_taken), 'Time taken to complete the method - varcall: {}'.format(method_time_taken), logger, 'info')

    ## 5. Stages: Variant Filteration
    def filter(gatk_depth_of_coverage_file):
        method_start_time = datetime.now()
        keep_logging('- START: Variant Filteration', 'START: Variant Filteration', logger, 'info')
        final_raw_vcf_mpileup = "%s/%s_aln_mpileup_raw.vcf" % (args.output_folder, args.analysis_name)
        #final_raw_indel_vcf = prepare_indel(final_raw_vcf_mpileup, args.output_folder, args.analysis_name, reference, logger, Config)
        if not os.path.isfile(gatk_depth_of_coverage_file):
            file_basename = os.path.basename(gatk_depth_of_coverage_file)
            keep_logging('- The input file {} does not exists. Please provide another file with full path or check the files path.\n'.format(file_basename), 'The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
            exit()
        Avg_dp_cmd = "grep \'^Total\' %s | awk -F\'\t\' \'{print $3}\'" % gatk_depth_of_coverage_file
        proc = sp.Popen([Avg_dp_cmd], stdout=sp.PIPE, shell=True)
        (out, err) = proc.communicate()
        Avg_dp = float(out)
        print "The Average Depth per reference genome base is: %s" % Avg_dp
        filter_variants(final_raw_vcf, args.output_folder, args.analysis_name, args.index, logger, Config, Avg_dp)
        final_raw_indel_vcf = final_raw_vcf_mpileup + "_indel.vcf"
        filter_indels(final_raw_indel_vcf, args.output_folder, args.analysis_name, args.index, logger, Config, Avg_dp)
        keep_logging('- END: Variant Filteration', 'END: Variant Filteration', logger, 'info')
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Time taken to complete the method - filter: {}'.format(method_time_taken), 'Time taken to complete the method - filter: {}'.format(method_time_taken), logger, 'info')

    ## 6. SNP annotation
    def annotation(vcf_file):
        logs_folder = (args.output_folder).replace(args.analysis_name, '') + "/Logs"
        vc_logs_folder = logs_folder + "/variant_calling"
        method_start_time = datetime.now()
        keep_logging('- START: Variant Annotation', 'START: Variant Annotation', logger, 'info')
        variant_annotation(vcf_file, args.index, vc_logs_folder, Config, logger)
        indel_annotation(vcf_file, args.index, vc_logs_folder, Config, logger)
        keep_logging('- END: Variant Annotation', 'END: Variant Annotation', logger, 'info')
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Time taken to complete the method - filter: {}'.format(method_time_taken), 'Time taken to complete the method - filter: {}'.format(method_time_taken), logger, 'info')
    
    def tabix_vcf():
        files_for_tabix = glob.glob("%s/*.vcf" % args.output_folder)
        tabix(files_for_tabix, "vcf", logger, Config)
    
    ## 7. Stages: Statistics
    def stats():
        method_start_time = datetime.now()
        keep_logging('- START: Generating Statistics Reports', 'START: Generating Statistics Reports', logger, 'info')
        alignment_stats_file = alignment_stats(out_sorted_bam, args.output_folder, args.analysis_name, logger, Config)
        vcf_stats_file = vcf_stats(final_raw_vcf, args.output_folder, args.analysis_name, logger, Config)
        picard_stats_file = picardstats(out_sorted_bam, args.output_folder, args.analysis_name, args.index, logger, Config)
        keep_logging('- END: Generating Statistics Reports', 'END: Generating Statistics Reports', logger, 'info')
        method_time_taken = datetime.now() - method_start_time
        keep_logging('- Time taken to complete the method - stats: {}'.format(method_time_taken), 'Time taken to complete the method - stats: {}'.format(method_time_taken), logger, 'info')

    if args.downsample == "yes":
        read1, read2 = downsample(args, logger)
        args.forward_raw = read1
        args.reverse_raw = read2

    if len(steps_list) == 1:
        if steps_list[0] == "All":
            clean()
            out_sam = align_reads()
            out_sorted_bam = post_align(out_sam)
            out_sorted_bam = "%s/%s_aln_sort.bam" % (args.output_folder, args.analysis_name)
            gatk_depth_of_coverage_file = "%s/%s_depth_of_coverage.sample_summary" % (args.output_folder, args.analysis_name)
            if not os.path.exists(gatk_depth_of_coverage_file):
                gatk_depth_of_coverage_file = coverage_depth_stats()
            final_raw_vcf, final_raw_indel_vcf = varcall()
            final_raw_vcf = "%s/%s_aln_mpileup_raw.vcf_5bp_indel_removed.vcf" % (args.output_folder, args.analysis_name)
            filter(gatk_depth_of_coverage_file)
            annotation((final_raw_vcf).replace('_aln_mpileup_raw.vcf_5bp_indel_removed.vcf', '_filter2_final.vcf_no_proximate_snp.vcf'))
            tabix_vcf()
            stats()

## Sanity checks and directory structure maintenance methods
# Validate Filenames for any unsupported characters
def Validate_filename( name ):
    pattern_strings = ['\.', '\&', '\>', 'aaa', '\*']
    pattern_string = '|'.join(pattern_strings)
    searchobj = re.search(pattern_string, name, flags=0)
    if searchobj:
        print "- The file " + name + " contains unsupported characters such as quotes, spaces, or &:%?*><\$. \nPlease Provide another file name.\n"
        exit()

def file_exists(path1, path2, reference):
    if not os.path.isfile(path1):
        file_basename = os.path.basename(path1)
        keep_logging('- The input file {} does not exists. Please provide another file with full path or check the files path.\n'.format(file_basename), 'The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
        exit()
    if path2 is not None:
        if not os.path.isfile(path2):
            file_basename = os.path.basename(path2)
            keep_logging('- The input file {} does not exists. Please provide another file with full path or check the files path.\n'.format(file_basename), 'The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
            exit()
    if not os.path.isfile(reference):
        file_basename = os.path.basename(reference)
        keep_logging('- The reference fasta file {} does not exists. Please provide another with full path file with full path or check the files path.\n'.format(file_basename), 'The reference fasta file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
        exit()
    if ConfigSectionMap("pipeline", Config)['aligner'] == "bwa":
        ref_index_suffix1 = reference + ".bwt"
        ref_index_suffix2 = reference + ".amb"
        ref_index_suffix3 = reference + ".ann"
        ref_index_suffix4 = reference + ".sa"
        ref_index_suffix5 = reference + ".pac"
    elif ConfigSectionMap("pipeline", Config)['aligner'] == "bowtie":
        ref_index_suffix1 = reference + ".1.bt2"
        ref_index_suffix2 = reference + ".2.bt2"
        ref_index_suffix3 = reference + ".3.bt2"
        ref_index_suffix4 = reference + ".4.ebwt"
        ref_index_suffix5 = reference + ".rev.1.bt2"
        ref_index_suffix6 = reference + ".rev.2.bt2"
    if not os.path.isfile(ref_index_suffix1):
        keep_logging('- The reference index files given below does not exists:\n {}\n {}\n {}\n {}\n {}'.format(ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5), 'The reference index files given below does not exists:\n {}\n {}\n {}\n {}\n {}'.format(ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5), logger, 'warning')
        create_index(reference, ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5)
    else:
        keep_logging('- Index file already exists.', 'Index file already exists.', logger, 'info')

    ref_fai_index = reference + ".fai"
    if not os.path.isfile(ref_fai_index):
        keep_logging('- The reference fai index file {} required for samtools does not exists.'.format(ref_fai_index), 'The reference fai index file {} required for samtools does not exists.'.format(ref_fai_index), logger, 'warning')
        create_fai_index(reference, ref_fai_index)
    else:
        keep_logging('- Samtools fai Index file already exists.', 'Samtools fai Index file already exists.', logger, 'info')

    dict_name = os.path.splitext(os.path.basename(reference))[0] + ".dict"
    if not os.path.isfile(ConfigSectionMap(args.index, Config)['ref_path'] + "/" + dict_name):
        keep_logging('- The reference seq dict file {} required for GATK and PICARD does not exists.'.format(dict_name), 'The reference seq dict file {} required for GATK and PICARD does not exists.'.format(dict_name), logger, 'warning')
        picard_seqdict(dict_name, reference)
    else:
        keep_logging('- The reference seq dict file required for GATK and PICARD exists.', 'The reference seq dict file required for GATK and PICARD exists.', logger, 'info')

def java_check():
    jd = sp.check_output(["java", "-version"], stderr=sp.STDOUT)
    jd_version = jd.split('\n', 1)[0]
    if len(jd) < 1:
        keep_logging('- Unable to find a java runtime environment. The pipeline requires java 6 or later.', 'Unable to find a java runtime environment. The pipeline requires java 6 or later.', logger, 'exception')
    else:
        keep_logging('- Java Availability Check {}'.format(jd_version), 'Java Availability Check {}'.format(jd_version), logger, 'info')

def fileformat(file1, file2, final_out):
    print "- Checking File format....\n"
    if not file1.endswith('.fastq.gz'):
        base = os.path.basename(file1)
        os.path.splitext(base)
        file_1 = os.path.splitext(base)[0]
        cmdstring = "gzip -d " + file1 + " > " + final_out + file_1
        print "- Compressing input file " + base
        os.system(cmdstring)

    if not file2.endswith('.fastq.gz'):
        base = os.path.basename(file2)
        os.path.splitext(base)
        file_2 = os.path.splitext(base)[0]
        cmdstring = "gzip -d " + file2 + " > " + final_out + file_2
        print "Compressing input file " + base
        os.system(cmdstring)

def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('- Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()

def create_index(reference,ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5):
    aligner = ConfigSectionMap("pipeline", Config)['aligner']
    keep_logging('- Creating Index of reference fasta file for {} aligner.'.format(aligner), 'Creating Index of reference fasta file for {} aligner'.format(aligner), logger, 'info')
    if aligner == "bwa":
        cmd = "%s %s %s" % (ConfigSectionMap("bwa", Config)['base_cmd'], ConfigSectionMap("bwa", Config)['index'], reference)
        keep_logging(cmd, cmd, logger, 'debug')
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
                keep_logging('- Error in {} Indexer. Exiting.'.format(aligner), 'Error in {} Indexer. Exiting.'.format(aligner), logger, 'exception')
                sys.exit(1)
        if not os.path.isfile(ref_index_suffix1):
            keep_logging('- The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), 'The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), logger, 'exception')
    elif aligner == "bowtie":
        cmd = "%s %s %s" % ( ConfigSectionMap("bowtie", Config)['build_cmd'], reference, reference)
        keep_logging(cmd, cmd, logger, 'debug')
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
                keep_logging('- Error in {} Indexer. Exiting.'.format(aligner), 'Error in {} Indexer. Exiting.'.format(aligner), logger, 'exception')
                sys.exit(1)
        if not os.path.isfile(ref_index_suffix1):
            keep_logging('- The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), 'The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), logger, 'exception')
    else:
        print "- Different Aligner in config file"

def create_fai_index(reference, ref_fai_index):
    keep_logging('- Creating FAI Index using Samtools.', 'Creating FAI Index using Samtools.', logger, 'info')
    cmd = "%s %s %s" % (ConfigSectionMap("samtools", Config)['base_cmd'], ConfigSectionMap("samtools", Config)['faiindex'], reference)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('- Error in Samtools FAI Indexing step. Exiting.', 'Error in Samtools FAI Indexing step. Exiting.', logger, 'exception')
        sys.exit(1)


    if not os.path.isfile(ref_fai_index):
        keep_logging('- The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), 'The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), logger, 'exception')
    else:
        keep_logging('- Samtools Fai Index file created.', 'Samtools Fai Index file created.', logger, 'info')

def picard_seqdict(dict_name, reference):
    #dict_name = os.path.splitext(os.path.basename(reference_filename))[0] + ".dict"
    keep_logging('- Creating Sequence Dictionary using Picard.', 'Creating Sequence Dictionary using Picard.', logger, 'info')
    cmd = "%s CreateSequenceDictionary R=%s O=%s/%s" % (ConfigSectionMap("picard", Config)['base_cmd'], reference, ConfigSectionMap(args.index, Config)['ref_path'], dict_name)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('- Error in Picard Sequence Dictionary creation step. Exiting.', 'Error in Picard Sequence Dictionary creation step. Exiting.', logger, 'exception')
        sys.exit(1)

def cleanup(args, logger):
    if args.clean:
        keep_logging('- Removing Intermediate Temporary files', 'Removing Intermediate Temporary files', logger, 'info')
        os.chdir(args.output_folder)
        os.system("rm %s/*.fq.gz %s/*.sam %s/*_aln.bam %s/*_marked.bam %s/*_marked.bai %s/*_unmapped.bed %s/*_temp_paste_file.txt" % (args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder))
        make_sure_path_exists("%s/%s_stats_results" % (args.output_folder, args.analysis_name))
        os.system("mv %s/*depth_of_coverage* %s/*_stats %s/*_markduplicates_metrics %s/*_values.txt %s/*_INFO.txt %s/*_collect_alignment_metrics.txt %s/*_collect_wgs_metrics.txt  %s/*_gc_bias_metrics.txt %s/%s_stats_results" % (args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.analysis_name))
        make_sure_path_exists("%s/%s_logs" % (args.output_folder, args.analysis_name))
        os.system("mv %s/*.log.txt %s/%s_logs" % (args.output_folder, args.output_folder, args.analysis_name))
        make_sure_path_exists("%s/%s_vcf_results" % (args.output_folder, args.analysis_name))
        os.system("mv %s/*_unmapped.bed_positions %s/header.txt %s/*.vcf* %s/%s_vcf_results" % (args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.analysis_name))


def downsample(args, logger):
    method_start_time = datetime.now()
    if args.coverage_depth:
        keep_logging('- Downsampling Coverage Depth to: %s' % args.coverage_depth, 'Downsampling Coverage Depth to: %s' % args.coverage_depth, logger, 'info')

    # Run Mash to estimate Genome size
    mash_cmd = "mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout" % args.forward_raw

    keep_logging('- Running Mash Sketch',
                 'Running: %s' % mash_cmd, logger, 'info')

    if args.genome_size:
        gsize = int(args.genome_size)
        keep_logging('- Using Genome Size: %s' % gsize,
                     'Using Genome Size: %s' % gsize, logger, 'info')
    else:
        try:
            call(mash_cmd, logger)
        except sp.CalledProcessError:
            keep_logging('- Error running Mash for estimating genome size.', 'Error running Mash for estimating genome size', logger, 'exception')
            sys.exit(1)

        with open("/tmp/sketch_stdout", 'rU') as file_open:
            for line in file_open:
                if line.startswith('Estimated genome size:'):
                    gsize = float(line.split(': ')[1].strip())
                if line.startswith('Estimated coverage:'):
                    est_cov = float(line.split(': ')[1].strip())
        file_open.close()

        keep_logging('- Estimated Genome Size from Mash Sketch: %s' % gsize,
                     'Estimated Genome Size from Mash Sketch: %s' % gsize, logger, 'info')

    # Extract basic fastq reads stats with seqtk
    seqtk_check = "seqtk fqchk -q3 %s > /tmp/%s_fastqchk.txt" % (args.forward_raw, os.path.basename(args.forward_raw))

    keep_logging('- Running seqtk to extract Fastq statistics',
                 'Running seqtk to extract Fastq statistics: %s' % seqtk_check, logger, 'info')

    try:
        call(seqtk_check, logger)
    except sp.CalledProcessError:
        keep_logging('- Error running seqtk for extracting fastq statistics.', 'Error running seqtk for extracting fastq statistics.', logger, 'exception')
        sys.exit(1)

    with open("/tmp/%s_fastqchk.txt" % os.path.basename(args.forward_raw), 'rU') as file_open:
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

    keep_logging('- Average Read Length: %s' % avg_len,
                 'Average Read Length: %s' % avg_len, logger, 'info')

    keep_logging('- Total number of bases in fastq: %s' % total_bases,
                 'Total number of bases in fastq: %s' % total_bases, logger, 'info')

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)

    keep_logging('- Raw Coverage Depth: %s x' % ori_coverage_depth,
                 'Raw Coverage Depth: %s x' % ori_coverage_depth, logger, 'info')

    proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    (nproc, err) = proc.communicate()
    nproc = nproc.strip()

    if ori_coverage_depth > 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        # print factor
        r1_sub = "/tmp/%s" % os.path.basename(args.forward_raw)

        # Downsample using seqtk
        try:
            keep_logging("- Downsampling reads with seqtk",
                         "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                             args.forward_raw, factor, nproc, os.path.basename(args.forward_raw)), logger, 'info')
            call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                args.forward_raw, factor, nproc, os.path.basename(args.forward_raw)), logger)
        except sp.CalledProcessError:
            keep_logging('- Error running seqtk for downsampling raw fastq reads.',
                         'Error running seqtk for downsampling raw fastq reads.', logger, 'info')
            sys.exit(1)

        if args.reverse_raw:
            r2_sub = "/tmp/%s" % os.path.basename(args.reverse_raw)

            try:
                # keep_logging("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                #     args.reverse_raw, factor, nproc, os.path.basename(args.reverse_raw)),
                #              "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                #                  args.reverse_raw, factor, nproc, os.path.basename(args.reverse_raw)), logger, 'info')
                call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                    args.reverse_raw, factor, nproc, os.path.basename(args.reverse_raw)), logger)
            except sp.CalledProcessError:
                keep_logging('- Error running seqtk for downsampling raw fastq reads.',
                             'Error running seqtk for downsampling raw fastq reads.', logger, 'exception')
                sys.exit(1)
        else:
            r2_sub = "None"
    else:
        r1_sub = args.forward_raw
        if args.reverse_raw:
            r2_sub = args.reverse_raw
        else:
            r2_sub = "None"

    method_time_taken = datetime.now() - method_start_time
    keep_logging('- Downsampling completed in {}'.format(method_time_taken), 'Downsampling completed in {}'.format(method_time_taken), logger, 'info')
    return r1_sub, r2_sub

# Start of Main Method/Pipeline
if __name__ == '__main__':
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    args = parser().parse_args()
    global config_file
    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
    global logger
    if args.output_folder != '':
        args.output_folder += '/'
    make_sure_path_exists(args.output_folder)
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    logger = generate_logger(args.output_folder, args.analysis_name, log_unique_time)
    global Config
    global files_to_delete
    files_to_delete = []
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    pipeline(args, logger)
    cleanup(args, logger)
    keep_logging('- End: Pipeline', 'End: Pipeline', logger, 'info')
    time_taken = datetime.now() - start_time_2
    keep_logging('- Total Time taken for the pipeline: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')

