__author__ = 'alipirani'
# Back up file
import sys
import os
import argparse
from config_settings import ConfigSectionMap
from modules.check_subroutines import *
from modules.stages import *

################################################### Command Line Argument Parsing ########################################################################################################################
parser = argparse.ArgumentParser(description='Assembly pipeline for Illumina PE data')
parser.add_argument('-PE1', action='store', dest="forward_raw", help='Paired End file 1', required=True)
parser.add_argument('-PE2', action='store', dest="reverse_raw", help='Paired End file 2', required=True)
parser.add_argument('-o', action='store', dest="output_folder", help='Output Path ending with output directory name to save the results', required=True)
parser.add_argument('-analysis', action='store', dest="analysis_name", help='Unique analysis name to save the results', required=True)
parser.add_argument('-c', action='store', dest="croplength", help='Crop Length in case needed')
parser.add_argument('-f', action='store', dest="bam_input", help='Input Bam')
parser.add_argument('-index_name', action='store', dest="index", help='Reference Index Name', required=True)
parser.add_argument('-r', action='store', dest="reference", help='Reference Fasta File to build index in case index not built')
args = parser.parse_args()
############################################################################ End #########################################################################################################################

################################################### Check Subroutines: Arguments, Input files, Reference Index ###########################################################################################
reference = ConfigSectionMap(args.index)['ref_path'] + "/" + ConfigSectionMap(args.index)['ref_name']
out_path = args.output_folder + "/"
make_sure_path_exists(out_path)
file_exists(args.forward_raw, args.reverse_raw, reference)
#fileformat(args.forward_raw, args.reverse_raw, out_path)
java_check()
############################################################################ End #########################################################################################################################



#
# ################################################### Stages: Raw data Pre-processing using Trimmomatic ####################################################################################################
# trimmomatic(args.forward_raw, args.reverse_raw, out_path, args.croplength)
# ############################################################################ End #########################################################################################################################
#
#
# # ################################################### Stages: Alignment ####################################################################################################################################
# split_field = prepare_readgroup(args.forward_raw)
# out_sam = align(args.bam_input, args.output_folder, args.index, split_field, args.analysis_name)
# #out_sam = "~/Desktop/test/aln.sam"
# out_sorted_bam = prepare_bam(out_sam, out_path, args.analysis_name)
# #out_sorted_bam = "/Users/alipirani/Desktop/sample_19/aln_sort.bam"
# # ############################################################################ End #########################################################################################################################
#
#

# ################################################### Stages: Variant Calling ##############################################################################################################################

out_sorted_bam = out_path + "/" + args.analysis_name + "_aln_sort_bam"
caller = ConfigSectionMap("pipeline")['variant_caller']
#print caller
if caller == "samtoolswithpostalignbam":
    print "\n################## Variant Calling using Samtools and post-align bam input files ##################\n"
    out_finalbam = post_align_bam(out_sorted_bam, out_path, args.index, args.analysis_name)
    final_raw_vcf = variant_calling(out_finalbam, out_path, args.index, args.analysis_name)
    print "The final raw VCF file: %s" % final_raw_vcf
    print "\n############################ END: samtoolswithpostalignbam Variant Calling. #######################\n"
elif caller == "gatkhaplotypecaller":
    print "\n################## Variant Calling using GATK haplotyper and post-align bam input files ##################\n"
    out_finalbam = post_align_bam(out_sorted_bam, out_path, args.index, args.analysis_name)
    final_raw_vcf = variant_calling(out_finalbam, out_path, args.index, args.analysis_name)
    print "The final raw VCF file: %s" % final_raw_vcf
    print "\n############################# END: gatkhaplotypecaller Variant Calling. ###################################\n"
elif caller == "samtools":
    print "\n################## Variant Calling using Samtools without post-align bam input files.##################\n"
    final_raw_vcf = variant_calling(out_sorted_bam, out_path, args.index, args.analysis_name)
    print "The final raw VCF file: %s" % final_raw_vcf
    print "\n################################ END: samtools Variant Calling. ########################################\n"
else:
    print "\nPlease provide Variant Caller name in config file under the section [pipeline].\n Options for Variant caller:\n samtools\nsamtoolswithpostalignbam\ngatkhaplotypecaller\n"
    exit()
############################################################################ End #########################################################################################################################

################################################### Stages: Statistics ##############################################################################################################################
alignment_stats_file = alignment_stats(out_sorted_bam, out_path)
#final_vcf_stats = vcf_stats(final_raw_vcf, out_path)
############################################################################ End #########################################################################################################################

################################################### Stages: VCF Statistics ##############################################################################################################################
#final_raw_vcf = "/Users/alipirani/Desktop/testing_analysis_name///test_purpose_aln_mpileup_raw.vcf"
#print final_raw_vcf
vcf_stats_file = vcf_stats(final_raw_vcf, out_path)
############################################################################ End #########################################################################################################################


################################################### Stages: Remove Unwanted Intermediate files ##############################################################################################################################
#out_sam = "/Users/alipirani/Desktop/testing_analysis_name///test_purpose_aln.sam"
#remove_files(args.analysis_name, out_path, out_sam, out_sorted_bam)
############################################################################ End #########################################################################################################################









