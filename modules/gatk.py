from __future__ import division
__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from modules.samtools import *
from modules.picard import *
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from sys import platform as _platform



############################################################### GATK: Indel Realignment ####################################################################################################
def indel_realign(out_marked_sort_bam_rename, reference, out_path, analysis):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    print "\n################## GATK: Indel Realignment. ##################\n"
    #require fai index of reference
    #require seq dict of reference
    reference_filename = ConfigSectionMap(reference)['ref_path'] + "/" + ConfigSectionMap(reference)['ref_name']
    ref_fai_index(reference_filename)
    picard_seqdict(reference_filename, reference)
    cmd = "java -jar %s -T RealignerTargetCreator -R %s -o %s/%s_aln_sort_marked.bam.list -I %s " % (base_cmd, reference_filename, out_path, analysis, out_marked_sort_bam_rename)
    os.system(cmd)
    cmd = "java -jar %s -I %s -R %s -T IndelRealigner -targetIntervals %s/%s_aln_sort_marked.bam.list -o %s/%s_aln_realigned.bam" % (base_cmd, out_marked_sort_bam_rename, reference_filename, out_path, analysis, out_path, analysis)
    print "\nRunning:\n [%s] \n" % cmd
    os.system(cmd)
    out_indel_realigned = "%s/%s_aln_realigned.bam" % (out_path, analysis)
    print "\n################## END: Indel Realignment. ##################\n"
    return out_indel_realigned

############################################################### END #########################################################################################################################


############################################################### GATK: Variant Calling #######################################################################################################
def gatkhaplotypecaller(out_finalbam, out_path, reference, analysis):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    reference_filename = ConfigSectionMap(reference)['ref_path'] + "/" + ConfigSectionMap(reference)['ref_name']
    cmd = "java -jar %s %s -R %s -I %s -o %s/%s_aln_gatk_raw.vcf" % (base_cmd, ConfigSectionMap("gatk")['haplotype_parameters'], reference_filename, out_finalbam, out_path, analysis)
    print "Running Command: [%s]" % cmd
    os.system(cmd)
    final_raw_vcf =  "%s/%s_aln_gatk_raw.vcf" % (out_path, analysis)
    return final_raw_vcf
############################################################### END #########################################################################################################################



############################################################### GATK: Variant Filtering #######################################################################################################
def gatk_filter1(final_raw_vcf, out_path, analysis, reference):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    gatk_filter1_parameter_expression = ConfigSectionMap("gatk")['gatk_filter1_parameter_expression']
    gatk_filter1_command = "java -jar %s -T VariantFiltration -R %s -o %s/%s_filter1_gatk.vcf --variant %s --filterExpression \"%s\" --filterName PASS_filter1" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter1_parameter_expression)
    print "Running Command: [%s]" % gatk_filter1_command
    os.system(gatk_filter1_command)
    filter_flag_command = "grep '#\|PASS_filter1' %s/%s_filter1_gatk.vcf > %s/%s_filter1_final.vcf" % (out_path, analysis, out_path, analysis)
    os.system(filter_flag_command)
    gatk_filter1_final_vcf = "%s/%s_filter1_final.vcf" % (out_path, analysis)
    return gatk_filter1_final_vcf
############################################################### END #########################################################################################################################

############################################################### GATK: Variant Filtering #######################################################################################################
def gatk_filter2(final_raw_vcf, out_path, analysis, reference, logger, Config, Avg_dp):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']
    if ConfigSectionMap(filter_criteria, Config)['avg_depth'] == "yes":
        print "The average depth filter is on."
        low_Dp = float(Avg_dp) / 2
        high_Dp = float(Avg_dp) * 5
        DP_filter = "DP > %s && DP < %s" % (int(low_Dp), int(high_Dp))
    else:
        DP_filter = "DP > %s" % ConfigSectionMap(filter_criteria, Config)['dp']
    MQ_filter = "MQ > %s" % ConfigSectionMap(filter_criteria, Config)['mq']
    FQ_filter = "FQ < %s" % ConfigSectionMap(filter_criteria, Config)['fq']
    FQ_filter2 = "FQ < %s" % ConfigSectionMap(filter_criteria, Config)['fq2']
    QUAL_filter = "QUAL > %s" % ConfigSectionMap(filter_criteria, Config)['qual']

    gatk_filter2_parameter_expression = "%s && %s && %s && %s && %s" % (FQ_filter, MQ_filter, QUAL_filter, DP_filter, FQ_filter2)
    #gatk_filter2_parameter_expression = ConfigSectionMap("gatk", Config)['gatk_filter2_parameter_expression']
    if os.path.exists(final_raw_vcf):
        gatk_filter2_command = "java -jar %s -T VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s --filterExpression \"%s\" --filterName PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
    else:
        gatk_filter2_command = "java -jar %s -T VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s.gz --filterExpression \"%s\" --filterName PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)

    filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_gatk.vcf > %s/%s_filter2_final.vcf" % (out_path, analysis, out_path, analysis)
    keep_logging(gatk_filter2_command, gatk_filter2_command, logger, 'debug')
    keep_logging(filter_flag_command, filter_flag_command, logger, 'debug')
    try:
        call(gatk_filter2_command, logger)
        call(filter_flag_command, logger)
        #print ""
    except sp.CalledProcessError:
        keep_logging('Error in GATK filter step. Exiting.', 'Error in GATK filter step. Exiting.', logger, 'exception')
        sys.exit(1)
    gatk_filter2_final_vcf = "%s/%s_filter2_final.vcf" % (out_path, analysis)
    return gatk_filter2_final_vcf
############################################################### END #########################################################################################################################



############################################################### GATK: VCF2Fasta #######################################################################################################
def gatk_vcf2fasta_filter1(only_snp_filter1_vcf_file, out_path, analysis, reference):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    vcf2fasta_filter1_cmd = "java -jar %s -R %s -T FastaAlternateReferenceMaker -o %s_filter1.fasta --variant %s" % (base_cmd, reference, only_snp_filter1_vcf_file, only_snp_filter1_vcf_file)
    #print vcf2fasta
    print "\n\nRunning Command: [%s]\n\n" % vcf2fasta_filter1_cmd
    os.system(vcf2fasta_filter1_cmd)
    if _platform == "darwin":
        change_header_cmd = "sed -i '' 's/>.*/>%s/g' %s_filter1.fasta" % (analysis, only_snp_filter1_vcf_file)
        os.system(change_header_cmd)
    else:
        change_header_cmd = "sed -i 's/>.*/>%s/g' %s_filter1.fasta" % (analysis, only_snp_filter1_vcf_file)
        os.system(change_header_cmd)
    gatk_vcf2fasta_filter1_file = "%s_filter1.fasta" % (only_snp_filter1_vcf_file)
    return gatk_vcf2fasta_filter1_file
############################################################### END #########################################################################################################################

############################################################### GATK: VCF2Fasta #######################################################################################################
def gatk_vcf2fasta_filter2(only_snp_filter2_vcf_file, out_path, analysis, reference, logger, Config):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    vcf2fasta_filter2_cmd = "java -jar %s -R %s -T FastaAlternateReferenceMaker -o %s_filter2.fasta --variant %s" % (base_cmd, reference, only_snp_filter2_vcf_file, only_snp_filter2_vcf_file)
    #print vcf2fasta
    keep_logging(vcf2fasta_filter2_cmd, vcf2fasta_filter2_cmd, logger, 'debug')
    try:
        call(vcf2fasta_filter2_cmd, logger)
        #print ""
    except sp.CalledProcessError:
        keep_logging('Error in GATK vcf 2 fasta step. Exiting.', 'Error in GATK vcf 2 fasta step. Exiting.', logger, 'exception')
        sys.exit(1)
    if _platform == "darwin":
        change_header_cmd = "sed -i '' 's/>.*/>%s/g' %s_filter2.fasta" % (analysis, only_snp_filter2_vcf_file)
        try:
            call(change_header_cmd, logger)
        #print ""
        except sp.CalledProcessError:
            keep_logging('Error in sed header change step. Exiting.', 'Error in sed header change step. Exiting.', logger, 'exception')
            sys.exit(1)
    else:
        change_header_cmd = "sed -i 's/>.*/>%s/g' %s_filter2.fasta" % (analysis, only_snp_filter2_vcf_file)
        try:
            call(change_header_cmd, logger)
        #print ""
        except sp.CalledProcessError:
            keep_logging('Error in sed header change step. Exiting.', 'Error in sed header change step. Exiting.', logger, 'exception')
            sys.exit(1)
    gatk_vcf2fasta_filter2_file = "%s_filter2.fasta" % (only_snp_filter2_vcf_file)
    return gatk_vcf2fasta_filter2_file
############################################################### END #########################################################################################################################



############################################################### GATK: DepthOfCoverage #######################################################################################################
def gatk_DepthOfCoverage(out_sorted_bam, out_path, analysis_name, reference, logger, Config):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    cmd = "java -jar %s -T DepthOfCoverage -R %s -o %s/%s_depth_of_coverage -I %s --summaryCoverageThreshold 5 --summaryCoverageThreshold 10 --summaryCoverageThreshold 15 --summaryCoverageThreshold 20 --summaryCoverageThreshold 25 --minBaseQuality 15" % (base_cmd, reference, out_path, analysis_name, out_sorted_bam)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
        #print ""
    except sp.CalledProcessError:
        keep_logging('Error in GATK Depth of Coverage step. Exiting.', 'Error in GATK Depth of Coverage step. Exiting.', logger, 'exception')
        sys.exit(1)
    gatk_depth_of_coverage_file = "%s/%s_depth_of_coverage.sample_summary" % (out_path, analysis_name)
    keep_logging('GATK Depth of Coverage file: {}'.format(gatk_depth_of_coverage_file), 'GATK Depth of Coverage file: {}'.format(gatk_depth_of_coverage_file), logger, 'debug')
    return gatk_depth_of_coverage_file
############################################################### GATK: DepthOfCoverage #######################################################################################################




