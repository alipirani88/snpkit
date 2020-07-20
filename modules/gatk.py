from __future__ import division
__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from modules.picard import *
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from sys import platform as _platform

def gatk_filter(final_raw_vcf, out_path, analysis, reference, logger, Config, Avg_dp):
    if ConfigSectionMap("pipeline", Config)['variant_caller'] == "samtools":
        base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
        filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']
        print "Using variant filter parameters from: %s" % filter_criteria
        if ConfigSectionMap(filter_criteria, Config)['avg_depth'] == "yes":
            keep_logging('The average depth filter is turned on.', 'The average depth filter is turned on.', logger, 'info')
            low_Dp = float(Avg_dp) / 2
            high_Dp = float(Avg_dp) * 5
            DP_filter = "DP > %s && DP < %s" % (int(low_Dp), int(high_Dp))
        else:
            DP_filter = "DP > %s" % ConfigSectionMap(filter_criteria, Config)['dp']
        MQ_filter = "MQ > %s" % ConfigSectionMap(filter_criteria, Config)['mq']
        FQ_filter = "FQ < %s" % ConfigSectionMap(filter_criteria, Config)['fq']
        FQ_filter2 = "FQ < %s" % ConfigSectionMap(filter_criteria, Config)['fq2']
        QUAL_filter = "QUAL > %s" % ConfigSectionMap(filter_criteria, Config)['qual']
        AF_filter = "AF1 > %s" % float(ConfigSectionMap(filter_criteria, Config)['af'])
        gatk_filter2_parameter_expression = "%s && %s && %s && %s && %s && %s" % (FQ_filter, MQ_filter, QUAL_filter, DP_filter, FQ_filter2, AF_filter)
        if os.path.exists(final_raw_vcf):
            gatk_filter2_command = "%s VariantFiltration -R %s -O %s/%s_filter2_gatk.vcf --variant %s --filter-expression \"%s\" --filter-name PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
        else:
            gatk_filter2_command = "%s VariantFiltration -R %s -O %s/%s_filter2_gatk.vcf --variant %s.gz --filter-expression \"%s\" --filter-name PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
        filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_gatk.vcf > %s/%s_filter2_final.vcf" % (out_path, analysis, out_path, analysis)
        keep_logging(gatk_filter2_command, gatk_filter2_command, logger, 'debug')
        keep_logging(filter_flag_command, filter_flag_command, logger, 'debug')
        try:
            call(gatk_filter2_command, logger)
            call(filter_flag_command, logger)
        except sp.CalledProcessError:
            keep_logging('Error in GATK filter step. Exiting.', 'Error in GATK filter step. Exiting.', logger, 'exception')
            sys.exit(1)
        gatk_filter2_final_vcf = "%s/%s_filter2_final.vcf" % (out_path, analysis)
        return gatk_filter2_final_vcf
    elif ConfigSectionMap("pipeline", Config)['variant_caller'] == "gatkhaplotypecaller":
        base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
        filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']
        keep_logging("Using variant filter parameters from: %s" % filter_criteria, "Using variant filter parameters from: %s" % filter_criteria, logger, 'info')
        if ConfigSectionMap(filter_criteria, Config)['avg_depth'] == "yes":
            keep_logging('The average depth filter is turned on.', 'The average depth filter is turned on.', logger,
                         'info')
            low_Dp = float(Avg_dp) / 2
            high_Dp = float(Avg_dp) * 5
            DP_filter = "DP > %s && DP < %s" % (int(low_Dp), int(high_Dp))
        else:
            DP_filter = "DP > %s" % ConfigSectionMap(filter_criteria, Config)['dp']
        MQ_filter = "MQ > %s" % ConfigSectionMap(filter_criteria, Config)['mq']
        QUAL_filter = "QD > %s" % ConfigSectionMap(filter_criteria, Config)['qual']
        AF_filter = "AF > %s" % float(ConfigSectionMap(filter_criteria, Config)['af'])
        gatk_filter2_parameter_expression = "%s && %s && %s && %s" % (
        MQ_filter, QUAL_filter, DP_filter, AF_filter)
        if os.path.exists(final_raw_vcf):
            gatk_filter2_command = "%s VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s --filter-expression \"%s\" --filter-name PASS_filter2" % (
            base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
        else:
            gatk_filter2_command = "%s VariantFiltration -R %s -o %s/%s_filter2_gatk.vcf --variant %s.gz --filter-expression \"%s\" --filter-name PASS_filter2" % (
            base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
        filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_gatk.vcf > %s/%s_filter2_final.vcf" % (
        out_path, analysis, out_path, analysis)
        keep_logging(gatk_filter2_command, gatk_filter2_command, logger, 'debug')
        keep_logging(filter_flag_command, filter_flag_command, logger, 'debug')
        try:
            call(gatk_filter2_command, logger)
            call(filter_flag_command, logger)
        except sp.CalledProcessError:
            keep_logging('Error in GATK filter step. Exiting.', 'Error in GATK filter step. Exiting.', logger,
                         'exception')
            sys.exit(1)
        gatk_filter2_final_vcf = "%s/%s_filter2_final.vcf" % (out_path, analysis)
        return gatk_filter2_final_vcf

def gatk_filter_contamination(final_raw_vcf, out_path, analysis, reference, logger, Config, Avg_dp):
    if ConfigSectionMap("pipeline", Config)['variant_caller'] == "samtools":
        base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
        filter_criteria = "contamination_filters"
        if ConfigSectionMap(filter_criteria, Config)['avg_depth'] == "yes":
            keep_logging('The average depth filter is turned on.', 'The average depth filter is turned on.', logger, 'info')
            low_Dp = float(Avg_dp) / 2
            high_Dp = float(Avg_dp) * 5
            DP_filter = "DP > %s && DP < %s" % (int(low_Dp), int(high_Dp))
        else:
            DP_filter = "DP > %s" % float(ConfigSectionMap(filter_criteria, Config)['dp'])
        MQ_filter = "MQ > %s" % float(ConfigSectionMap(filter_criteria, Config)['mq'])
        FQ_filter = "FQ > %s" % float(ConfigSectionMap(filter_criteria, Config)['fq'])
        QUAL_filter = "QUAL > %s" % float(ConfigSectionMap(filter_criteria, Config)['qual'])
        AF_filter = "AF1 < %s" % float(ConfigSectionMap(filter_criteria, Config)['af'])
        gatk_filter2_parameter_expression = "%s && %s && %s && %s && %s" % (FQ_filter, MQ_filter, QUAL_filter, DP_filter, AF_filter)
        if os.path.exists(final_raw_vcf):
            gatk_filter2_command = "%s VariantFiltration -R %s -o %s/%s_filter2_contamination.vcf --variant %s --filter-expression \"%s\" --filter-name PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
        else:
            gatk_filter2_command = "%s VariantFiltration -R %s -o %s/%s_filter2_contamination.vcf --variant %s.gz --filter-expression \"%s\" --filter-name PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
        filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_contamination.vcf > %s/%s_filter2_final_contamination.vcf" % (out_path, analysis, out_path, analysis)
        keep_logging(gatk_filter2_command, gatk_filter2_command, logger, 'debug')
        keep_logging(filter_flag_command, filter_flag_command, logger, 'debug')
        try:
            call(gatk_filter2_command, logger)
            call(filter_flag_command, logger)
        except sp.CalledProcessError:
            keep_logging('Error in GATK filter step. Exiting.', 'Error in GATK filter step. Exiting.', logger, 'exception')
            sys.exit(1)
        gatk_filter2_final_contamination_vcf = "%s/%s_filter2_final_contamination.vcf" % (out_path, analysis)
        #extract_dp = "egrep -v \"^#\" %s | cut -f 8 | sed 's/^.*DP=\([0-9]*\);.*$/\1/' > %s/%s_depth_values.txt" % (gatk_filter2_final_contamination_vcf, out_path, analysis)
        extract_dp = "egrep -v \"^#\" %s | cut -f 8 | grep -Po 'DP=[0-9]*;?' | sed 's/DP=//g' | sed 's/;//g' > %s/%s_depth_values.txt" % (
        gatk_filter2_final_contamination_vcf, out_path, analysis)
        extract_pos = "grep -v '^#' %s | awk -F'\t' '{print $2}' > %s/%s_POS_values.txt" % (gatk_filter2_final_contamination_vcf, out_path, analysis)
        extract_fq = "awk -F'\t' '{print $8}' %s  | grep -o 'FQ=.*' | sed 's/FQ=//g' | awk -F';' '{print $1}' > %s/%s_FQ_values.txt" % (gatk_filter2_final_contamination_vcf, out_path, analysis)
        extract_mq = "egrep -v \"^#\" %s | cut -f 8 | sed 's/^.*MQ=\([0-9]*\);.*$/\1/' > %s/%s_MQ_values.txt" % (gatk_filter2_final_contamination_vcf, out_path, analysis)
        extract_af = "awk -F'\t' '{print $8}' %s  | grep -o 'AF1=.*' | sed 's/AF1=//g' | awk -F';' '{print $1}' > %s/%s_AF1_values.txt" % (gatk_filter2_final_contamination_vcf, out_path, analysis)
        try:
            call(extract_dp, logger)
            call(extract_pos, logger)
            call(extract_fq, logger)
            call(extract_mq, logger)
            call(extract_af, logger)
            keep_logging(extract_dp, filter_flag_command, logger, 'debug')
            keep_logging(extract_pos, filter_flag_command, logger, 'debug')
            keep_logging(extract_fq, filter_flag_command, logger, 'debug')
            keep_logging(extract_mq, filter_flag_command, logger, 'debug')
            keep_logging(extract_af, filter_flag_command, logger, 'debug')
        except sp.CalledProcessError:
            keep_logging('Error in GATK contamination filter step. Exiting.', 'Error in GATK contamination filter step. Exiting.', logger, 'exception')
            sys.exit(1)
        header = "pos,af"
        header_cmd = "echo \"%s\" > %s/header.txt" % (header, out_path)
        call(header_cmd, logger)
        paste_command = "paste -d, %s/%s_POS_values.txt %s/%s_AF1_values.txt > %s/%s_temp_paste_file.txt" % (out_path, analysis, out_path, analysis, out_path, analysis)
        call(paste_command, logger)
        combine_file_cmd = "cat %s/header.txt %s/%s_temp_paste_file.txt > %s/%s_INFO.txt" % (out_path, out_path, analysis, out_path, analysis)
        call(combine_file_cmd, logger)
        return gatk_filter2_final_contamination_vcf
    elif ConfigSectionMap("pipeline", Config)['variant_caller'] == "gatkhaplotypecaller":
        print "filter"

def gatk_filter_indel(final_raw_vcf, out_path, analysis, reference, logger, Config, Avg_dp):
    # if ConfigSectionMap("pipeline", Config)['variant_caller'] == "samtools":
    #     base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
    #     filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']
    #     if ConfigSectionMap(filter_criteria, Config)['avg_depth'] == "yes":
    #         keep_logging("Using variant filter parameters from: %s" % filter_criteria,
    #                      "Using variant filter parameters from: %s" % filter_criteria, logger, 'info')
    #         low_Dp = float(Avg_dp) / 2
    #         high_Dp = float(Avg_dp) * 5
    #         DP_filter = "DP > %s && DP < %s" % (int(low_Dp), int(high_Dp))
    #     else:
    #         DP_filter = "DP > %s" % ConfigSectionMap(filter_criteria, Config)['dp']
    #     MQ_filter = "MQ > %s" % ConfigSectionMap(filter_criteria, Config)['mq']
    #     FQ_filter = "FQ < %s" % ConfigSectionMap(filter_criteria, Config)['fq']
    #     FQ_filter2 = "FQ < %s" % ConfigSectionMap(filter_criteria, Config)['fq2']
    #     QUAL_filter = "QUAL > %s" % ConfigSectionMap(filter_criteria, Config)['qual']
    #     AF_filter = "AF1 > %s" % float(ConfigSectionMap(filter_criteria, Config)['af'])
    #
    #     gatk_filter2_parameter_expression = "%s && %s && %s && %s && %s && %s" % (FQ_filter, MQ_filter, QUAL_filter, DP_filter, FQ_filter2, AF_filter)
    #     if os.path.exists(final_raw_vcf):
    #         gatk_filter2_command = "java -jar %s -T VariantFiltration -R %s -o %s/%s_filter2_indel_gatk.vcf --variant %s --filterExpression \"%s\" --filterName PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
    #     else:
    #         gatk_filter2_command = "java -jar %s -T VariantFiltration -R %s -o %s/%s_filter2_indel_gatk.vcf --variant %s.gz --filterExpression \"%s\" --filterName PASS_filter2" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
    #
    #     filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_indel_gatk.vcf > %s/%s_filter2_indel_final.vcf" % (out_path, analysis, out_path, analysis)
    #     keep_logging(gatk_filter2_command, gatk_filter2_command, logger, 'debug')
    #     keep_logging(filter_flag_command, filter_flag_command, logger, 'debug')
    #     try:
    #         call(gatk_filter2_command, logger)
    #         call(filter_flag_command, logger)
    #     except sp.CalledProcessError:
    #         keep_logging('Error in GATK filter step. Exiting.', 'Error in GATK filter step. Exiting.', logger, 'exception')
    #         sys.exit(1)
    #     gatk_filter2_final_vcf = "%s/%s_filter2_indel_final.vcf" % (out_path, analysis)
    #     return gatk_filter2_final_vcf
    # elif ConfigSectionMap("pipeline", Config)['variant_caller'] == "gatkhaplotypecaller":
    base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
    filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']
    keep_logging("Using variant filter parameters from: %s" % filter_criteria, "Using variant filter parameters from: %s" % filter_criteria, logger, 'info')
    if ConfigSectionMap(filter_criteria, Config)['avg_depth'] == "yes":
        keep_logging('The average depth filter is turned on.', 'The average depth filter is turned on.', logger,
                     'info')
        low_Dp = float(Avg_dp) / 2
        high_Dp = float(Avg_dp) * 5
        DP_filter = "DP > %s && DP < %s" % (int(low_Dp), int(high_Dp))
    else:
        DP_filter = "DP > %s" % float(ConfigSectionMap(filter_criteria, Config)['dp'])
    MQ_filter = "MQ > %s" % float(ConfigSectionMap(filter_criteria, Config)['mq'])
    QUAL_filter = "QD > %s" % float(ConfigSectionMap(filter_criteria, Config)['qd'])
    AF_filter = "AF > %s" % float(ConfigSectionMap(filter_criteria, Config)['af'])

    #gatk_filter2_parameter_expression = "%s && %s && %s && %s" % (MQ_filter, QUAL_filter, DP_filter, AF_filter)
    gatk_filter2_parameter_expression = "%s && %s && %s && %s" % (MQ_filter, QUAL_filter, DP_filter, AF_filter)

    if os.path.exists(final_raw_vcf):
        gatk_filter2_command = "%s VariantFiltration -R %s -O %s/%s_filter2_indel_gatk.vcf --variant %s --filter-expression \"%s\" --filter-name PASS_filter2" % (
        base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)
    else:
        gatk_filter2_command = "%s VariantFiltration -R %s -O %s/%s_filter2_indel_gatk.vcf --variant %s.gz --filter-expression \"%s\" --filter-name PASS_filter2" % (
        base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter2_parameter_expression)

    filter_flag_command = "grep '#\|PASS_filter2' %s/%s_filter2_indel_gatk.vcf > %s/%s_filter2_indel_final.vcf" % (
    out_path, analysis, out_path, analysis)
    keep_logging(gatk_filter2_command, gatk_filter2_command, logger, 'debug')
    keep_logging(filter_flag_command, filter_flag_command, logger, 'debug')
    try:
        call(gatk_filter2_command, logger)
        call(filter_flag_command, logger)
    except sp.CalledProcessError:
        keep_logging('Error in GATK filter step. Exiting.', 'Error in GATK filter step. Exiting.', logger,
                     'exception')
        sys.exit(1)
    gatk_filter2_final_vcf = "%s/%s_filter2_indel_final.vcf" % (out_path, analysis)
    return gatk_filter2_final_vcf

def gatk_DepthOfCoverage(out_sorted_bam, out_path, analysis_name, reference, logger, Config):
    base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
    # cmd = "java -Xmx8G -jar %s/GenomeAnalysisTK.jar -T DepthOfCoverage -R %s -o %s/%s_depth_of_coverage -I %s --summaryCoverageThreshold 1 --summaryCoverageThreshold 5 --summaryCoverageThreshold 9 --summaryCoverageThreshold 10 --summaryCoverageThreshold 15 --summaryCoverageThreshold 20 --summaryCoverageThreshold 25 --ignoreDeletionSites --fix_misencoded_quality_scores" % (os.path.dirname(os.path.dirname(os.path.abspath(__file__))), reference, out_path, analysis_name, out_sorted_bam)

    cmd = "java -Xmx5G -jar %s/GenomeAnalysisTK.jar -T DepthOfCoverage -R %s -o %s/%s_depth_of_coverage -I %s --summaryCoverageThreshold 1 --summaryCoverageThreshold 5 --summaryCoverageThreshold 9 --summaryCoverageThreshold 10 --summaryCoverageThreshold 15 --summaryCoverageThreshold 20 --summaryCoverageThreshold 25 --ignoreDeletionSites" % (
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), reference, out_path, analysis_name, out_sorted_bam)

    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in GATK Depth of Coverage step. Exiting.', 'Error in GATK Depth of Coverage step. Exiting.', logger, 'exception')
        sys.exit(1)
    gatk_depth_of_coverage_file = "%s/%s_depth_of_coverage.sample_summary" % (out_path, analysis_name)
    keep_logging('GATK Depth of Coverage file: {}'.format(gatk_depth_of_coverage_file), 'GATK Depth of Coverage file: {}'.format(gatk_depth_of_coverage_file), logger, 'debug')
    return gatk_depth_of_coverage_file

def gatk_vcf2fasta_filter2(only_snp_filter2_vcf_file, out_path, analysis, reference, logger, Config):
    base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
    vcf2fasta_filter2_cmd = "%s FastaAlternateReferenceMaker -R %s -O %s_filter2.fasta --variant %s" % (base_cmd, reference, only_snp_filter2_vcf_file, only_snp_filter2_vcf_file)
    keep_logging(vcf2fasta_filter2_cmd, vcf2fasta_filter2_cmd, logger, 'debug')
    try:
        call(vcf2fasta_filter2_cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in GATK vcf 2 fasta step. Exiting.', 'Error in GATK vcf 2 fasta step. Exiting.', logger, 'exception')
        sys.exit(1)
    if _platform == "darwin":
        change_header_cmd = "sed -i '' 's/>.*/>%s/g' %s_filter2.fasta" % (analysis, only_snp_filter2_vcf_file)
        try:
            call(change_header_cmd, logger)
        except sp.CalledProcessError:
            keep_logging('Error in sed header change step. Exiting.', 'Error in sed header change step. Exiting.', logger, 'exception')
            sys.exit(1)
    else:
        change_header_cmd = "sed -i 's/>.*/>%s/g' %s_filter2.fasta" % (analysis, only_snp_filter2_vcf_file)
        try:
            call(change_header_cmd, logger)
        except sp.CalledProcessError:
            keep_logging('Error in sed header change step. Exiting.', 'Error in sed header change step. Exiting.', logger, 'exception')
            sys.exit(1)
    gatk_vcf2fasta_filter2_file = "%s_filter2.fasta" % (only_snp_filter2_vcf_file)
    return gatk_vcf2fasta_filter2_file

def gatkhaplotypecaller(out_finalbam, out_path, reference, analysis, logger, Config):
    base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
    reference_filename = ConfigSectionMap(reference, Config)['ref_path'] + "/" + ConfigSectionMap(reference, Config)['ref_name']
    cmd = "%s %s -R %s -I %s -O %s/%s_aln_mpileup_raw.vcf" % (base_cmd, ConfigSectionMap("gatk", Config)['haplotype_parameters'], reference_filename, out_finalbam, out_path, analysis)
    keep_logging('Running Command: [%s]' % cmd, 'Running Command: [%s]' % cmd, logger, 'info')
    #os.system(cmd)
    call(cmd, logger)
    final_raw_vcf =  "%s/%s_aln_mpileup_raw.vcf" % (out_path, analysis)
    return final_raw_vcf

""" Unused methods """
def gatk_filter1(final_raw_vcf, out_path, analysis, reference):
    base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
    gatk_filter1_parameter_expression = ConfigSectionMap("gatk")['gatk_filter1_parameter_expression']
    gatk_filter1_command = "%s VariantFiltration -R %s -O %s/%s_filter1_gatk.vcf --variant %s --filter-expression \"%s\" --filter-name PASS_filter1" % (base_cmd, reference, out_path, analysis, final_raw_vcf, gatk_filter1_parameter_expression)
    keep_logging('Running Command: [%s]' % gatk_filter1_command, 'Running Command: [%s]' % gatk_filter1_command, logger, 'info')
    os.system(gatk_filter1_command)
    filter_flag_command = "grep '#\|PASS_filter1' %s/%s_filter1_gatk.vcf > %s/%s_filter1_final.vcf" % (out_path, analysis, out_path, analysis)
    os.system(filter_flag_command)
    gatk_filter1_final_vcf = "%s/%s_filter1_final.vcf" % (out_path, analysis)
    return gatk_filter1_final_vcf

def indel_realign(out_marked_sort_bam_rename, reference, out_path, analysis):
    base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
    #require fai index of reference
    #require seq dict of reference
    reference_filename = ConfigSectionMap(reference)['ref_path'] + "/" + ConfigSectionMap(reference)['ref_name']
    ref_fai_index(reference_filename)
    picard_seqdict(reference_filename, reference)
    cmd = "%s RealignerTargetCreator -R %s -O %s/%s_aln_sort_marked.bam.list -I %s " % (base_cmd, reference_filename, out_path, analysis, out_marked_sort_bam_rename)
    os.system(cmd)
    cmd = "%s IndelRealigner -I %s -R %s -targetIntervals %s/%s_aln_sort_marked.bam.list O %s/%s_aln_realigned.bam" % (base_cmd, out_marked_sort_bam_rename, reference_filename, out_path, analysis, out_path, analysis)
    os.system(cmd)
    out_indel_realigned = "%s/%s_aln_realigned.bam" % (out_path, analysis)
    return out_indel_realigned

# def gatkhaplotypecaller(out_finalbam, out_path, reference, analysis):
#     base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("gatk", Config)['gatk_bin'] + "/" + ConfigSectionMap("gatk", Config)['base_cmd']
#     reference_filename = ConfigSectionMap(reference)['ref_path'] + "/" + ConfigSectionMap(reference)['ref_name']
#     cmd = "java -jar %s %s -R %s -I %s -o %s/%s_aln_gatk_raw.vcf" % (base_cmd, ConfigSectionMap("gatk")['haplotype_parameters'], reference_filename, out_finalbam, out_path, analysis)
#     keep_logging('Running Command: [%s]' % cmd, 'Running Command: [%s]' % cmd, logger, 'info')
#     os.system(cmd)
#     final_raw_vcf =  "%s/%s_aln_gatk_raw.vcf" % (out_path, analysis)
#     return final_raw_vcf

def gatk_vcf2fasta_filter1(only_snp_filter1_vcf_file, out_path, analysis, reference):
    base_cmd = ConfigSectionMap("gatk", Config)['base_cmd']
    vcf2fasta_filter1_cmd = "%s FastaAlternateReferenceMaker -R %s -o %s_filter1.fasta --variant %s" % (base_cmd, reference, only_snp_filter1_vcf_file, only_snp_filter1_vcf_file)
    keep_logging('Running Command: [%s]' % vcf2fasta_filter1_cmd, 'Running Command: [%s]' % vcf2fasta_filter1_cmd, logger, 'info')
    os.system(vcf2fasta_filter1_cmd)
    if _platform == "darwin":
        change_header_cmd = "sed -i '' 's/>.*/>%s/g' %s_filter1.fasta" % (analysis, only_snp_filter1_vcf_file)
        os.system(change_header_cmd)
    else:
        change_header_cmd = "sed -i 's/>.*/>%s/g' %s_filter1.fasta" % (analysis, only_snp_filter1_vcf_file)
        os.system(change_header_cmd)
    gatk_vcf2fasta_filter1_file = "%s_filter1.fasta" % (only_snp_filter1_vcf_file)
    return gatk_vcf2fasta_filter1_file








