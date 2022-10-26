__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap

def samtobam(out_sam, out_path, analysis, files_to_delete, logger, Config):
    base_cmd = ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s view -Sb %s > %s/%s_aln.bam" % (base_cmd, out_sam, out_path, analysis)
    # keep_logging('SAM to BAM Conversion', 'SAM to BAM Conversion', logger, 'info')
    # keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in SAM-to-BAM Conversion step. Exiting.', 'Error in SAM-to-BAM Conversion step. Exiting.', logger, 'exception')
        sys.exit(1)
    out_bam = "%s/%s_aln.bam" % (out_path, analysis)
    #files_to_delete.append(out_bam)
    if not os.path.isfile(out_bam):
        keep_logging('Error in SAM-to-BAM Conversion step. Exiting.', 'Error in SAM-to-BAM Conversion step. Exiting.', logger, 'exception')
        exit()
    else:
        return out_bam

def sort_bam(out_bam, out_path, analysis, logger, Config):
    base_cmd = ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s sort %s -m 500M -@ 0 -o %s/%s_aln_sort.bam -T %s/%s_aln_sort_temp" % (base_cmd, out_bam, out_path, analysis, out_path, analysis)
    # keep_logging('Sorting BAM file', 'Sorting BAM file', logger, 'info')
    # keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in BAM Sorting step. Exiting.', 'Error in BAM sorting step. Exiting.', logger, 'exception')
        sys.exit(1)
    sort_bam = "%s/%s_aln_sort.bam" % (out_path, analysis)
    if not os.path.isfile(sort_bam):
        keep_logging('Error in BAM Sorting step. Exiting.', 'Error in BAM Sorting step. Exiting.', logger, 'exception')
        exit()
    else:
        #os.remove(out_bam)
        return sort_bam

def index_bam(out_sort_bam, out_path, logger, Config):
    base_cmd = ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s index %s" % (base_cmd, out_sort_bam)
    # keep_logging(cmd, cmd, logger, 'info')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Samtools Indexing step. Exiting.', 'Error in Samtools Indexing step. Exiting.', logger, 'exception')
        sys.exit(1)

def samtoolswithpostalignbam(out_finalbam, out_path, reference_filename, analysis):
    mpileup_parameters = ConfigSectionMap("samtools")['mpileup_parameters']
    reference = ConfigSectionMap(reference_filename)['ref_path'] + "/" + ConfigSectionMap(reference_filename)['ref_name']
    cmd = "%s mpileup %s %s %s > %s/%s_aln_mpileup_postalign_raw.vcf" % (base_cmd, mpileup_parameters, reference, out_finalbam, out_path, analysis)
    os.system(cmd)
    final_raw_vcf =  "%s/%s_aln_mpileup_postalign_raw.vcf" % (out_path, analysis)
    return final_raw_vcf

def samtools(out_finalbam, out_path, reference_filename, analysis, logger, Config):
    base_cmd = ConfigSectionMap("samtools", Config)['base_cmd']
    mpileup_parameters = ConfigSectionMap("samtools", Config)['mpileup_parameters']
    reference = ConfigSectionMap(reference_filename, Config)['ref_path'] + "/" + ConfigSectionMap(reference_filename, Config)['ref_name']
    bcf_base_cmd = ConfigSectionMap("bcftools", Config)['base_cmd']
    cmd = "%s mpileup %s %s %s | %s call -O v -v -c -o %s/%s_aln_mpileup_raw.vcf" % (base_cmd, mpileup_parameters, reference, out_finalbam, bcf_base_cmd, out_path, analysis)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Samtools Variant Calling step. Exiting.', 'Error in Samtools Variant Calling step. Exiting.', logger, 'exception')
        sys.exit(1)
    final_raw_vcf =  "%s/%s_aln_mpileup_raw.vcf" % (out_path, analysis)
    return final_raw_vcf

def flagstat(out_sorted_bam, out_path, analysis, logger, Config):
    base_cmd = ConfigSectionMap("samtools", Config)['base_cmd']
    cmd = "%s flagstat %s > %s/%s_alignment_stats" % (base_cmd, out_sorted_bam, out_path, analysis)
    # keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('- Error in Samtools Alignment Stats step. Exiting.', 'Error in Samtools Alignment Stats step. Exiting.', logger, 'exception')
        sys.exit(1)
    alignment_stats_file = "%s/%s_alignment_stats" % (out_path, analysis)
    return alignment_stats_file

# Samtools latest version has a bug for rmdup command.
# def rmdup(out_sort_bam, out_path, analysis):
#     cmd = "%s rmdup %s %s/%s_rmdup.bam" % (base_cmd, out_sort_bam, out_path, analysis)
#     out_rmdup_bam = os.system(cmd)
#     return out_rmdup_bam




