__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap

def markduplicates(out_sorted_bam, out_path, analysis, files_to_delete, logger, Config):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("picard", Config)['picard_bin'] + "/" + ConfigSectionMap("picard", Config)['base_cmd']
    keep_logging('Removing PCR duplicates using PICARD', 'Removing PCR duplicates using PICARD', logger, 'info')
    cmd = "java -jar %s MarkDuplicates REMOVE_DUPLICATES=true INPUT=%s OUTPUT=%s/%s_aln_marked.bam METRICS_FILE=%s/%s_markduplicates_metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES=500" % (base_cmd, out_sorted_bam, out_path, analysis, out_path, analysis)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
            keep_logging('Error in Picard Duplicates Removal step. Exiting.', 'Error in Picard Duplicates Removal step. Exiting.', logger, 'exception')
            sys.exit(1)
    out_marked_bam = "%s/%s_aln_marked.bam" % (out_path, analysis)
    #files_to_delete.append(out_marked_bam)
    if not os.path.isfile(out_marked_bam):
        keep_logging('Problem in Picard MarkDuplicate Step', 'Problem in Picard MarkDuplicate Step', logger, 'exception')
        exit()
    else:
        return out_marked_bam

def picard_seqdict(reference_filename, reference, logger, Config):
    dict_name = os.path.splitext(os.path.basename(reference_filename))[0] + ".dict"
    cmd = "java -jar %s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s/%s" % (base_cmd, reference_filename, ConfigSectionMap(reference, Config)['ref_path'],dict_name)
    os.system(cmd)



def picardstats(out_sorted_bam, out_path, analysis, reference, logger, Config):
    reference_filename = ConfigSectionMap(reference, Config)['ref_path'] + "/" + ConfigSectionMap(reference, Config)['ref_name']
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("picard", Config)[
        'picard_bin'] + "/" + ConfigSectionMap("picard", Config)['base_cmd']

    cmd = "java -jar %s CollectWgsMetrics I=%s O=%s/%s_collect_wgs_metrics.txt R=%s" % (
    base_cmd, out_sorted_bam, out_path, analysis, reference_filename)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Picard CollectWgsMetrics step. Exiting.',
                     'Error in Picard CollectWgsMetrics step. Exiting.', logger, 'exception')
        sys.exit(1)

    cmd = "java -jar %s CollectAlignmentSummaryMetrics I=%s O=%s/%s_collect_alignment_metrics.txt R=%s" % (
        base_cmd, out_sorted_bam, out_path, analysis, reference_filename)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Picard CollectWgsMetrics step. Exiting.',
                     'Error in Picard CollectWgsMetrics step. Exiting.', logger, 'exception')
        sys.exit(1)

    cmd = "java -jar %s CollectGcBiasMetrics I=%s O=%s/%s_gc_bias_metrics.txt R=%s S=%s/%s_summary_metrics.txt CHART=%s/%s_gc_bias_metrics.pdf " % (
        base_cmd, out_sorted_bam, out_path, analysis, reference_filename, out_path, analysis, out_path, analysis)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Picard CollectWgsMetrics step. Exiting.',
                     'Error in Picard CollectWgsMetrics step. Exiting.', logger, 'exception')
        sys.exit(1)
