__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from sys import platform as _platform

def bamqc(out_sorted_bam, out_path, analysis, logger, Config):
    qualimap_base_command = ConfigSectionMap("qualimap", Config)['base_cmd']
    qualimap_bam_qc_cmd = "%s bamqc -bam %s -outdir %s -outfile %s_report.pdf -outformat pdf" % (qualimap_base_command, out_sorted_bam, out_path, analysis)
    keep_logging(qualimap_bam_qc_cmd, qualimap_bam_qc_cmd, logger, 'debug')
    try:
        call(qualimap_bam_qc_cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Qualimap step. Exiting.', 'Error in Qualimap step. Exiting.', logger, 'exception')
        sys.exit(1)
    qualimap_report_file = "%s/%s_report.pdf" % (out_path, analysis)
    keep_logging('Qualimap Report: {}'.format(qualimap_report_file), 'Qualimap Report: {}'.format(qualimap_report_file), logger, 'debug')
    return qualimap_report_file
