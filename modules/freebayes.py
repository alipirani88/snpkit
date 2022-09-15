__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap


def freebayes(out_finalbam, out_path, reference_filename, analysis, logger, Config):
    base_cmd = ConfigSectionMap("freebayes", Config)['base_cmd']
    freebayes_parameters = ConfigSectionMap("freebayes", Config)['freebayes_parameters']
    reference = ConfigSectionMap(reference_filename, Config)['ref_path'] + "/" + ConfigSectionMap(reference_filename, Config)['ref_name']
    cmd = "%s %s 4 %s -f %s %s > %s/%s_aln_freebayes_raw.vcf" % (base_cmd, reference.replace('.fasta', '.txt'), freebayes_parameters, reference, out_finalbam, out_path, analysis)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Freebayes Variant Calling step. Exiting.', 'Error in Freebayes Variant Calling step. Exiting.', logger, 'exception')
        sys.exit(1)
    final_raw_vcf_freebayes =  "%s/%s_aln_freebayes_raw.vcf" % (out_path, analysis)
    return final_raw_vcf_freebayes