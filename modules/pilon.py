__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from modules.log_modules import keep_logging
from modules.logging_subprocess import *

def pilon(bam, reference, out_path, analysis, logger, Config):
    out_path = out_path + "/pilon"
    reference_path = ConfigSectionMap(reference, Config)['ref_path'] + "/" + ConfigSectionMap(reference, Config)['ref_name']
    os.system("mkdir %s" % out_path)
    cmd = "pilon --genome %s --bam %s --output %s --outdir %s --changes --vcf --tracks --variant --dumpreads --verbose 2>%s/pilon.log" % (reference_path, bam, analysis, out_path, out_path)
    keep_logging('Running Pilon: [%s]' % cmd, 'Running Pilon: [%s]' % cmd, logger, 'info')
    call(cmd, logger)
    keep_logging('- Finished Running Pilon: [%s]' % cmd, '- Finished Running Pilon: [%s]' % cmd, logger, 'info')