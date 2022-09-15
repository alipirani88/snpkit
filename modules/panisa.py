__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from modules.log_modules import keep_logging
from modules.logging_subprocess import *

def panisa(bam, out_path, analysis, logger, Config):
    cmd = "panISa.py %s -o %s/%s_panISa.txt" % (bam, out_path, analysis)
    keep_logging('Running Panisa: [%s]' % cmd, 'Running Panisa: [%s]' % cmd, logger, 'info')
    call(cmd, logger)
    panISa_output = "%s/%s_panISa.txt" % (out_path, analysis)
    return panISa_output

def ISFinder(panISa_output, out_path, analysis, logger, Config):
    cmd = "ISFinder_search.py %s -o %s/%s_ISFinder.txt" % (panISa_output, out_path, analysis)
    keep_logging('Running ISFinder: [%s]' % cmd, 'Running ISFinder: [%s]' % cmd, logger, 'info')
    call(cmd, logger)
    ISFinder_output = "%s/%s_ISFinder.txt" % (out_path, analysis)
    return ISFinder_output
