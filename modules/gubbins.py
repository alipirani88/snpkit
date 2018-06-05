__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *

def gubbins(gubbins_dir, input_fasta, logger, Config):
    keep_logging('\nRunning Gubbins on input: %s\n' % input_fasta, '\nRunning Gubbins on input: %s\n' % input_fasta, logger,
                 'info')
    call("cd %s" % ConfigSectionMap("gubbins", Config)['gubbins_bin'], logger)
    gubbins_cmd = "%s/%s --prefix %s/%s %s" % (ConfigSectionMap("gubbins", Config)['gubbins_bin'], ConfigSectionMap("gubbins", Config)['base_cmd'], gubbins_dir, (os.path.basename(input_fasta)).replace('.fa', ''), input_fasta)
    #call(gubbins_cmd, logger)
    keep_logging('\nRunning Gubbins: %s' % input_fasta, '\nRunning Gubbins: %s\n' % input_fasta,
                 logger,
                 'info')