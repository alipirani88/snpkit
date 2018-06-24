__author__ = 'alipirani'

import os
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *
from modules.log_modules import *

def bioawk_make_reference_size(reference, logger, Config):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bioawk", Config)['bioawk_bin'] + "/" + ConfigSectionMap("bioawk", Config)['base_cmd']
    command = base_cmd + " -c fastx '{ print $name, length($seq) }' < %s > %s.size" % (reference, reference)    
    try:
        call(command, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Bioawk step. Exiting.', 'Error in Bioawk step. Exiting.', logger, 'exception')
        sys.exit(1)
    reference_size_file = "%s.size" % reference
    return reference_size_file
