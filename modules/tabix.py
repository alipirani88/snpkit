__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *

def tabix(file_array, preset, logger, Config):
    for file in file_array:
        # Make this parallel
        call("bgzip -c %s > %s%s" % (file, file, ".gz"), logger)
        call("tabix -p %s -f %s.gz" % (preset, file),logger)
