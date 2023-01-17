__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *

def tabix(file_array, preset, logger, Config):
    for file in file_array:
        # keep_logging("bgzip -c %s > %s.gz" % (file, file), "bgzip -c %s > %s.gz" % (file, file), logger, 'debug')
        # keep_logging("tabix -p %s -f %s.gz" % (preset, file), "tabix -p %s -f %s.gz" % (preset, file), logger, 'debug')
        call("bgzip -c %s > %s.gz" % (file, file), logger)
        call("tabix -p %s -f %s.gz" % (preset, file),logger)
