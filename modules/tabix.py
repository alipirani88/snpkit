__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *


def tabix(file_array, preset, logger, Config):
    for file in file_array:
        call("%s/%s/bgzip -c %s > %s%s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], file, file, ".gz"), logger)
        call("%s/%s/tabix -p %s -f %s.gz" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], preset, file),logger)
