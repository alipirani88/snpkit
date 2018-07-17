__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *

def tabix(file_array, preset, logger, Config):
    keep_logging('Tabix indexing vcf files', 'Tabix indexing vcf files', logger,
                 'info')
    for file in file_array:
        keep_logging('Tabix indexing vcf files: %s\n' % file, 'Tabix indexing vcf files: %s\n' % file, logger,
                     'info')
        # Make this parallel
        call("%s/%s/bgzip -c %s > %s%s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], file, file, ".gz"), logger)
        call("%s/%s/tabix -p %s -f %s.gz" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("vcftools", Config)['tabix_bin'], preset, file),logger)
