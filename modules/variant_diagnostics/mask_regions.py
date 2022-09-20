__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
import json
import subprocess

def mask_regions(mask_file, outdir, logger, Config):
    keep_logging('- Masking custom region from bed file: %s' % mask_file,
                 '- Masking custom region from bed file: %s' % mask_file, logger,
                 'info')
    mask_positions_file = outdir + "mask_positions.txt"
    f1 = open(mask_positions_file, 'w+')
    with open(mask_file, 'rU') as fp:
        for line in fp:
            line_array = line.split('\t')
            lower_index = int(line_array[1])
            upper_index = int(line_array[2])
            for positions in range(lower_index, upper_index):
                f1.write(str(positions) + '\n')
    f1.close()
    return mask_positions_file