__author__ = 'alipirani'
import os
import sys
if sys.version_info < (3, 2):
    import subprocess32 as sp
else:
    import subprocess as sp
import errno
import gzip
import re
from modules.log_modules import keep_logging
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *

def trim(input1, input2, out_path, crop, logger, Config):
    if input2 != "None":
        adapter_file = ConfigSectionMap("Trimmomatic", Config)['adaptor_filepath']
        clean_filenames = out_path + ConfigSectionMap("Trimmomatic", Config)['f_p'] + " " + out_path + ConfigSectionMap("Trimmomatic", Config)['f_up'] + " " + out_path + ConfigSectionMap("Trimmomatic", Config)['r_p'] + " " + out_path + ConfigSectionMap("Trimmomatic", Config)['r_up']
        illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['seed_mismatches'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['simple_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] +  ConfigSectionMap("Trimmomatic", Config)['minadapterlength'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['keep_both_reads']
        sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic", Config)['window_size'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['window_size_quality']
        minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic", Config)['minlength']
        headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic", Config)['headcrop_length']
        if not crop:
            cmdstring = "trimmomatic PE -phred33 " + input1 + " " + input2 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string + " 2> %s/%s_trim_out.log" % (out_path, os.path.basename(os.path.dirname(out_path)))
            keep_logging(cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
            except sp.CalledProcessError:
                    keep_logging(' - Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                    sys.exit(1)
        else:
            crop_string = 'CROP:' + crop
            cmdstring = "trimmomatic PE " + input1 + " " + input2 + " " + clean_filenames + " " + crop_string + " " + illumina_string + " " + sliding_string + " " + minlen_string + " 2> %s/%s_trim_out.log" % (out_path, os.path.basename(os.path.dirname(out_path)))
            try:
                call(cmdstring, logger)
            except sp.CalledProcessError:
                keep_logging(' - Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                sys.exit(1)
    else:
        adapter_file = ConfigSectionMap("Trimmomatic", Config)['adaptor_filepath']
        clean_filenames = out_path + ConfigSectionMap("Trimmomatic", Config)['f_p']
        illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['seed_mismatches'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['simple_clipthreshold']
        sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic", Config)['window_size'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['window_size_quality']
        minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic", Config)['minlength']
        headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic", Config)['headcrop_length']
        if not crop:
            cmdstring = "trimmomatic SE " + input1 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string + " 2> %s/%s_trim_out.log" % (out_path, os.path.basename(os.path.dirname(out_path)))
            keep_logging(cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
            except sp.CalledProcessError:
                keep_logging(' - Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                sys.exit(1)

        else:
            crop_string = 'CROP:' + crop
            cmdstring = "trimmomatic SE " + input1 + " " + clean_filenames + " " + crop_string + " " + illumina_string + " " + sliding_string + " " + minlen_string + + " 2> %s/%s_trim_out.log" % (out_path, os.path.basename(os.path.dirname(out_path)))
            keep_logging(cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
            except sp.CalledProcessError:
                    keep_logging(' - Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                    sys.exit(1)