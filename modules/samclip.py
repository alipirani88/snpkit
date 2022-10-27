__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap

def samclip(out_sam, out_path, analysis, reference, logger, Config):
    cmd = "samclip --ref %s --max 10 < %s > %s/%s_clipped.sam" % (reference, out_sam, out_path, analysis)
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in clipping soft and hard clipped alignments using samclip. Exiting.', 'Error in clipping soft and hard clipped alignments using samclip. Exiting.',
                     logger, 'exception')
        exit()
    out_clipped_sam = "%s/%s_clipped.sam" % (out_path, analysis)
    return out_clipped_sam