__author__ = 'alipirani'

import os
from modules.log_modules import keep_logging
from modules.logging_subprocess import *

def align_bowtie(base_cmd, forward_clean, reverse_clean, forward_unpaired, reverse_unpaired, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type, parameters):
    if type == "PE":
        cmd = "%s -x %s -1 %s -2 %s -S %s/%s_PE_aln.sam -t -p 8 %s %s" % (base_cmd, reference, forward_clean, reverse_clean, out_path, analysis, parameters, split_field)
        out_sam = "%s/%s_PE_aln.sam" % (out_path, analysis)
    else:
        cmd = "%s -x %s -U %s -S %s/%s_SE_aln.sam -t -p 8 %s %s" % (base_cmd, reference, forward_clean, out_path, analysis, parameters, split_field)
        out_sam = "%s/%s_SE_aln.sam" % (out_path, analysis)
    keep_logging("COMMAND: " + cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging(' - Error in Alignment step. Exiting.', 'Error in Alignment step. Exiting.', logger, 'exception')
        sys.exit(1)

    files_to_delete.append(out_sam)
    if not os.path.isfile(out_sam):
        keep_logging(' - Problem in BWA alignment. SAM file was not generated.', 'Problem in BWA alignment. SAM file was not generated', logger, 'exception')
        exit()
    else:
        return out_sam
