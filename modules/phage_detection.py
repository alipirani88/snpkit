__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
import json
import subprocess

def run_phaster(reference_genome, outdir, logger, Config):
    out_name = (os.path.basename(reference_genome)).split('.')
    if not os.path.isfile("%s/%s_phaster_post.json" % (outdir, str(out_name[0]))):
        if not os.path.isfile("%s/summary.txt" % outdir) or os.stat("%s/summary.txt") == 0:
            keep_logging('Running Phaster on %s' % reference_genome, 'Running Phaster on %s' % reference_genome, logger,
                     'info')
            phaster_post_cmd = "wget --post-file=\"%s\" \"http://phaster.ca/phaster_api\" -O %s/%s_phaster_post.json" % (reference_genome, outdir, str(out_name[0]))
            keep_logging('Running: %s' % phaster_post_cmd, 'Running: %s' % phaster_post_cmd, logger, 'debug')
            call(phaster_post_cmd, logger)
    else:
        keep_logging("Phaster Post Json file %s/%s_phaster_post.json exists" % (outdir, str(out_name[0])), "Phaster Post Json: %s/%s_phaster_post.json exists" % (outdir, str(out_name[0])), logger,
                     'info')
def parse_phaster(reference_genome, outdir, logger, Config):
    out_name = (os.path.basename(reference_genome)).split('.')
    if os.path.isfile("%s/%s_phaster_post.json" % (outdir, str(out_name[0]))) and os.stat("%s/%s_phaster_post.json" % (outdir, str(out_name[0]))).st_size != 0:
        with open('%s/%s' % (outdir, str(out_name[0]) + "_phaster_post.json")) as json_data:
            data = json.load(json_data)
            keep_logging("Phaster Post Json file exists... The status of Phaster job id %s is %s" % (data["job_id"], data["status"]), "Phaster Post Json file exists... The status of Phaster job id %s is %s" % (data["job_id"], data["status"]), logger,
                         'info')
            phaster_get_cmd = "wget \"http://phaster.ca/phaster_api?acc=%s\" -O %s/%s" % (
            data["job_id"], outdir, str(out_name[0]) + "_phaster_get.json")
            if not os.path.isfile("%s/%s_phaster_get.json" % (outdir, str(out_name[0]))):
                keep_logging('Running: %s' % phaster_get_cmd, 'Running: %s' % phaster_get_cmd, logger, 'debug')
                call(phaster_get_cmd, logger)
            else:
                call(phaster_get_cmd, logger)
                with open('%s/%s' % (outdir, str(out_name[0]) + "_phaster_get.json")) as json_get_data:
                    get_data = json.load(json_get_data)
                    phaster_get_zip_cmd = "wget \"http://%s\" -O %s/%s_phaster_get.zip" % (str(get_data["zip"]), outdir, str(out_name[0]))
                    phaster_unzip_cmd = "unzip -o %s/%s_phaster_get.zip" % (outdir, str(out_name[0]))
                json_get_data.close()
                keep_logging('Running: %s' % phaster_get_zip_cmd, 'Running: %s' % phaster_get_zip_cmd, logger, 'debug')
                keep_logging('Running: %s' % phaster_unzip_cmd, 'Running: %s' % phaster_get_cmd, logger, 'debug')
                call(phaster_get_zip_cmd, logger)
                call(phaster_unzip_cmd, logger)
    if os.path.isfile("%s/summary.txt" % outdir):
        keep_logging('Extracting Phage region information from %s/summary.txt' % outdir, 'Extracting Phage region information from %s/summary.txt' % outdir, logger, 'info')
        get_phage_regions = "sed -n -e '/REGION/,$p' %s/summary.txt | awk 'NR>2' | awk -F' ' '{print $5}'" % outdir
        proc = subprocess.Popen([get_phage_regions], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = (out.strip()).split('\n')
        #print out
        phage_positions = []
        for i in out:
            i_range = i.split('-')
            end_range = int(i_range[1]) + 1
            phage_positions.extend(list(range(int(i_range[0]), end_range)))

        f_open = open("%s/phage_region_positions.txt" % outdir, 'w+')
        for pos in phage_positions:
            f_open.write(str(pos) + "\n")
    else:
        keep_logging('Phaster output file %s/summary.txt not found' % outdir,
                     'Phaster output file %s/summary.txt not found' % outdir, logger, 'exception')
        exit()

    keep_logging('Number of phage Positions: %s' % len(phage_positions),
                 'Number of phage Positions: %s' % len(phage_positions),
                 logger, 'info')
    keep_logging('The Phage region positions in this file %s/phage_region_positions.txt will be filtered out' % outdir,
                 'The Phage region positions in this file %s/phage_region_positions.txt will be filtered out' % outdir, logger, 'info')
    return "%s/phage_region_positions.txt" % outdir