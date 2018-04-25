__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import *
import json
import subprocess

def nucmer_repeat(reference, outdir, logger, Config):
    keep_logging('\nFinding repeat region in reference genome: %s\n' % reference, '\nFinding repeat region in reference genome: %s\n' % reference, logger,
                 'info')
    prefix = str(reference.split('.')[0]) + "_repeat"
    nucmer_repeat_cmd = "%s/%s/%s --maxmatch --nosimplify --prefix=%s %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("mummer", Config)['mummer_bin'], ConfigSectionMap("mummer", Config)['nucmer_base_cmd'], prefix, reference, reference)
    keep_logging('Running: %s' % nucmer_repeat_cmd, 'Running: %s' % nucmer_repeat_cmd, logger, 'debug')
    call(nucmer_repeat_cmd, logger)
    showcoords_cmd = "%s/%s/show-coords -I %s -r %s.delta > %s.coords" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("mummer", Config)['mummer_bin'], ConfigSectionMap("mummer", Config)['percent_id'], prefix, prefix)
    keep_logging('Running: %s' % showcoords_cmd, 'Running: %s' % showcoords_cmd, logger, 'debug')
    call(showcoords_cmd, logger)
    repeat_match_cmd = "%s/%s/repeat-match %s > %s.repeat_match" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("mummer", Config)['mummer_bin'], reference, prefix)
    tandem_repeats_cmd = "%s/%s/exact-tandems %s %s > %s_tandem_repeats_file" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("mummer", Config)['mummer_bin'], reference, ConfigSectionMap("mummer", Config)['min_tandem_repeat_length'], prefix)
    keep_logging('Running: %s' % tandem_repeats_cmd, 'Running: %s' % tandem_repeats_cmd, logger, 'debug')
    keep_logging('Running: %s' % repeat_match_cmd, 'Running: %s' % repeat_match_cmd, logger, 'debug')
    call(tandem_repeats_cmd, logger)
    call(repeat_match_cmd, logger)
    inexact_repeat_positions = []
    with open("%s.coords" % prefix) as fp:
        for i in xrange(6):
            fp.next()
        for line in fp:
            line = line.strip()
            line_split = line.split('|')
            range_str = str('-'.join(line_split[0].strip().split()))
            i_range = range_str.split('-')
            end_range = int(i_range[1]) + 1
            inexact_repeat_positions.extend(list(range(int(i_range[0]), end_range)))
            range_str = str('-'.join(line_split[1].strip().split()))
            i_range = range_str.split('-')
            end_range = int(i_range[1]) + 1
            inexact_repeat_positions.extend(list(range(int(i_range[0]), end_range)))
    fp.close()
    keep_logging('No. of inexact repeat matches positions: %s' % len(set(sorted(inexact_repeat_positions))), 'No. of inexact repeat matches: %s' % len(set(sorted(inexact_repeat_positions))), logger, 'info')
    tandem_repeats = []
    with open("%s_tandem_repeats_file" % prefix) as fp:
        for i in xrange(5):
            fp.next()
        for line in fp:
            line = line.strip()
            line_split = line.split()
            end_coords = int(line_split[0]) + int(line_split[1])
            tandem_repeats.extend(list(range(int(line_split[0]), end_coords)))
    keep_logging('No. of Tandem repeat matches positions: %s' % len(set(sorted(tandem_repeats))),
                 'No. of Tandem repeat matches positions: %s' % len(set(sorted(tandem_repeats))), logger, 'info')
    All_repeats = sorted(set(inexact_repeat_positions + tandem_repeats))
    f_open = open("%s/repeat_region_positions.txt" % outdir, 'w+')
    for pos in All_repeats:
        f_open.write(str(pos) + '\n')
    f_open.close()
    return "%s/repeat_region_positions.txt" % outdir