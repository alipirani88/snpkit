from __future__ import division
__author__ = 'alipirani'


import sys
import os
import argparse
import errno
import ConfigParser
import glob
import re
import multiprocessing
import subprocess as subprocess
from datetime import datetime
from joblib import Parallel, delayed
from config_settings import ConfigSectionMap
from argparse import RawTextHelpFormatter


""" Command Line Argument Parsing """
def parser():
    parser = argparse.ArgumentParser(description='\nGenerate a README pieplien for every pipeline run.\n', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    # Input this file at the start of run. It should be submitted by the user or submitter
    required.add_argument('-readme_meta', action='store', dest="readme_meta", help='Input metadata file', required=True)
    required.add_argument('-out_dir', action='store', dest="out_dir", help='Output Directory name where this README will be moved.', required=True)
    return parser

# Set up logging modules and config file
start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
start_time_2 = datetime.now()

# Pass arguments to args object
args = parser().parse_args()

global config_file
global log_unique_time
global Config_readme
global logger

log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')


Config_readme = ConfigParser.ConfigParser()
Config_readme.read(args.readme_meta)


readme_file = args.out_dir + "/README.md"

print readme_file

if not os.path.isfile(readme_file):
    f=open(readme_file, 'w+')
    f.write("Request submitted by: %s" % ConfigSectionMap("Main", Config_readme)['submitter'])
    f.close()

else:
    print "README file already exists: Overwriting this file"
    f = open(readme_file, 'w+')
    f.write("Request submitted by: %s\n" % ConfigSectionMap("Main", Config_readme)['submitter'])
    f.write("Project Name: %s\n" % ConfigSectionMap("Main", Config_readme)['project_name'])
    f.write("Date when pipeline was run: %s\n" % ConfigSectionMap("Main", Config_readme)['date'])
    f.write("Piepline Version: %s\n" % ConfigSectionMap("Main", Config_readme)['version'])
    f.write("Comments: %s\n" % ConfigSectionMap("Description", Config_readme)['comments'])
    f.write("Parameters used/changed: %s\n" % ConfigSectionMap("Description", Config_readme)['parameters'])
    f.write("Filters used for this pipeline run: %s\n" % ConfigSectionMap("Description", Config_readme)['filters'])
    f.close()












