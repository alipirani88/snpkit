__author__ = 'alipirani'

import sys
import os
import argparse
import errno
from datetime import datetime
import ConfigParser
from config_settings import ConfigSectionMap
if sys.version_info < (3, 2):
    import subprocess32 as sp
else:
    import subprocess as sp


def submit_job
