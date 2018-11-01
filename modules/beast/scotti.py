# PIPELINE TO RUN SCOTTI ON FLUX

# Import modules
import sys
import os
import argparse
from subprocess import call
import re

# Get command line arguments
parser = argparse.ArgumentParser(description='Generate SCOTTI XMLs and run SCOTTI.')
parser.add_argument('xmls', nargs='*', help='input SCOTTI BEAST xml file(s)')
parser.add_argument('--nReps', '-n', metavar='N', type=int, help='number of beast runs for each xml', default=3)
parser.add_argument('--invarSites', '-i', metavar='FILE', help='''
                    option 1: fasta file to  get number of invariant As, Cs, Gs, and Ts not in alignment
                    option 2: text file with number of invariant As, Cs, Gs, and Ts in alignment''')
parser.add_argument('--gubbinsGFF', '-g', metavar='GFF', help='''gubbins output GFF file to mask recombinant sites in alignment before counting invariant sites''')
parser.add_argument('--tree', '-t', help='starting tree file in newick format')
parser.add_argument('-scriptsDir', '-d', metavar='DIR', help='directory where all of the scripts are housed', default='/nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/beast')
args = parser.parse_args()
# Rename variables
xmls = args.xmls
nReps = args.nReps
invarSites = args.invarSites
gubbinsGFF = args.gubbinsGFF
tree = args.tree
scriptsDir = args.scriptsDir

# Set paths
python = '/nfs/esnitkin/bin_group/anaconda3/bin/python'
maskVarsPath = '/nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/mask_gubbins_variants.py'
