<<<<<<< HEAD
# PIPELINE TO RUN BEAST ON FLUX
=======
# PIPELINE TO RUN SCOTTI ON FLUX
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be

# Import modules
import sys
import os
import argparse
from subprocess import call
import re
<<<<<<< HEAD
from fasta_functions import *
from scotti_functions import *

# Get command line arguments
parser = argparse.ArgumentParser(description='Generate SCOTTI xmls and run in BEAST2.')
parser.add_argument('fasta', nargs='*', help='fasta file(s) of sequences to run SCOTTI on')
parser.add_argument('--subsets','-s', help='csv file where each row is the names of the subset of  sequences in the fasta file to run SCOTTI on')
parser.add_argument('--overwrite',"-ov", action='store_true', help="Overwrite output files.")
parser.set_defaults(overwrite=False)
parser.add_argument('--dates',"-d", help='Input csv file with sampling dates associated to each sample name (same name as fasta file). If not specified, looks for SCOTTI_dates.csv',default=os.getcwd()+"SCOTTI_dates.csv")
parser.add_argument('--hosts',"-ho", help='Input csv file with sampling hosts associated to each sample name (same name as fasta file). If not specified, looks for SCOTTI_hosts.csv',default=os.getcwd()+"SCOTTI_hosts.csv")
parser.add_argument('--hostTimes',"-ht", help='Input csv file with the earliest times in which hosts are infectable, and latest time when hosts are infective. Use same host names as those used for the sampling hosts file. If not specified, looks for SCOTTI_hostTimes.csv',default=os.getcwd()+"SCOTTI_hostTimes.csv")
parser.add_argument('--maxHosts',"-m", help='Maximum number of hosts (between sampled, non-sampled, and non-observed) allowed.', type=int, default=-1)
parser.add_argument('--unlimLife',"-u", action='store_true', help="non-sampled, non-observed hosts have unlimited life span (unlimited contribution in time to the outbreak). By default non-sampled non-observed hosts have limited duration. Usually asymptomaic patients have limited contributions, while to model environmental contamination it's better to use unlimited hosts. Unlimited life span can also decrease computational demand.")#, dest='limL', action='store_true'
parser.set_defaults(unlimLife=False)
parser.add_argument('--penalizeMigration',"-p", action='store_true', help="penalize lineages that are still in a deme after deme closure. By default, don't (original version).")#, dest='limL', action='store_true'
parser.set_defaults(penalizeMigration=False)
parser.add_argument('--numIter',"-N", help='Number of iterations of the MCMC. Default=100,000,000.', type=int, default=100000000)
parser.add_argument('--mutationModel',"-mu", help='String represnting the mutation model to be used. Possibilities are \"HKY\" or \"JC\". Further alternaties can be specified by modifying the output xml manually.', default="HKY")
parser.add_argument('--n_reps', '-n', metavar='N', type=int, help='number of beast runs for each xml', default=3)
parser.add_argument('--invar_sites', '-i', metavar='FILE', help='''
                    option 1: fasta file to  get number of invariant As, Cs, Gs, and Ts not in alignment
                    option 2: text file with number of invariant As, Cs, Gs, and Ts in alignment''')
parser.add_argument('--gubbins_gff', '-g', metavar='GFF', help='''gubbins output GFF file to mask recombinant sites in alignment before counting invariant sites''')
parser.add_argument('--scripts_dir', '-sd', metavar='DIR', help='directory where all of the scripts are housed', default='/nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/beast')
args = parser.parse_args()
# Rename variables
fasta = args.fasta
subsets = args.subsets
n_reps = args.n_reps
invar_sites = args.invar_sites
gubbins_gff = args.gubbins_gff
scripts_dir = args.scripts_dir
overwrite = args.overwrite
dates = args.dates
hosts = args.hosts
hostTimes = args.hostTimes
maxHosts = args.maxHosts
unlimLife = args.unlimLife
penalizeMigration = args.penalizeMigration
numIter = args.numIter
mutationModel = args.mutationModel

# Change names of sample ids and host names if needed
new_input_files = add_prefixes(hosts, dates, hostTimes, fasta)

if new_input_files is not None:
    hosts = new_input_files[0]
    dates = new_input_files[1]
    hostTimes = new_input_files[2]
    fasta = new_input_files[3]

# Make subsets of fasta file
if subsets is not None:
    print('Getting subsets of sequences from ' + fasta)
    fas = get_fasta_subsets(subsets,fasta)
else:
    fas = fasta

# Get invariant site counts
if invar_sites is not None:
    if not invar_sites.endswith(tuple(['fa','fasta'])):
        # Path to file with invariant site counts
        invar_counts_file = invar_sites
    else:
        # Count invariant sites in whole genome alignment
        invar_counts_file = count_invar_sites(invar_sites, gubbins_gff)
    with open(invar_counts_file) as f:
        invar_counts = f.readline().strip().split(' ')
    fixedAs = int(invar_counts[0])
    fixedCs = int(invar_counts[1])
    fixedGs = int(invar_counts[2])
    fixedTs = int(invar_counts[3])
else:
    fixedAs = 0
    fixedCs = 0
    fixedGs = 0
    fixedTs = 0

# For each subset, generate the xml
xmls = []
for fa in fas:
    print('Processing ' + fa + '.')
    output = fa.split('.')[0]
    xml = generate_scotti_xml(fa, dates, hosts, hostTimes, output, overwrite, maxHosts, unlimLife, penalizeMigration, numIter, mutationModel, fixedAs, fixedCs, fixedGs, fixedTs)
    print(xml)
    xmls.append(xml)

# Generate input file for beast script
with open('input_beast.txt', 'w') as f:
    for x in xmls:
        for i in range(n_reps):
            f.write(x + '\n')

print('input_beast.txt generated.')

# Generate PBS scripts and submit jobs to flux
print('Generating PBS script(s) and submitting jobs to flux.')
call(['perl', scripts_dir + '/beast_pbs_flux.pl', 'input_beast.txt'])

=======

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
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be
