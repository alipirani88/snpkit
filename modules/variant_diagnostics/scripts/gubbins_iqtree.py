# script to run gubbins and make tree
# written in python3

# modules needed: gubbins

# import modules
import argparse
import sys
import re
import os
# set path - commenting out zena's settings
#sys.path.insert(1, '/nfs/esnitkin/bin_group/pipeline/Github/beast/')
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fasta_functions import mask_positions

# parse command line arguments
parser = argparse.ArgumentParser(description='''Make tree using iqtree;
                                 optionally run gubbins first.''')
parser.add_argument('alignment', metavar='ALN',
                    help='Alignment file in FASTA format.')
parser.add_argument('-ng', '--nogubbins', dest='ng',
                    action='store_true',
                    help='Don\'t do recombination filtering using gubbins.')
parser.add_argument('-w', '--whole_genome', dest='w',
                    action='store_true',
                    help='Make tree using whole genome alignment.')
parser.add_argument('-v', '--variants_only', dest='v',
                    action='store_true',
                    help='Make tree using just variant sites.')
parser.add_argument('-m', '--model', metavar='MODEL', default='GTR+G', action='store', dest='m',
                    help='Nucleotide substitution model to use in iqtree.')
parser.add_argument('-o', '--outgroup', metavar='OG', default=None, nargs='*', action='store', dest='o', help='Outgroup sample(s) in alignment file.')

args = parser.parse_args()

if not args.w and not args.v:
    parser.error('At least one of -w and -v required.')


# get prefix for output files
pref = (os.path.basename(args.alignment)).replace('.fa', '')


# change working directory to where alignment file is
wd = os.path.dirname(args.alignment)
print (wd)
os.chdir(wd)

# if perform recombination filtering with gubbins
if args.ng:
    print('Not running gubbins.')
    # use unmasked alignment
    fasta_wga = args.alignment
else:
    # create new alignment if outgroups present
    if args.o is not None:
        no_og_fasta = re.sub('.fa*', '_no-outgroup.fa', args.alignment)
        subset_fasta(args.o, args.alignment, no_og_fasta)
        fasta = no_og_fasta
    else:
        fasta = args.alignment
    # run gubbins
    gub = 'run_gubbins.py --prefix ' + pref + ' --threads 12 ' + fasta
    print('Gubbins command: ' + gub)
    os.system(gub)
    # mask recombinant variants in whole genome alignment
    fasta_wga = mask_positions(args.alignment, pref + '.recombination_predictions.gff', mask_all=args.o)
    print ('Masked whole genome alignment - %s' % fasta_wga)

# if build tree with only variants
if args.v:
    fasta_vars = rm_invar_sites(fasta_wga)
    if not os.path.exists('iqtree_var'):
        # Create target Directory
        os.mkdir(wd + '/iqtree_var')

    os.chdir(wd + '/iqtree_var')
    # run iqtree
    iqtree_var = '/nfs/esnitkin/bin_group/anaconda3/bin/iqtree -s' + \
                 fasta_vars + '-nt AUTO -bb 1000 -m ' + args.m + '-pre ' + pref
    print('iqtree variant sites command: ' + iqtree_var)
    os.system(iqtree_var)

# if build tree with whole genome alignment
if args.w:
    # make iqtree_wga directory

    if not os.path.exists('iqtree_masked_wga'):
        # Create target Directory
        os.mkdir(wd + '/iqtree_masked_wga')

    os.chdir(wd + '/iqtree_masked_wga')
    # run iqtree
    iqtree_wga = 'iqtree -s ' + wd + '/' + os.path.basename(fasta_wga) + ' -nt AUTO -bb 1000 -m ' + args.m + ' -pre ' + wd + '/iqtree_masked_wga/' + pref + ' -t PARS -ninit 2'
    print('iqtree WGA command: ' + iqtree_wga)
    os.system(iqtree_wga)

