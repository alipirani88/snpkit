# script to run gubbins and make tree
# written in python3

# modules needed: gubbins

# import modules
import argparse
import sys
import re
import os

# set path
sys.path.insert(1, '/nfs/esnitkin/bin_group/pipeline/Github/beast/')

# import functions
import fasta_functions

# parse command line arguments
parser = argparse.ArgumentParser(description='''Make tree using iqtree;
                                 optionally run gubbins first.''')
parser.add_argument('alignment', metavar='ALN',
                    help='Alignment file in FASTA format.')
parser.add_argument('-ng', '--nogubbins',
                    help='Don\'t do recombination filtering using gubbins.', dest="ng")
parser.add_argument('-w', '--whole_genome',
                    action='store_true',
                    help='Make tree using whole genome alignment.', dest="w")
parser.add_argument('-v', '--variants_only',
                    action='store_true',
                    help='Make tree using just variant sites.', dest="v")
parser.add_argument('-m', '--model', metavar='MODEL', default='MFP',
                    help='Nucleotide substitution model to use in iqtree.', dest="m")
parser.add_argument('-o', '--outgroup', metavar='OG', default=None, nargs='*',
                    help='Outgroup sample(s) in alignment file.', dest="o")

args = parser.parse_args()



if not args.w and not args.v:
    parser.error('At least one of -w and -v required.')

# get prefix for output files
pref = args.alignment.split('/')[-1]
pref = pref.rsplit('.', 1)[-1]
print('prefix: ' + pref)

# change working directory to where alignment file is
wd = args.alignment.split('/')[0]
os.chdir(wd)

# modules to load
# modules = 'ml gubbins'
# os.system(modules)



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
    fasta_wga = mask_positions(args.alignment,
                               pref + '.recombination_predictions.gff',
                               mask_all=args.o)

# if build tree with only variants
if args.v:
    fasta_vars = rm_invar_sites(fasta_wga)
    # make iqtree_var directory
    os.mkdir('iqtree_var')
    os.chdir('iqtree_var')
    # run iqtree
    iqtree_var = '/nfs/esnitkin/bin_group/anaconda3/bin/iqtree -s' + \
                 fasta_vars + '-nt AUTO -bb 1000 -m ' + args.m + '-pre ' + pref
    print('iqtree variant sites command: ' + iqtree_var)
    os.system(iqtree_var)

# if build tree with whole genome alignment
if args.w:
    # make iqtree_wga directory
    #os.mkdir('iqtree_wga')
    os.system('mkdir iqtree_wga')
    os.chdir('iqtree_wga')
    # run iqtree
    iqtree_wga = '/nfs/esnitkin/bin_group/anaconda3/bin/iqtree -s' + \
                 fasta_wga + '-nt AUTO -bb 1000 -m ' + args.m + '-pre ' + pref
    print('iqtree WGA command: ' + iqtree_wga)
    os.system(iqtree_wga)

