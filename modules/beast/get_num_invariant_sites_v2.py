# Get number of invariant A's, C's, G's, and T's from a whole-genome alignment.
# These numbers can be used in BEAST analyses.
# V2: uses snp-sites to get VCF to find variant positions instead of SNP matrix.
# Input:
#   (1) Whole genome alignment file in FASTA format
#       (Note: This should be the one used to create the variant position
#       FASTA file used in BEAST.)
# Output:
#   (1) VCF file of variants (created by snp-sites)
#   (2) Text file with number of invariant A's, C's, G's, and T's in alignment (order: A C G T)

# import libraries
import argparse
import os
import re
from Bio import AlignIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
# import os
# from Bio import SeqIO
from collections import Counter
import pandas as pd

# parse command line arguments (input files)
parser = argparse.ArgumentParser(
                description='''Get number of invariant A's, C's, G's, and T's
                    from a whole-genome alignment.
                    Input:
                      (1) Whole genome alignment file in FASTA format
                          (The one used to create the variant position
                          FASTA file used in BEAST, or the recombination-masked version..)
                    Output:
		      (1) VCF file of variants (created by snp-sites)
                      (2) Text file with number of invariant
                          A's, C's, G's, and T's in alignment (order: A C G T)
                    ''')
parser.add_argument('fasta', metavar='FASTA', type=str,
                    help='whole genome alignment file in FASTA format')
parser.add_argument('-o', '--outdir', metavar='OUTDIR', type=str,
                    default='.',
                    help='''output directory
                    (default: .)''')
args = parser.parse_args()

aln_path = args.fasta
outdir = args.outdir

# read in alignment to be able to get invariant base counts
print('Reading ', aln_path, '.', sep='')
aln = AlignIO.read(aln_path, 'fasta')

# get variant positions
var_site_outfile = outdir + '/' + re.split('/|\.', aln_path)[-2] + \
                '_snp-sites.vcf'
print('Getting variant positions using snp-sites.')
cmd = '/nfs/esnitkin/bin_group/anaconda3/bin/snp-sites ' + aln_path + \
      ' -v ' + ' -o ' + var_site_outfile

os.system(cmd)

positions = []
with open(var_site_outfile) as f:
    for line in f:
        li=line.strip()
        if not li.startswith("#"):
            positions.append(line.split('\t')[1])

# get ref allele for invariant sites
for record in aln:
    seq_str = list(str(record.seq))
    for index in positions:
        index = int(index)
        seq_str[index] = ''
    invar_sites = ''.join(seq_str)
    break

# get invariant site count for each base
print('Counting bases.')
invar_counts = Counter(invar_sites)

# write base counts to files
invar_base_counts_file = outdir + '/' + re.split('/|\.', aln_path)[-2] + \
                '_invar_site_counts.txt'
print('Writing ', invar_base_counts_file, ' (order: A C G T).', sep='')
with open(invar_base_counts_file, 'w') as f:
    for base, count in sorted(invar_counts.items()):
        f.write('%s ' % (count))

print('Done.')
