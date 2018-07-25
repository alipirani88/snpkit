# Add gubbins filtered masks to full alignment.
# Note: the full alignment should be the same one used as input to gubbins.
# Input:
#   (1) Whole genome alignment file in FASTA format
#   (2) Gubbins GFF output file
# Output:
#   (1) Masked whole genome alignment file in FASTA format
#   (2) Text file with list of masked positions in each genome
#   (3) Masked variant-only alignment file in FASTA format

# import libraries
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from collections import defaultdict
import re
import os

# parse command line arguments (input files)
parser = argparse.ArgumentParser(
                description='''Add gubbins filtered masks to full alignment.
                    Output:
                    (1) Masked whole genome alignment file in FASTA format
                    (2) Text file with list of masked positions in each genome
                    ''')
parser.add_argument('-a', '--alignment', metavar='FASTA', type=str,
                    required=True,
                    help='whole genome alignment file in FASTA format')
parser.add_argument('-g', '--gff', metavar='GFF', type=str,
                    required=True,
                    help='gubbins GFF output file')
parser.add_argument('-o', '--outdir', metavar='OUTFILE', type=str,
                    default='.',
                    help='''output fasta file directory
                    (default: .)''')
args = parser.parse_args()

aln_path = args.alignment
gubbins_gff_path = args.gff
outdir = args.outdir

# read in alignment and gubbins gff file
print('Reading ', gubbins_gff_path, '.', sep='')
gff = pd.read_csv(gubbins_gff_path, sep='\t', skiprows=2, header=None)
print('Reading ', aln_path, '.', sep='')
aln = AlignIO.read(aln_path, 'fasta')

# get indices/positions of recombinant regions identified by gubbins
print('Getting recombinant positions.')
recomb_regions = defaultdict(list)

for row in gff.iterrows():
    start = row[1][3]
    end = row[1][4]
    region = list(range(start, end))
    taxa = row[1][8].split(';')[2]
    taxa = taxa.replace('taxa=\"', '')
    taxa = taxa.replace('\"', '')
    taxa = list(taxa.split())
    for isolate in taxa:
        for position in region:
            recomb_regions[isolate].append(position)

# mask indices/positions of recombinant regions identified by gubbins
print('Masking recombinant positions in whole genome alignment.')
sample_masked_indices = defaultdict(list)
new_aln = list()

for record in aln:
    seq_str = list(str(record.seq))
    masked_indices = recomb_regions.get(record.id, [])
    for index in masked_indices:
        seq_str[index] = 'N'
    seq_str = ''.join(seq_str)
    new_record = SeqRecord(Seq(seq_str), id=record.id, description='')
    sample_masked_indices[record.id] = masked_indices
    new_aln.append(new_record)

# write new FASTA file with recombinant regions masked
fasta_outfile = outdir + '/' + re.split('/|\.', aln_path)[-2] + \
                '_gubbins_masked.fa'
text_outfile = outdir + '/' + re.split('/|\.', aln_path)[-2] + \
               '_masked_recomb_positions.txt'
var_site_outfile = outdir + '/' + re.split('/|\.', aln_path)[-2] + \
                '_gubbins_masked_var_sites.fa'

print('Writing', fasta_outfile)
with open(fasta_outfile, 'w') as handle:
    SeqIO.write(new_aln, handle, 'fasta')

# Write text file with list of recombinant sites for each genome
print('Writing', text_outfile)
with open(text_outfile, 'w') as handle:
    for sample, positions in sample_masked_indices.items():
        line = str(sample) + '\t' + ','.join(map(str, positions)) + '\n'
        handle.write(line)

# Get variant sites and write to fasta file using snp-sites
print('Getting variant sites using snp-sites.')
cmd = '/nfs/esnitkin/bin_group/anaconda3/bin/snp-sites ' + fasta_outfile + \
      ' -m ' + ' -o ' + var_site_outfile

os.system(cmd)
print(var_site_outfile, 'written.')

print('Done.')
