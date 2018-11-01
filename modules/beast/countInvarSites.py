# Count invariant sites in alignment file

# Import modules
from subprocess import call
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from collections import defaultdict
import re
import os

# Set paths
python = '/nfs/esnitkin/bin_group/anaconda3/bin/python'
maskVarsPath = '/nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/mask_gubbins_variants.py'


# Add gubbins filtered masks to full alignment.
# Note: the full alignment should be the same one used as input to gubbins.
# Input:
#    (1) Whole genome alignment file in FASTA format
#    (2) Gubbins GFF output file
#    (3) 
# Output:
#    (1) Masked whole genome alignment file in FASTA format
#    (2) Text file with list of masked positions in each genome
#    (3) Masked variant-only alignment file in FASTA format
def mask_positions(fasta, gff):



# Count invariant sites in an alignment file (fasta format)
# Input: 
#    (1) Path to alignment file in fasta format
#    (2) (optional) Path to [gubbins] GFF file to mask recombinant sites first
#    (3)
#    (4) (optional) Path to mask_gubbins_variants.py (default: '/nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/mask_gubbins_variants.py')    
# Output: 
#    (1) Text file (*_invar_site_counts.txt) with invariant site counts in the following order: A,C,G,T
#    (2) If GFF path given, masked fasta file (*_gubbins_masked.fa)
#    (3) Returns name of text file with invariant site counts (*_invar_site_counts.txt)
def count_invar_sites(fasta,gff=None,mask_vars_path=mask_vars_path)
    # Mask recombinant regions before counting invariant sites
    if gff is not None:
        call([python, mask_vars_path, fasta, gff])
        aln_file = re.split('/|\.', fasta)[-2] + \
                  '_gubbins_masked.fa'
    else:
        aln_file = fasta
    # Count invariant sites in whole genome alignment
    call([python, scripts_dir + '/get_num_invariant_sites_v2.py', aln_file])
    invar_counts_file = re.split('/|\.', alnFile)[-2] + \
                      '_invar_site_counts.txt'
    return(invar_counts_file)
   
