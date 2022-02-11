# System wide imports
from __future__ import division
import sys
import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from collections import Counter
import errno
from pyfasta import Fasta
from datetime import datetime
import time
from subprocess import call
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Parse Command line Arguments
parser = argparse.ArgumentParser(
    description='Mask QC Filtered Variant positions in pre/post gubbins alignments - The script takes path to whole genome multi-sequence alignment, path to a SNP code matrix, list of Sample Ids to subset the matrix/alignment (Default is using all the Sample Ids in alignment and matrix)')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-alignment', action='store', dest="alignment",
                      help='Path to Multi-sequence Alignment to mask filtered positions')
required.add_argument('-matrix', action='store', dest="matrix",
                      help='Path to SNP code Matrix')
optional.add_argument('-sample_ids', action='store', dest="sample_ids",
                      help='List of Sample Ids to subset matrix and alignment. Default is All')
optional.add_argument('-reference_id', action='store', dest="reference_id",
                      help='Reference Genome Fasta ID to retain in masked alignment')
args = parser.parse_args()

def subset_fasta(ids, fastain, fastaout):
    """Creates a new fasta file with a subset of original sequences.

    Args:
        idnames: List of sequence names to be extracted.
        fastain: Path to fasta file containing all sequences of interest.
        fastaout: Name of output fasta file containing only the sequences of interest.

    Output:
        Fasta file containing only the sequences of interest.
    """

    with open(fastaout, 'w') as f:
        for rec in SeqIO.parse(fastain, 'fasta'):
            if str(rec.id) in ids:
                f.write('>' + ids[ids.index(rec.id)] + '\n')
                f.write(str(rec.seq) + '\n')

def get_fasta_subsets(csv, fasta):
    """Creates new fasta files with subsets of the original sequences.

    Args:
        csv: csv file where each line contains a different subset of sequence names to be extracted.
        fasta: Path to fasta file containing all sequences of interest.

    Output:
        Fasta files containing only the sequences of interest.

    Returns:
        List of paths to subsetted fasta files.
    """

    fas = []

    with open(csv) as f:
        csv_reader = csv.reader(f, delimiter=',')
        count = 0
        for subset in csv_reader:
            count += 1
            fastaout = fasta.split('.')[0] + '_subset' + str(count) + '.fasta'
            fa = subset_fasta(subset, fasta, fastaout)
            fas.append(fa)

    return(fas)

def mask_positions(fasta, gff, outdir='.', masked_sites_file=False, mask_all=None):
    """Masks positions in GFF file in alignment file in fasta format.

    If fasta is the whole genome alignment fasta used in gubbins
    and gff is the gubbins output GFF file, then recombinant
    regions identified by gubbins are masked (on a per-genome basis).
    Optionally returns text file with list of masked positions in each genome.

    Args:
        fasta: Alignment fasta file in which sites will be masked
               (ex. whole genome alignment that was input to gubbins).
        gff: GFF file containing sites that will be masked in certain genomes
             (ex. gubbins output GFF).
        outdir: Output file directory (default: current working directory).
        masked_sites_file: If true, generates a text file with a list of masked
                           positions in each genome.
        mask_all: list of taxa for which to mask positions that are recombinant
                  in any genome (ex. outgroup(s))

    Output:
        Masked whole genome alignment file in FASTA format.
        (Optional) Text file with list of masked positions in each genome.

    Returns:
        Path to masked fasta file.
    """

    # Read in alignment and gff file
    gff = pd.read_csv(gff, sep='\t', skiprows=2, header=None)
    aln = AlignIO.read(fasta, 'fasta')

    # Get indices/positions of recombinant regions identified by gubbins
    print 'Getting recombinant positions.'
    recomb_regions = defaultdict(list)
    all_recomb_regions = set()

    for row in gff.iterrows():
        start = row[1][3]
        end = row[1][4]
        region = list(range(start, end+1))
        taxa = row[1][8].split(';')[2]
        taxa = taxa.replace('taxa=\"', '')
        taxa = taxa.replace('\"', '')
        taxa = list(taxa.split())
        for isolate in taxa:
            for position in region:
                recomb_regions[isolate].append(position)
                all_recomb_regions.add(position)

    # Mask indices/positions of recombinant regions identified by gubbins
    print 'Masking recombinant positions in whole genome alignment.'
    sample_masked_indices = defaultdict(list)
    new_aln = list()

    for record in aln:
        seq_str = list(str(record.seq))
        if mask_all is None:
            masked_indices = recomb_regions.get(record.id, [])
        elif record.id in mask_all:  # mask all recombinant positions
            print 'Masking all positions in ' + record.id
            masked_indices = all_recomb_regions  # mask only recomb pos in tax
        else:
            masked_indices = recomb_regions.get(record.id, [])
        for index in masked_indices:
            seq_str[index-1] = 'N'
        seq_str = ''.join(seq_str)
        new_record = SeqRecord(Seq(seq_str), id=record.id, description='')
        sample_masked_indices[record.id] = masked_indices
        new_aln.append(new_record)

    # Write new FASTA file with recombinant regions masked
    fasta_outfile = outdir + '/' + re.split('/|\.', fasta)[-2] + \
        '_gubbins_masked.fa'

    print 'Writing %s' % fasta_outfile
    with open(fasta_outfile, 'w') as handle:
        SeqIO.write(new_aln, handle, 'fasta')
    handle.close()

    if masked_sites_file:
        # Write text file with list of recombinant sites for each genome
        text_outfile = outdir + '/' + re.split('/|\.', fasta)[-2] + \
                       '_masked_recomb_positions.txt'
        print 'Writing %s' % text_outfile
        with open(text_outfile, 'w') as handle:
            for sample, positions in sample_masked_indices.items():
                line = str(sample) + '\t' + ','.join(map(str, positions)) + \
                       '\n'
                handle.write(line)

    return(fasta_outfile)

def subset_SNP_matrix(matrix, sample_ids):
    print "Subset SNP code matrix - %s based on Sample Ids - %s" % (matrix, sample_ids)
    data = pd.read_csv(matrix, sep='\t')
    rowname = "Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product"
    id_list_array = []
    id_list_array.append(rowname)
    count = 0
    with open(sample_ids) as fp:
        for line in fp:
            line = line.strip()
            count = count + 1
            id_list_array.append(line)
    # print data[id_list_array]
    print "Number of Sample Ids - %s" % count
    data[id_list_array].to_csv('SNP_code_matrix_%s.tsv' % sample_ids.replace('.tsv', ''), index=False, sep='\t')
    subset_matrices = "SNP_code_matrix_%s.tsv" % sample_ids.replace('.tsv', '')
    return subset_matrices, id_list_array

def extract_filtered_positions(subset_matrices):
    filtered_positions = []
    filter_codes = ['4', '2', '-1', '-3', '-4']
    with open("%s" % subset_matrices, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        next(csv_reader, None)
        for row in csv_reader:
            if "2" in row[1:] or "4" in row[1:] or "-1" in row[1:] or "-3" in row[1:] or "-4" in row[1:]:
                filtered_positions.append(row[0].split(' ')[3])
    print "Number of QC filtered variant positions - %s" % len(filtered_positions)
    with open("%s_filtered_positions.txt" % subset_matrices.replace('.tsv', ''), 'w') as out:
        for i in filtered_positions:
            out.write(i + '\n')
    out.close()
    return filtered_positions

def mask_filtered_variant_positions(fasta, list_of_positions):
    """Masks positions in list_of_positions array in alignment file in fasta format.

    Args:
        fasta: Alignment fasta file in which sites will be masked
               (ex. whole genome alignment that is input to gubbins).
        list_of_positions: list of positions that will be masked in all genomes


    Output:
        Masked whole genome alignment file in FASTA format.
        (Optional) Text file with list of masked positions in each genome.

    Returns:
        Path to masked fasta file.
    """

    # Read in alignment
    aln = AlignIO.read(fasta, 'fasta')


    # Mask indices/positions of recombinant regions identified by gubbins
    print 'Masking QC filtered variant positions in whole genome alignment.'
    new_aln = list()

    for record in aln:
        seq_str = list(str(record.seq))
        for index in list_of_positions:
            #print len(seq_str)
            #print (int(index)-1)
            # seq_str[0] = 'N'
            seq_str[int(index)-1] = 'N'
        seq_str = ''.join(seq_str)
        new_record = SeqRecord(Seq(seq_str), id=record.id, description='')
        new_aln.append(new_record)

    fasta_outfile = fasta.replace('_temp.fa', '_QC_filtered_positions_masked.fa')

    print 'Writing %s' % fasta_outfile
    with open(fasta_outfile, 'w') as handle:
        SeqIO.write(new_aln, handle, 'fasta')
    handle.close()
    return(fasta_outfile)

if __name__ == '__main__':
    if args.sample_ids:
        subset_matrices, id_list_array = subset_SNP_matrix(args.matrix, args.sample_ids)
        print "Using subset matrices %s for sample IDs - %s" % (subset_matrices, args.sample_ids)
    else:
        subset_matrices = args.matrix
        print "Using subset matrices %s for sample IDs - All" % (subset_matrices)

    filtered_positions = extract_filtered_positions(subset_matrices)

    id_list_array = [sub.replace('R1_001.fastq.gz', '') for sub in id_list_array]
    id_list_array.append(args.reference_id)
    temp_subset_fasta = ((args.alignment).replace('.fa', '_%s' % args.sample_ids)).replace('.tsv', '_temp.fa')
    #print id_list_array
    #print filtered_positions
    #print temp_subset_fasta
    subset_fasta(id_list_array, args.alignment, temp_subset_fasta)
    fasta_outfile = mask_filtered_variant_positions(temp_subset_fasta, filtered_positions)

    gubbins_iqtree_script = "python %s/gubbins_iqtree.py -w %s" % (
                os.path.dirname(os.path.abspath(__file__)), fasta_outfile)

    job_file_name = fasta_outfile.replace('.fa', '.sh')

    with open(job_file_name, 'w') as out:
        out.write(gubbins_iqtree_script + '\n')
    out.close()

    # Chesk Masking
    # aln = AlignIO.read(fasta_outfile, 'fasta')
    # for record in aln:
    #     seq_str = list(str(record.seq))
    #     print seq_str[1]
