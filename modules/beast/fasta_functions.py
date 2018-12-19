# Functions to manipulate fasta files

# Functions:
#    subset_fasta
#    mask_positions
#    count_invar_sites

# Import modules
from subprocess import call
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from collections import defaultdict
from collections import Counter
import re
import os

<<<<<<< HEAD
def subset_fasta(ids, fastain, fastaout):
    """Creates a new fasta file with a subset of original sequences.

    Args:
        idnames: List of sequence names to be extracted.
=======
def subset_fasta(idnames, fastain, fastaout):
    """Creates a new fasta file with a subset of original sequences.

    Args:
        idnames: Text file with list of names of sequences to be extracted.
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be
        fastain: Path to fasta file containing all sequences of interest.
        fastaout: Name of output fasta file containing only the sequences of interest.

    Output:
        Fasta file containing only the sequences of interest.
<<<<<<< HEAD
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
=======

    Returns:
        Path to subsetted fasta file.
    """    

    ids = []

    for line in open(idnames, 'r'):
        if line.strip() == 'x':
            pass
        else:
            ids.append(line.strip())

    with open(fastaout, 'w') as f:
        for rec in SeqIO.parse(fastain, 'fasta'):
            if str(rec.id) in ids:
                f.write('>' + names[ids.index(rec.id)] + '\n')
                f.write(str(rec.seq) + '\n')

    return(fastaout)
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be


def rm_invar_sites(fasta, outfile=None, outfmt=['fasta','vcf'], outdir='.', 
                   path={'snpsites':'/nfs/esnitkin/bin_group/anaconda3/bin/snp-sites'}):
    """Removes invariant sites from a multifasta file using snp-sites.

    Args:
        fasta: Alignment fasta file from which to get variant sites.
        outfile: (Optional) Output file name.
        outfmt: Format of output file (fasta or vcf).
        outdir: (Optional) Output directory (default: current working directory).
        path: Dictionary of paths including the path to snp-sites (key: snpsites, value: path)

    Output:
        Fasta or VCF file of variant sites.

    Returns:
        Path to output file.
    """

<<<<<<< HEAD
    if len(outfmt) > 1:
=======
   if len(outfmt) > 1:
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be
       outfmt = outfmt[0]

    # Get output fasta file name
    if outfile == None:
        outfile = outdir + '/' + re.split('/|\.', fasta)[-2] + 'var_sites.' + outfmt

    # Get output format
    if outfmt == 'fasta':
        flag = ' -m '
    elif outfmt == 'vcf':
        flag = ' -v '

    # Get variant sites and write to fasta file using snp-sites
    print('Getting variant sites using snp-sites.')
    cmd = path['snpsites'] + ' ' + fasta + \
          flag + ' -o ' + outfile

    os.system(cmd)
    print(outfile, 'written.')

    return(outfile)


def mask_positions(fasta, gff, outdir='.', masked_sites_file=False):
    """Masks positions in GFF file in alignment file in fasta format.

    If fasta is the whole genome alignment fasta used in gubbins 
    and gff is the gubbins output GFF file, then recombinant 
    regions identified by gubbins are masked (on a per-genome basis).
    Optionally returns a text file with list of masked positions in each genome.

    Args:
        fasta: Alignment fasta file in which sites will be masked.
        gff: GFF file containing sites that will be masked in certain genomes
            (ex. gubbins output GFF).
        outdir: Output file directory (default: current working directory).
        masked_sites_file: If true, generates a text file with a list of masked positions in each genome.

    Output:
        Masked whole genome alignment file in FASTA format.
        (Optional) Text file with list of masked positions in each genome.

    Returns:
        Path to masked fasta file.
    """

    # Read in alignment and gff file
    print('Reading ', gff, '.', sep='')
    gff = pd.read_csv(gff, sep='\t', skiprows=2, header=None)
    print('Reading ', fasta, '.', sep='')
    aln = AlignIO.read(fasta, 'fasta')

    # Get indices/positions of recombinant regions identified by gubbins
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

    # Mask indices/positions of recombinant regions identified by gubbins
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

    # Write new FASTA file with recombinant regions masked
    fasta_outfile = outdir + '/' + re.split('/|\.', fasta)[-2] + \
                    '_gubbins_masked.fa'

    print('Writing', fasta_outfile)
    with open(fasta_outfile, 'w') as handle:
        SeqIO.write(new_aln, handle, 'fasta')

    if masked_sites_file:
        # Write text file with list of recombinant sites for each genome
        text_outfile = outdir + '/' + re.split('/|\.', fasta)[-2] + \
                       '_masked_recomb_positions.txt'
        print('Writing', text_outfile)
        with open(text_outfile, 'w') as handle:
            for sample, positions in sample_masked_indices.items():
                line = str(sample) + '\t' + ','.join(map(str, positions)) + '\n'
                handle.write(line)

    return(fasta_outfile)

def count_invar_sites(fasta,gff=None,outdir='.',path={'snpsites':'/nfs/esnitkin/bin_group/anaconda3/bin/snp-sites'}):
    """Counts invariant sites in an alignment file (fasta format).
    
    Gets invariant site count for As, Cs, Gs, and Ts from an alignment file. 
    If gff is not None, positions in the GFF file will be masked before 
    invariant sites are counted.

    Args:
        fasta: Path to alignment file in fasta format.
        gff: [Optional] GFF file of sections of genomes to mask (ex. gubbins output GFF). 
        outdir: Output file directory (default: current working directory).
        path: Dictionary of paths including the path to snp-sites (key: snpsites, value: path)

    Output:
        Text file (*_invar_site_counts.txt) with invariant site counts in the following order: A,C,G,T.
        VCF file of variants (created by snp-sites)
        If GFF path given, masked fasta file (*_gubbins_masked.fa).
    
    Returns:
        Name of text file with invariant site counts (*_invar_site_counts.txt).

    """

    # Mask recombinant regions before counting invariant sites
    if gff is not None:
        aln_file = mask_positions(fasta, gff)
    else:
        aln_file = fasta

    # Count invariant sites in whole genome alignment

    # Read in alignment
    print('Reading ', aln_file, '.', sep='')
    aln = AlignIO.read(aln_file, 'fasta')

    # Get variant positions
    var_site_outfile = outdir + '/' + re.split('/|\.', aln_file)[-2] + \
                       '_snp-sites.vcf'
    print('Getting variant positions using snp-sites.')
    cmd = path['snpsites'] + ' ' + aln_file + \
      ' -v ' + ' -o ' + var_site_outfile

    os.system(cmd)

    positions = []
    with open(var_site_outfile) as f:
        for line in f:
            li=line.strip()
            if not li.startswith("#"):
                positions.append(line.split('\t')[1])

<<<<<<< HEAD
    # Get allele for invariant sites
    invar_sites = []
=======
    # Get ref allele for invariant sites
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be
    for record in aln:
        seq_str = list(str(record.seq))
        for index in positions:
            index = int(index)
            seq_str[index] = ''
<<<<<<< HEAD
        #tmp = ''.join(seq_str)
        if len(invar_sites) == 0:
            invar_sites = seq_str
        else:
            for i,b in enumerate(invar_sites):
                if b is 'N' or b is 'n':
                    invar_sites[i] = seq_str[i]
        invar_counts = Counter(invar_sites)
        if invar_counts['N'] == 0 and invar_counts['n'] == 0:
            break

    del invar_counts['']

    # Get invariant site count for each base
    #print('Counting bases.')
    #invar_counts = Counter(invar_sites)
=======
        invar_sites = ''.join(seq_str)
        break

    # Get invariant site count for each base
    print('Counting bases.')
    invar_counts = Counter(invar_sites)
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be

    # Write base counts to files
    invar_counts_file = outdir + '/' + re.split('/|\.', aln_file)[-2] + \
                          '_invar_site_counts.txt'
    print('Writing ', invar_counts_file, ' (order: A C G T).', sep='')
    with open(invar_counts_file, 'w') as f:
        for base, count in sorted(invar_counts.items()):
<<<<<<< HEAD
            #if base in ['A','a','C','c','G','g','T','t']:
            print(base + ' ' + str(count))
            f.write('%s ' % (count))
            #else:
            #    print(base + ' ' + str(count))
=======
            f.write('%s ' % (count))

>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be
    return(invar_counts_file)
