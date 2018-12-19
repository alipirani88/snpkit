# PIPELINE TO RUN BEAST ON FLUX

# Import modules
import sys
import os
import argparse
from subprocess import call
import re
<<<<<<< HEAD
from fasta_functions import *
=======
>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be

# Get command line arguments
parser = argparse.ArgumentParser(description='Run BEAST pipeline.')
parser.add_argument('xmls', nargs='*', help='input BEAST xml file(s)')
parser.add_argument('--n_reps', '-n', metavar='N', type=int, help='number of beast runs for each xml', default=3)
parser.add_argument('--invar_sites', '-i', metavar='FILE', help='''
                    option 1: fasta file to  get number of invariant As, Cs, Gs, and Ts not in alignment
                    option 2: text file with number of invariant As, Cs, Gs, and Ts in alignment''')
parser.add_argument('--gubbins_gff', '-g', metavar='GFF', help='''gubbins output GFF file to mask recombinant sites in alignment before counting invariant sites''')
parser.add_argument('--tree', '-t', help='starting tree file in newick format')
parser.add_argument('-scripts_dir', '-d', metavar='DIR', help='directory where all of the scripts are housed', default='/nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/beast')
args = parser.parse_args()
# Rename variables
xmls = args.xmls
n_reps = args.n_reps
invar_sites = args.invar_sites
gubbins_gff = args.gubbins_gff
tree = args.tree
scripts_dir = args.scripts_dir

<<<<<<< HEAD
=======
# Set paths
python = '/nfs/esnitkin/bin_group/anaconda3/bin/python'
mask_vars_path = '/nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/mask_gubbins_variants.py'

>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be
# Get invariant site counts
if invar_sites is not None:
    if not invar_sites.endswith(tuple(['fa','fasta'])):
        # Path to file with invariant site counts
        invar_counts_file = invar_sites
    else:
<<<<<<< HEAD
        # Count invariant sites in whole genome alignment
        invar_counts_file = count_invar_sites(invar_sites, gubbins_gff)
    with open(invar_counts_file) as f:
        invar_counts = f.readline().strip()

=======
        # Mask recombinant regions before counting invariant sites
        if gubbins_gff is not None:
            call([python, mask_vars_path, invar_sites, gubbins_gff])
            aln_file = re.split('/|\.', invar_sites)[-2] + \
                      '_gubbins_masked.fa'
        else:
            aln_file = invar_sites
        # Count invariant sites in whole genome alignment
        call([python, scripts_dir + '/get_num_invariant_sites_v2.py', aln_file])
        invar_counts_file = re.split('/|\.', aln_file)[-2] + \
                          '_invar_site_counts.txt'
    # Read file with invariant site counts
    with open(invar_counts_file) as f:
        invar_counts = f.readline().strip()


>>>>>>> 7b2c92a5aa76895d2a815d185b3d861945add2be
# Final xmls for beast
xmls_for_beast = []

# For each xml input file, modify the xml as desired and start beast
for xml in xmls:

    print('Processing ' + xml + '.')

    # Change prefixes in BEAST xml (more informative when looking at trace files in Tracer)
    print('Adding xml file name as prefix to sequence data in BEAST xml.')
    call([scripts_dir + '/change_prefixes_beast_xml.sh', xml])

    # Get newest xml name
    xml2 = re.split('/|\.', xml)[-2] + '_renamed.xml'

    # Add starting tree, if desired
    if tree is not None:
        print('Adding starting tree to BEAST xml.')
        call([scripts_dir + '/insert_starting_tree_beast.sh', xml2, tree])
        # Get newest xml name
        xml3 = re.split('/|\.', xml2)[-2] + '_st.xml'
    else:
        xml3 = xml2

    # Get invariant site counts
    if invar_sites is not None:
        # Add invariant site counts to xml
        call([scripts_dir + '/add_invariant_sites_beast.sh', xml3, invar_counts]) #[int(x) for x in invar_counts.split(' ')])
        # Get newest xml name
        xml4 = re.split('/|\.', xml3)[-2] + '_invSites.xml'
    else:
        xml4 = xml3

    # Final xml file
    xml_final = xml4
    print('Final xml for input into beast: ' + xml_final + '.')

    # Delete intermediate xmls
    to_remove = []
    for fname in os.listdir('.'):
        b = re.search('.*_renamed_.*xml', fname)
        s = fname == xml_final
        if b and not s:
            print('Removing intermediate file ' + fname + '.')
            os.remove(fname)   

    # Add xml file name to list for beast
    xmls_for_beast.append(xml_final) 

# Generate input file for beast script
with open('input_beast.txt', 'w') as f:
    for x in xmls_for_beast:
        for i in range(n_reps):
            f.write(x + '\n')

#call([scripts_dir + '/generate_input_file_beast_pbs.sh', ' '.join(xmls_for_beast), n_reps])
print('input_beast.txt generated.')

# Generate PBS scripts and submit jobs to flux
print('Generating PBS script(s) and submitting jobs to flux.')
call(['perl', scripts_dir + '/beast_pbs_flux.pl', 'input_beast.txt'])

