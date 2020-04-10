#!/bin/bash

# $1 = path to fasta file
# $2 = gubbins (1) or no gubbins (0)
# $3 = flux account to submit to (default: esnitkin_fluxod)
#
# Usage:
# Gubbins:
# gubbins_iqtree_raxml.sh /path/to/whole/genome/alignment/fasta/file 1
# No gubbins:
# gubbins_iqtree_raxml.sh /path/to/snps/only/fasta/file 0

# get prefix for output files
#pref=$(echo $1 | cut -d. -f1 | rev | cut -d/ -f1 | rev)
pref1="${1%.*}"
pref="${pref1##*/}"
echo prefix: $pref

# get working directory
wd=$(dirname $1)
cd $wd

# modules to load (some of these might not be necessary)
modules=$(echo python-anaconda2/201607 biopython fasttree dendropy reportlab RAxML raxml bioperl fastml/gub gubbins openmpi/1.10.2/gcc/4.8.5 gcc/4.8.5)

# Get Directory name where the bash script is located: Changes by Ali 31 Oct 2019
DIRECTORY=`dirname $0`

# get account
if [ -z "$3" ]; then
  echo Will submit jobs to esnitkin_flux.
  acct=esnitkin_flux
else
  echo Will submit jobs to $3.
  acct=$3
fi

mkdir iqtree_results
mkdir raxml_results

# if gubbins
if [ $2 = 1 ]; then
  echo Will run gubbins.
  # gubbins command
  gub=$(echo run_gubbins.py --prefix $pref --threads 12 $1)
  echo $gub > ${pref}_gubbins_command.sh

  # raxml command
  raxml=$(echo mpirun -np 2 raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -N autoMRE -m ASC_GTRGAMMA --asc-corr=lewis -s ../$pref.filtered_polymorphic_sites.fasta -n $pref.filtered_polymorphic_sites_raxML -T 6)
echo $raxml > ${pref}_raxml_command.sh

  # iqtree command
  iqtree=$(echo iqtree -s ../$pref.filtered_polymorphic_sites.fasta -nt AUTO -bb 1000 -m MFP -pre $pref)
  echo $iqtree > ${pref}_iqtree_command.sh
  
  # generate pbs scripts for gubbins, raxml, iqtree
  python $DIRECTORY/pbs_script_maker.py -c ${pref}_gubbins_command.sh -o ${pref}_gubbins.pbs -M "$modules" -a $acct -wd $wd
  python $DIRECTORY/pbs_script_maker.py -c ${pref}_iqtree_command.sh -o ${pref}_iqtree.pbs -M "$modules" -a $acct -wd $wd/iqtree_results
  python $DIRECTORY/pbs_script_maker.py -c ${pref}_raxml_command.sh -o ${pref}_raxml.pbs -M "$modules" -a $acct -wd $wd/raxml_results

  # start gubbins, iqtree, raxml jobs
  echo qsub ${pref}_gubbins.pbs
  gub_qsub=$(qsub ${pref}_gubbins.pbs)
  gub_id=$(echo $gub_qsub | cut -d. -f1)
  echo gubbins job id: $gub_id

  mv ${pref}_iqtree.pbs iqtree_results
  mv ${pref}_iqtree_command.sh iqtree_results
  cd iqtree_results
  echo qsub -W depend=afterok:$gub_id ${pref}_iqtree.pbs
  iqtree_qsub=$(qsub -W depend=afterok:$gub_id ${pref}_iqtree.pbs)

  cd $wd
  mv ${pref}_raxml.pbs raxml_results
  mv ${pref}_raxml_command.sh raxml_results
  cd raxml_results
  echo qsub -W depend=afterok:$gub_id ${pref}_raxml.pbs
  raxml_qsub=$(qsub -W depend=afterok:$gub_id ${pref}_raxml.pbs)

  cd $wd

else
  echo Will not run gubbins.

  echo Finding and removing invariant sites.
  snp-sites -o ${pref}_varSites.fa $1

  # raxml command
  raxml=$(echo mpirun -np 2 raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -N autoMRE -m ASC_GTRGAMMA --asc-corr=lewis -s ../${pref}_varSites.fa -n ${pref}_raxML -T 6)
  echo $raxml > ${pref}_raxml_command.sh

  # iqtree command
  iqtree=$(echo iqtree -s ../${pref}_varSites.fa -nt AUTO -bb 1000 -m MFP+ASC -pre ${pref}_varSites)
  echo $iqtree > ${pref}_iqtree_command.sh

  # generate pbs scripts for iqtree, raxml
  python $DIRECTORY/pbs_script_maker.py -c ${pref}_iqtree_command.sh -o ${pref}_iqtree.pbs -M "$modules" -a $acct -wd $wd/iqtree_results
  python $DIRECTORY/pbs_script_maker.py -c ${pref}_raxml_command.sh -o ${pref}_raxml.pbs -M "$modules" -a $acct -wd $wd/raxml_results

  # start iqtree, raxml jobs

  mv ${pref}_iqtree.pbs iqtree_results
  mv ${pref}_iqtree_command.sh iqtree_results
  cd iqtree_results
  echo qsub ${pref}_iqtree.pbs
  iqtree_qsub=$(qsub ${pref}_iqtree.pbs)
  
  cd $wd
  mv ${pref}_raxml.pbs raxml_results
  mv ${pref}_raxml_command.sh raxml_results
  cd raxml_results
  echo qsub ${pref}_raxml.pbs
  raxml_qsub=$(qsub ${pref}_raxml.pbs)

  cd $wd

fi

