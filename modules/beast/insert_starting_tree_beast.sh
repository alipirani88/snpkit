#!/bin/bash

# Add starting tree to BEAST xml file (instead of starting with a random tree)
# Arguments: 1) xml file
#	     2) tree file

tree=$(less $2)

for i in $(ls $1); do

outfile=$(echo $i | cut -d. -f1)
outfile=${outfile}_st.xml
echo $outfile

sed -e 's/id="RandomTree/id="NewickTree/g' -e "s/evolution.tree.RandomTree\" estimate=\"false\"/util.TreeParser\" IsLabelledNewick=\"true\" offset=\"0\" newick=\""$tree"\"/g" -e '/<populationModel id="ConstantPopulation/,/populationModel>/d' $i > $outfile

done
