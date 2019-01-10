#!/bin/bash

# input: 
# (1) beast xml file (ex. beauti output) 
# (2) number of invariant A's, C's, G's, and T's (in that order)
# ex. bash add_invariant_sites_beast.sh beast.xml 100 120 130 140

fout=$(echo $(echo $1 | cut -d. -f1)_invSites.xml)
dataf=$(awk 'NR==1,/data id=\"/' $1 | tail -n1 | cut -d'"' -f2)
awk 'NR==1,/data id=\"/{sub(/data id=\"/, "data id=\"orig_")} 1' $1 > $fout

sed -ibak "1,/data>/s/data>/&\n\n<data id='${dataf}' \n spec='FilteredAlignment' \n filter='-' \n data='@orig_${dataf}' \n constantSiteWeights='${2} \n ${3} ${4} ${5}'\/> <!-- A C G T -->/" $fout

rm ${fout}bak
