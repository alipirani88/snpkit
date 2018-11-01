#!/bin/bash

## Changes prefixes in beast xml file
## Specifically, adds prefixes to the sequence id name and discrete trait name based on xml file name (the prefix is the xml file name)
## Takes as input: beast xml file
## Usage: change_prefixes_beast.sh beast.xml 
## Note: in example above, the prefix would be 'beast'

ids=$(grep -A1 '<data' $1 | grep 'id=' | grep -o -P '(?<=").*(?=")')
len_ids=$(echo $ids | wc -w)
seqs=$(echo $ids | cut -d' ' -f1)
if [ $len_ids == 2 ]; then 
  discrete=$(echo $ids | cut -d' ' -f2)
fi

prefix=$(echo $1 | grep -o -P '.*(?=\.xml)')

sed "s/${seqs}/${prefix}_${seqs}/g" $1 > ${prefix}_renamed.xml
if [ $len_ids == 2 ]; then
  sed -ibak "s/${discrete}/${prefix}_${discrete}/g" ${prefix}_renamed.xml
  rm ${prefix}_renamed.xmlbak
fi
sed -ibak "s/EBSP/${prefix}_EBSP/g" ${prefix}_renamed.xml
rm ${prefix}_renamed.xmlbak
