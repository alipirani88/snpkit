#!/bin/bash

# creates file with names of certain files in directory based on what would be listed with ls
# ex: generate_input_file_beast_pbs.sh '*.xml' 3 
# In the example above, an input_beast.txt file would be generated with all of the xml files in the directory listed 3 times. This will cause BEAST to run in triplicate for each xml file.
# if using wildcard, put it in quotes

num=$2

export num;

ls $1 | perl -ne 'print $_ x $ENV{num}' > input_beast.txt
