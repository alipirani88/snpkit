# Functions to manipulate beast xml files

# Functions:
#    add_starting_tree

# Import modules
import glob
import re

def add_starting_tree(xml,tree):
    """Creates a new beast xml file with a starting tree added.

    Args:
        xml: beast xml file (or regular expression such as *xml) to add starting tree to
        tree: tree file with tree to add to beast xml file

    Output:
        Beast xml file with starting tree added.
        Returns list of output xmls.
    """    
   
    # Read in tree 
    with open(tree, 'r') as f:
        st = f.read().replace('\n', '')
    
    # Get all xmls
    xmls = glob.glob(xml)
    
    # Initialize list of output xml names
    outxmls = []

    # For each xml, add starting tre
    for fi in xmls:
        # Read in xml
        with open(fi, 'r') as f:
            x = f.read()
        # Get name of output file
        outxml = re.sub('.xml','_st.xml',fi)
        # Add output file name to list of output xmls
        outxmls.append(outxml)
        # Remove random tree, add Newick tree
        x = re.sub('id=\"RandomTree','id=\"NewickTree',x)
        tosub = 'util.TreeParser\" IsLabelledNewick=\"true\" offset=\"0\" newick=\"' + st + '\"'
        x = re.sub('evolution.tree.RandomTree\" estimate=\"false\"',tosub,x)
        x = re.sub('<populationModel id=\"ConstantPopulation.*populationModel>','',x,flags=re.DOTALL)
        # Write xml with starting tree to new file
        with open(outxml,'w') as f:
            #print('Writing ' + outxml + '.') 
            f.write(x)
    # Return new xml file names
    return(outxmls)
        

