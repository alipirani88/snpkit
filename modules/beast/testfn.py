a=['b','1','a','2']
no_numbers = True
for i in a:
    try:
        int(i)
        no_numbers = False
        break
    except:
        continue

from scotti_functions import *

test = add_prefixes('one_facil_subtree_hosts.csv', 'one_facil_subtree_dates.csv', 'one_facil_subtree_hostTimes.csv')

#test = add_prefixes('one_facil_subtree_hosts.renamed.csv', 'one_facil_subtree_dates.renamed.csv', 'one_facil_subtree_hostTimes.renamed.csv')

print(test[0])
