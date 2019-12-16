#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION - parse out entries from list of parents.

On most/all systems, this should give same results as:
    grep -f parents.txt file.gff




    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--gff', '-g', '-i', '-f',
                   type= str, default=False, required=True,
                   help='''Path to GFF.''')

parser.add_argument('--parents', '-p',
                   type= str, default=False, required=True,
                   help='''Path to parents.txt.''')

parser.add_argument('--delim', '-d',
                   type= str, default='-', 
                   help='''Parent names for genes are exact matched. When looking for the gene name in parts like mRNA, the part is named genename-ETC.
                        i.e. "genename" and "ETC" are separated by "-" as a delimiter. If you expect it to be different, use this option to change it.''')

parser.add_argument('--namechanger','-n',
                    type= str, default='',
                    help='''Add given string to the front of all ID, Name, and Parent.''')

parser.add_argument('--use_sfx_tree','-T',
                    type= str, default='',
                    help='''Experimental: use suffix tree matching.''')

parser.add_argument('--stree_end_char','-E',
                    type= str, default='$',
                    help='''Experimental: use suffix tree matching. This sets the end character. default = $. It should be something that does not come up otherwise in character set.''')

args = parser.parse_args()


d = {}


with open(args.parents) as f:
    parents = list(set([line.strip() for line in f.readlines()]))



## Function to prune name down to original parent name - useful for the sfx tree approach -- or for both approaches really.
def prune_x(x):
    # The idea here was to anticipate the parts that follow the gene name
    # It seems the names can be broken on -mRNA- for all gene parts to yield the gene name as element 0.
    pass

    
if args.use_sfx_tree:
    from suffix_trees import STree
    # Add a $ at the end of each parents to prevent non-specific matching
    # Need to also add $ to end of query when querying the tree
    #   e.g. if parent = ABC and one sees AB, it would be found in stree,
    #       but not if ABC$ was entered and one saw AB$
    #   Note: $ can be made more specific if need be -- e.g. $$$ or 
    #   e.g. if parent = 
    sfxtree = STree.STree([e+'$' for e in parents])
    def is_in_parents(x, endchar='$'):
        x = prune_name(x)
        
        #stree returns first example found; use find_all to return list of all starts
        #value is -1 if not found or empty list with find_all
        ans = sfxtree.find( x+endchar ) 
        
        for parent in parents:
            if x == parent or x.startswith(parent+delim):
                return True
        return False
else:
    def is_in_parents(x, delim='-'):
        for parent in parents:
            if x == parent or x.startswith(parent+delim):
                return True
        return False
    
#OLD BUG - given parent AB, it will take anythingstarting with AB: AB, ABC, ABWTF, ABETC
##def is_in_parents(x):
##    for parent in parents:
##        if x.startswith(parent):
##            return True
##    return False


def is_in_parents(x, delim='-'):
    for parent in parents:
        if x == parent or x.startswith(parent+delim):
            return True
    return False

##def namechanger(line, d, prefix):
##    keys = d.keys()
##    for key in ('ID','Name','Parent'):
##        if key in keys:
##            d[key] =
##            
##    key = k.strip('_')
##                    if key in ('ID','Name','Parent'):
##                        d[key] = args.namechanger + v
##                    else:
##                        d[key] = v

def return_alt_kv(k, join, prefix, v):
    alt_v = []
    for sub_v in v.split(','):
        alt_v.append(prefix + sub_v)
    new_v = ','.join(alt_v)
    return k + join + new_v

with open(args.gff) as f:
    for line in f:
        d = {}
        if line:
            line = line.strip().split('\t')
            if line[0][0] != '#' and len(line)>3:
                desc = line[8].rstrip(';').split(';')
                altdesc = []
                for e in desc:
                    try:
                        k,v = e.split('=')
                        join = '='
                        #d[k.strip('_')] = v
                    except: #GO terms, etc sep by :
                        k,v = e.split(':')
                        join = ':'
                        #d[k.strip('_')] = v
                    d[k.strip('_')] = v
                    if k in ('ID','Name','Parent'):
                        #OLD BUG: altdesc.append( k + join + args.namechanger + v )
                        altdesc.append( return_alt_kv(k, join, args.namechanger, v) )
                    else:
                        altdesc.append( e )
                altdesc = ';'.join(altdesc)
                if is_in_parents(d['ID'], args.delim):
                    if args.namechanger:
                        line = line[:8] + [altdesc] + line[9:]
                        #pass #line = namechanger(line, d, args.namechanger)
                    print '\t'.join(line)

