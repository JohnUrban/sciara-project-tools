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

parser.add_argument('--namechanger','-n',
                    type= str, default='',
                    help='''Add given string to the front of all ID, Name, and Parent.''')
args = parser.parse_args()


d = {}


with open(args.parents) as f:
    parents = list(set([line.strip() for line in f.readlines()]))


def is_in_parents(x):
    for parent in parents:
        if x.startswith(parent):
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
                        altdesc.append( k + join + args.namechanger + v )
                    else:
                        altdesc.append( e )
                altdesc = ';'.join(altdesc)
                if is_in_parents(d['ID']):
                    if args.namechanger:
                        line = line[:8] + [altdesc] + line[9:]
                        #pass #line = namechanger(line, d, args.namechanger)
                    print '\t'.join(line)

