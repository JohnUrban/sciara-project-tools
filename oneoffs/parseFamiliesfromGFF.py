#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION - parse out entries from list of parents.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--gff', '-g', '-i', '-f',
                   type= str, default=False, required=True,
                   help='''Path to GFF.''')

parser.add_argument('--parents', '-p',
                   type= str, default=False, required=True,
                   help='''Path to parents.txt.''')
args = parser.parse_args()


d = {}


with open(args.parents) as f:
    parents = list(set([line.strip() for line in f.readlines()]))


def is_in_parents(x):
    for parent in parents:
        if x.startswith(parent):
            return True
    return False

with open(args.gff) as f:
    for line in f:
        if line:
            line = line.strip().split('\t')
            if line[0][0] != '#' and len(line)>3:
                desc = line[8].rstrip(';').split(';')
                for e in desc:
                    try:
                        k,v = e.split('=')
                        d[k.strip('_')] = v
                    except: #GO terms, etc sep by :
                        k,v = e.split(':')
                        d[k.strip('_')] = v
                if is_in_parents(d['ID']):
                    print '\t'.join(line)

