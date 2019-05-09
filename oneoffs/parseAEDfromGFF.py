#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION - parse out AED table from Maker GFF.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--gff', '-g', '-i', '-f',
                   type= str, default=False, required=True,
                   help='''Path to GFF.''')
args = parser.parse_args()


d = {}
with open(args.gff) as f:
    for line in f:
        if line:
            line = line.strip().split('\t')
            if line[0][0] != '#' and len(line)>3:
                if line[2] == 'mRNA':
                    desc = line[8].rstrip(';').split(';')
                    for e in desc:
                        try:
                            k,v = e.split('=')
                            d[k.strip('_')] = v
                        except: #GO terms, etc sep by :
                            k,v = e.split(':')
                            d[k.strip('_')] = v
                    out = [d['Parent'], d['ID'], d['AED'], d['eAED']]
                    print '\t'.join(out)

