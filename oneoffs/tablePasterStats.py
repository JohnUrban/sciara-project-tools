#!/usr/bin/env python2.7
import sys, os, argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Given set of tables that have shared or overlapping elements in column K,
    go through each table and make lists of values in column J for each element.
    Compute stats for each element.

Use case:
    Two-column tables (name and count) for three RNA-seq replicates.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('table', metavar='table', nargs='*',
                   type= str, 
                   help='''Paths to as many tables as you want. They will be treated as one giant table.
                        This can also be left empty for standard in or specified as - or stdin.''')

parser.add_argument('-n', '--name_column', type=int, default=1,
                    help='''Column where elements/names are found.''')

parser.add_argument('-s', '--score_column', type=int, default=2,
                    help='''Column where scores are found.''')


parser.add_argument('-d', '--delim', type=str, default='\t',
                    help='''Delimiter. Default = tab.''')

args = parser.parse_args()




if len(args.table) == 0 or args.table[0] in ('stdin', '-') or args.table[0].startswith('<('):
    args.table = [sys.stdin]
    def opentable(x):
        return x
    def closetable(x):
        return x
else:
    def opentable(x):
        return open(x)
    def closetable(x):
        return x.close()

def getstats(name, vals):
    return [name] + [np.mean(vals), np.std(vals, ddof=1), np.median(vals), np.min(vals), np.max(vals)]




score_col = args.score_column - 1
name_col = args.name_column - 1

# Collect values for each element
lists = defaultdict(list)
for table in args.table:
    tablefile = opentable(table)
    for line in tablefile:
        line = line.strip().split(args.delim)
        lists[line[name_col]].append( float(line[score_col]) )
    closetable(tablefile)

# Convert to numpy arrays and return stats
for key in lists.keys():
    lists[key] = np.array(lists[key])
    print '\t'.join([str(e) for e in getstats(key, lists[key])])

