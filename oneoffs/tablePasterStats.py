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

Returns:
    Name, mean, stdev, COV, Zscore_min, Zscore_median, Zscore_max, min, median, max.

    If 3 tables are given as in the use case, the min/med/max will essentially be the 3 input values in order from lowest to highest.

    COV = coefficent of variation = stdev/mean.
    Note when mean = 0, stdev/1 is returned.

    Zscores are (x-mean)/stdev.
    Note when stdev = 0, (x-mean)/1 is returned.

    COV and Zscores give you an idea of spread and outliers.
    
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
    mu = np.mean(vals)
    sd = np.std(vals, ddof=1)
    med = np.median(vals)
    m = np.min(vals)
    M = np.max(vals)
    sd1 = 1 if sd == 0 else sd
    mu1 = 1 if mu == 0 else mu
    
    return [name] + [mu, sd, sd/mu1, (m-mu)/sd1, (med-mu)/sd1, (M-mu)/sd1, m, med, M] 




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

