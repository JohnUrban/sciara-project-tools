#!/usr/bin/env python2.7
import sys, os, argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Given the output of "bedtools cluster" as input,

return only entries that are part of clusters with X or more elements.

""", formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('bed', metavar='bed', nargs='*',
                   type= str, 
                   help='''Path to BED file.
                        This can also be left empty for standard in or specified as - or stdin.''')

parser.add_argument('-k', '--clusterColumn', type=int, required=True,
                    help='''Specify the column that cluster IDs are in.''')

parser.add_argument('-m', '--min', type=int, default=-1,
                    help='''Specify minimum number of elements a cluster must have.''')

parser.add_argument('-M', '--max', type=int, default=1000000000000,
                    help='''Specify maximum number of elements a cluster must have.''')

args = parser.parse_args()




## Functions
def bed_to_dict(bed, col, keytype=int):
    # clusters are integers so default keytype is int... however, later usages may use strings for example so I left it as an option
    bedDict = defaultdict(list)
    for line in bed:
        line = line.strip().split()
        bedDict[keytype(line[col])].append( line )
    return bedDict
            


def read_bed(f, col):
    if f in ('stdin', '-') or f.startswith('<('):
        fh = sys.stdin
        return bed_to_dict(fh, col)
    else:
        with open(f) as fh:
            return bed_to_dict(fh, col)


def filterBed(bedDict, m, M):
    for cluster in sorted(bedDict.keys()):
        n = len(bedDict[cluster])
        if n >= m and n <= M:
            for line in bedDict[cluster]:
                print '\t'.join(line)


def main(args):
    filterBed( bedDict=read_bed(args.bed[0], args.clusterColumn - 1), m=args.min, M=args.max )




## Execute
try:
    main(args)
except IOError:
    pass
                


