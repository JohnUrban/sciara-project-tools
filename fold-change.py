#!/usr/bin/env python2.7
import os, sys, gzip
import argparse
from collections import defaultdict
from math import log10, log
import numpy as np




#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """




John Urban (2015, 2016, 2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)



parser.add_argument('-A', '--fileA', type=str, 
                    help = '''.''')

parser.add_argument('-B', '--fileB', type=str, 
                    help = '''.''')

parser.add_argument('-C', '--pseudocount', type=int, default=0,
                    help = '''.''')

parser.add_argument('-EA', '--excludeA', type=int, default=0,
                    help = '''Exclude kmers from a fileA with raw counts less than this. Default 0.
                            This operation is independent of fileB.
                            It is mainly to exclude error kmers from fileA.
                            Also see --excludeB.
                            To exclude kmers from both that only appear in 1, see --intersect.''')
parser.add_argument('-EB', '--excludeB', type=int, default=0,
                    help = '''Exclude kmers from a fileB with raw counts less than this. Default 0.
                            This operation is independent of fileA.
                            It is mainly to exclude error kmers from fileB.
                            Also see --excludeA.
                            To exclude kmers from both that only appear in 1, see --intersect.''')
parser.add_argument('-I', '--intersect', action='store_true', default=False,
                    help = '''The default behavior is to analyze the union of kmers sets between fileA and fileB.
                            Thus any kmer that appears in one and not the other will get a 0-count (or pseudocount)
                            for the file it was absent in.
                            Kmers that apear in only 1 file may be errors or may not be of interest for another reason.
                            This tells it to analyze only kmers that were in both files (after --excludeA/B options).''')


args = parser.parse_args()


def get(fh, mincount=0):
    d = defaultdict(int)
    stdin = False
    if fh in ('-','stdin') or fh.startswith('<('):
        f = sys.stdin
        stdin=True
    elif fh.endswith('.gz'):
        f = gzip.open(fh, 'rb')
    else:
        f = open(fh)
        
    for line in f:
        line = line.strip().split()
        count = int(line[1])
        if count >= mincount:
            d[line[0]] = int(line[1])
    if not stdin:
        f.close()
    return d

def median(d):
    return np.median(d.values())

def mad(d,med):
    return np.median(np.absolute(np.array(d.values()) - med))
    

def med_norm(d, med):
    newd = {}
    med = float(med)
    for kmer in d.keys():
        newd[kmer] = d[kmer]/med
    return newd

def mad_norm(d, med, mad):
    newd = {}
    med = float(med)
    for kmer in d.keys():
        newd[kmer] = (d[kmer]-med)/mad
    return newd

def table(d1,d2,id1,id2):
    for kmer in sorted(d1.keys()):
        print kmer, d1[kmer], d2[kmer], d1[kmer]/d2[kmer], np.log2(d1[kmer]/d2[kmer]), id1[kmer], id2[kmer], id1[kmer]/id2[kmer] 

if __name__ == "__main__":

    ## read in A
    A = get(args.fileA, mincount=args.excludeA)
    B = get(args.fileB, mincount=args.excludeB)

    ## all kmers
    if args.intersect:
        allkmers = set(A.keys()).intersection(set(B.keys()))
    else:
        allkmers = list(set(A.keys()).union(set(B.keys())))

    ## ensure each kmer in both
    for kmer in allkmers:
        A[kmer] += args.pseudocount
        B[kmer] += args.pseudocount

    ## obtain medians
    med_A = median(A)
    med_B = median(B)

    Anorm = med_norm(A, med_A)
    Bnorm = med_norm(B, med_B)

    ## print
    table(Anorm, Bnorm, A, B)
