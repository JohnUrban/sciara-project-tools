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
args = parser.parse_args()


def get(fh):
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

def table(d1,d2):
    for kmer in sorted(d1.keys()):
        print kmer, d1[kmer], d2[kmer], d1[kmer]/d2[kmer], np.log2(d1[kmer]/d2[kmer]) 

if __name__ == "__main__":

    ## read in A
    A = get(args.fileA)
    B = get(args.fileB)

    ## all kmers
    allkmers = list(set(A.keys()+B.keys()))


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
    table(Anorm, Bnorm)
