#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import sys, argparse, re
from collections import defaultdict
from itertools import product
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""
    January 7, 2020. John Urban.

    Merges histogram files I've made from NuPoP analysis.
    
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('nupop_hists', metavar='nupop_hists', nargs='+',
                   type= str, 
                   help='''Path to fastx.''')

#parser.add_argument('-n', '--nlines', type=int, required=True, help='''Number of lines a correct hist file should have....''')
parser.add_argument('-c', '--getcounts', action='store_true', default=False, help='''Get counts column in newer NuPoP hists...''')

args = parser.parse_args()




class NuPoP(object):
    #"pos","pstart","occ","nl","aff"
    def __init__(self, fh=None, getcounts=False):
        self.fh = fh
        self.getcounts = getcounts
        self.pstart = defaultdict(float)
        self.nucocc = defaultdict(float)
        self.nl = defaultdict(float)
        self.affinity = defaultdict(float)
        self.counts = defaultdict(float)
        if self.fh is not None:
            self._process()


        
    def haslines(self, n=10):
        return self.nlines > n
    def _process(self):
        with open(self.fh) as f:
            self.lines = f.readlines()
        self.nlines = len(self.lines)
        for line in self.lines:
            line = line.strip().split()
            pos = int(line[0])
            self.pstart[pos] += float(line[1])
            self.nucocc[pos] += float(line[2])
            self.nl[pos] += float(line[3])
            self.affinity[pos] += float(line[4])
            if self.getcounts:
                self.counts[pos] += float(line[5])
                
    def add(self, other):
        for pos in other.affinity.keys():
            self.pstart[pos] += other.pstart[pos]
            self.nucocc[pos] += other.nucocc[pos]
            self.nl[pos] += other.nl[pos]
            self.affinity[pos] += other.affinity[pos]
            if self.getcounts:
                self.counts[pos] += other.counts[pos]

    def __str__(self):
        out=''
        for pos in sorted(self.affinity.keys()):
            l = [pos, self.pstart[pos], self.nucocc[pos], self.nl[pos], self.affinity[pos]]
            if self.getcounts:
                l.append( self.counts[pos] )
            out += '\t'.join([str(e) for e in l]) + '\n'
        return out.strip()


nupopsum = NuPoP(fh=None, getcounts=args.getcounts)
for hist in args.nupop_hists:
    nupop = NuPoP(fh=hist, getcounts=args.getcounts)
    if nupop.haslines():
        nupopsum.add(nupop)

print nupopsum

