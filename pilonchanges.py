#!/usr/bin/env python
import gzip, sys, argparse
import numpy as np
from Bio import SeqIO


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in pilon.changes file. Reports shit.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument("-i", "--inputfile",
                   type= str, default=False, required=True,
                   help='''Input pilon.changes file.''')
args = parser.parse_args()


def complement(seq):
    ''' assumes uppercase '''
    c=""
    for b in seq:
            if b == "A": c += "T"
            elif b == "C": c += "G"
            elif b == "G": c += "C"
            elif b == "T": c += "A"
    return c    

def revcomp(seq):
    return complement(seq)[-1::-1]

def a(line):
    ''' line from changes file'''
    pass
    
insertions = 0
deletions = 0
subs = 0
inslens = {}
dellens = {}
sublens = {}

