#!/usr/bin/env python
import argparse
from Bio import SeqIO
from collections import defaultdict
from numpy.random import choice
parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Take a multi-fasta file (.fa file with >=1 entry) where sequences may contain IUPAC letters,
    convert non-ACGT into ACGT.

    ABYSS, for example, uses non-ACGT letters. This breaks some other software.

    non-ACGT IUPAC letters include:
    K = G or T
    M = A or C
    S = G or C
    R = A or G
    W = A or T
    Y = C or T
    B = C or G or T
    D = A or G or T
    H = A or C or T
    V = A or C or G
    N = A or C or G or T

    N will be kept as N by default as these usually indicate gaps.
    
    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--fasta', '-f', required=True,
                   type=str, default=False,
                   help='''Provide path to fasta file''')
parser.add_argument('--samplesize', '-n',
                   type=int, default=1000000, 
                   help='''Provide number of bases to sample to determine ACGT frequencies.
If 0, uniform probabilities are assumed.
If -1, the entire file will be used.
Default = 1 million. ''')
parser.add_argument('--pseudocount', '-c',
                   type=int, default=0, 
                   help='''You may want to provide a pseudocount (e.g. 1) if doing a small sample size.
Default is 0 to allow for correctly learning frequencies from sequences missing a certain base.''')
parser.add_argument('--convertN', '-N', 
                   action='store_true', default=False,
                   help='''Convert Ns into ACGT. Default: False.''')
args = parser.parse_args()





## functions
def probs(ltrs, freq):
    s = 0
    for b in ltrs:
        s += freq[b]
    p = []
    for b in ltrs:
        p.append( freq[b]/s )
    return p

if args.convertN:
    safe_ltrs = "AaCcGgTt"
    other_ltrs = "KkMmSsRrWwYyBbDdHhVvNn"
else:
    safe_ltrs = "AaCcGgTtNn"
    other_ltrs = "KkMmSsRrWwYyBbDdHhVv"


ctr = 0
pc = args.pseudocount
if args.samplesize == 0:
    freq = {"A":1.0, "C":1.0, "G":1.0, "T":1.0}
else:
    if args.samplesize < 0:
        args.samplesize = float('inf')
    freq = {"A":0.0+pc, "C":0.0+pc, "G":0.0+pc, "T":0.0+pc}
    for fa in SeqIO.parse(args.fasta, 'fasta'):
        for b in str(fa.seq):
            if ctr <= args.samplesize:
                ctr += 1
                B = b.upper()
                if B in "ACGT":
                    freq[B] += 1.0
            else:
                break


translate = {}
translate['K'] = lambda : choice(['G','T'], p=probs('GT',freq))
translate['k'] = lambda : choice(['g','t'], p = probs('GT',freq))
translate['M'] = lambda : choice(['A','C'], p = probs('AC',freq))
translate['m'] = lambda : choice(['a','c'], p = probs('AC',freq))
translate['S'] = lambda : choice(['G','C'], p = probs('GC',freq))
translate['s'] = lambda : choice(['g','c'], p = probs('GC',freq))
translate['R'] = lambda : choice(['A','G'], p = probs('AG',freq))
translate['r'] = lambda : choice(['a','g'], p = probs('AG',freq))
translate['W'] = lambda : choice(['A','T'], p = probs('AT',freq))
translate['w'] = lambda : choice(['a','t'], p = probs('AT',freq))
translate['Y'] = lambda : choice(['C','T'], p = probs('CT',freq))
translate['y'] = lambda : choice(['c','t'], p = probs('CT',freq))
translate['B'] = lambda : choice(['C','G','T'], p = probs('CGT',freq))
translate['b'] = lambda : choice(['c','g','t'], p = probs('CGT',freq))
translate['D'] = lambda : choice(['A','G','T'], p = probs('AGT',freq))
translate['d'] = lambda : choice(['a','g','t'], p = probs('AGT',freq))
translate['H'] = lambda : choice(['A','C','T'], p = probs('ACT',freq))
translate['h'] = lambda : choice(['a','c','t'], p = probs('ACT',freq))
translate['V'] = lambda : choice(['A','C','G'], p = probs('ACG',freq))
translate['v'] = lambda : choice(['a','c','g'], p = probs('ACG',freq))
translate['N'] = lambda : choice(['A','C','G', 'T'], p = probs('ACGT',freq))
translate['n'] = lambda : choice(['a','c','g','t'], p = probs('ACGT',freq))



for fa in SeqIO.parse(args.fasta, 'fasta'):
    newseq = ''
    for ltr in str(fa.seq):
        if ltr not in safe_ltrs:
            newseq += translate[ltr]()
        else:
            newseq += ltr
    print ">"+fa.description
    print newseq
        
        
