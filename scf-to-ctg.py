#!/usr/bin/env python
import sys
import argparse
import numpy as np
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Take in FASTA file.

Break sequence at Ns.

Return FASTA file with each scf broken into sequence between Ns (gaps).

Naming will be:
scfname_1
.
.
.
scfname_n

If there are n-1 such gaps.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fasta', '-f',
                   type=str, required=True,
                   help='''Path to input fasta file''')

args = parser.parse_args()


for fa in SeqIO.parse(args.fasta, 'fasta'):
    scf = str(fa.seq)
    ctglist = scf.split('N')
    ctgnum = 0
    Ncnt = 0
    for candidate in ctglist:
        if candidate:
            ctgnum += 1
            length = len(candidate)
            print ">" + fa.name + "_ctg" + str(ctgnum) + "len="+str(length) + " after " + str(Ncnt) + " Ns"
            print candidate
            Ncnt = 0
        else:
            Ncnt += 1
