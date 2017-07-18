#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
import random

parser = argparse.ArgumentParser(description="""

Take in FASTA file.

Add some Ns at random places.


    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fasta', '-f', type=str, required=True, help='''Path to input fasta file''')

parser.add_argument('--num', '-n', type=int, default=1, help='''How many gaps to add. Adds up to this many, but not necessarily this many. Default: 1.''')
parser.add_argument('--minlen', '-m', type=int, default=100, help='''Minimum gap length. Default: 100.''')
parser.add_argument('--maxlen', '-M', type=int, default=1000, help=''''Maximum gap length. Default: 1000.''')
parser.add_argument('--minctglen', '-c', type=int, default=100000, help='''Minimum length contig can be to introduce gaps. Default: 100000 (100kb).''')

args = parser.parse_args()

gapsadded=0

for fa in SeqIO.parse(args.fasta, 'fasta'):
    seq = str(fa.seq)
    length = len(seq)
    if length >= args.minctglen and gapsadded < args.num:
        start = random.randint(0,length-args.maxlen)
        end = start + random.randint(args.minlen, args.maxlen)
        gaplen = end-start
        newseq = seq[:start]
        newseq += 'N'*gaplen
        newseq += seq[end:]
        gapsadded += 1
    else:
        newseq = seq
    print ">" + fa.description 
    print newseq

