#!/usr/bin/env python

import sys, argparse
from Bio import SeqIO
##from Bio.SeqIO import FastaIO
parser = argparse.ArgumentParser(description="""

DESCRIPTION - format fasta file to have no line wrapping or desired line length.

    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str, default=False, required=True,
                   help='''Path to input file in FASTA format.''')
parser.add_argument('--wraplength', '-L',
                   type=int, default=0,
                   help='''Each fasta line should be <= this length. When given 0, it does no wrapping (entire sequence on 1 line). Default is 0 (no wrapping).''')
args = parser.parse_args()



if args.fasta == "-" or args.fasta == "stdin":
    args.fasta = sys.stdin



for fa in SeqIO.parse(args.fasta , "fasta"):
    print ">"+fa.description
    if args.wraplength <= 0:
        print str(fa.seq)
    else:
        for i in range(0, len(fa), args.wraplength):
            print str(fa.seq[i:i+args.wraplength])



