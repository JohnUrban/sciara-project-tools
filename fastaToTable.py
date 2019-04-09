#!/usr/bin/env python2.7

import sys, argparse
from Bio import SeqIO
##from Bio.SeqIO import FastaIO
parser = argparse.ArgumentParser(description="""

DESCRIPTION - Take in Fasta, output 2 or 3-column file of name, seq, [other].
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str, default=False, required=True,
                   help='''Path to input file in FASTA format.''')
parser.add_argument('--delim', '-d',
                   type= str, default='_', 
                   help='''Delimiter to join white-space-split fasta headers. Default = _''')

args = parser.parse_args()



if args.fasta == "-" or args.fasta == "stdin":
    args.fasta = sys.stdin



for fa in SeqIO.parse(args.fasta , "fasta"):
    a = fa.description.strip().split()
    name = a[0]
    seq = str(fa.seq)
    if len(a) > 1:
    	other = args.delim.join(a[1:])
	out = [name, seq, other]
    else:
	out = [name, seq]
    print '\t'.join(out)



