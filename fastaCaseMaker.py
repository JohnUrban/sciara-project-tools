#!/usr/bin/env python2.7

import sys, argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description="""

DESCRIPTION - change seqs in FASTA to completely uppercase (default) or lowercase.

    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str, default=False, required=True,
                   help='''Path to input file in FASTA format.''')
parser.add_argument('--lower', '-L',
                   action='store_true', default=False,
                   help='''Make sequences lower-case instead of upper-case.''')
args = parser.parse_args()



if args.fasta == "-" or args.fasta == "stdin":
    args.fasta = sys.stdin




for fa in SeqIO.parse(args.fasta, "fasta"):
    print ">"+fa.description
    if args.lower:
        print str(fa.seq).lower()
    else:
        print str(fa.seq).upper()


