#!/usr/bin/env python
from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Take a file with multiple fasta entries in it and split it into
    a file for each fasta entry with the name of the fa file being
    the name of the entry.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str,
                   help='''Input file in fasta format containing multiple sequences. ''',
                   required= True)


args = parser.parse_args()


for fa in SeqIO.parse(args.fasta, "fasta"):
    SeqIO.write(fa, fa.name+".fa", "fasta")

