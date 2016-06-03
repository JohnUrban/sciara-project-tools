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
parser.add_argument('--numnames', '-n', type=str, default=False,
                    help=''' By default, fasta filenames are the fasta entry names.
This option allows you to specify a word for the filename instead. Then each entry will be named WORD_N.fa
numbered according to its position in the fastafile.''')


args = parser.parse_args()


if args.numnames:
    i=0
    word = args.numnames+"_"
    for fa in SeqIO.parse(args.fasta, "fasta"):
        i+=1
        SeqIO.write(fa, word+str(i)+".fa", "fasta")
else:
    for fa in SeqIO.parse(args.fasta, "fasta"):
        SeqIO.write(fa, fa.name+".fa", "fasta")


