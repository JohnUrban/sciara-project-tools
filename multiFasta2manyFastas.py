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

parser.add_argument('--dir', '-d',
                   type= str, default=False,
                   help='''Directory to write output files to. Default: pwd. ''',)

parser.add_argument('--suffix', '-s',
                   type= str, default=False,
                   help='''Suffix word to add just before .fa. Default: "".
Note: if you want a dot or dash before it, it needs to be included.''',)


args = parser.parse_args()

infa = args.fasta
if infa == "-" or infa == "stdin":
    infa = sys.stdin
    
outdir = ""
if args.dir:
    outdir = args.dir
    if not outdir.endswith("/"):
        outdir += "/"

sfx = ""
if args.suffix:
    sfx = args.suffix
    
if args.numnames:
    i=0
    word = args.numnames+"_"
    for fa in SeqIO.parse(infa, "fasta"):
        i+=1
        SeqIO.write(fa, outdir+word+str(i)+sfx+".fa", "fasta")
else:
    for fa in SeqIO.parse(infa, "fasta"):
        SeqIO.write(fa, outdir+fa.name+sfx+".fa", "fasta")


