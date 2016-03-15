#!/usr/bin/env python

import sys, argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description="""

DESCRIPTION - version 24Mar2015
    
    Change/modify name of each fasta entry in a file with MANY fasta entries.

    Fasta names usually follow ">name" format.

    Some may follow ">name otherinfo".
    This can screw up some programs.
    A way to keep only the part up to the first white space is:
       $ fasta_name_changer.py -f file.fa -k 1
    Or keeping only the second block of info as the name:
       $ fasta_name_changer.py -f file.fa -k 2
    etc.

    Add nameinfo in front of name or selected part of name:
       $ fasta_name_changer.py -f file.fa -F prefix
       $ fasta_name_changer.py -f file.fa -k 1 -F prefix

    Add nameinfo in  back of name or selected part of name:
       $ fasta_name_changer.py -f file.fa -B prefix
       $ fasta_name_changer.py -f file.fa -k 1 -B prefix

    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str, default=False, required=True,
                   help='''Path to input file in FASTA format.''')
parser.add_argument('--keep', '-k',
                   type=int, default=False,
                   help='''Break name up by whitespace. Keep only kth chunk (e.g. 1).''')
parser.add_argument('--front', '-F',
                   type=str, default=False,
                   help='''Add this info to front of name (operation occurs after -k)''')
parser.add_argument('--back', '-B',
                   type=str, default=False,
                   help='''Add this info to back of name (operation occurs after -k)''')


args = parser.parse_args()
if args.keep:
    chunk = args.keep-1
## SeqIO automatically takes ">name" from ">name other info" with entry.name
## entry.description gives whole name "name other info" -- removes ">" in front
for entry in SeqIO.parse(args.fasta, "fasta"):
    name = entry.description 
    if args.keep:
        i = args.keep-1
        name = name.split()[chunk]
    if args.front:
        name = args.front+name
    if args.back:
        name = name+args.back
    name = ">"+name
    sys.stdout.write(name + "\n")
    sys.stdout.write(str(entry.seq) + "\n")


