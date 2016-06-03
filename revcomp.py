#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from helper_functions import *

parser = argparse.ArgumentParser(description="""

Given a fasta/fastq file, return reverse complements.

    """, formatter_class= argparse.RawTextHelpFormatter)


filetype = parser.add_mutually_exclusive_group()
filetype.add_argument('-fa', type=str, default="",
                   help='''Fasta input. If stdin, either: leave empty''')
filetype.add_argument('-fq', type=str, default="",
                   help='''Fastq input.''')
parser.add_argument("-out",
                    type=int, default=False,
                    help='''Output file prefix. Default: stdout.''')
special = parser.add_mutually_exclusive_group()
special.add_argument("-even",
                    action='store_true', default=False,
                    help='''Only revcomp even numbered entries.''')
special.add_argument("-odd",
                    action='store_true', default=False,
                    help='''Only revcomp odd numbered entries.''')
special.add_argument("-only",
                    type=str, default=False,
                    help='''Only revcomp entry numbers supplied as comma-separated list.''')
args = parser.parse_args()


if args.only:
    entrynums = set([int(e) for e in args.only.strip().split(',')])
    
if args.fa:
    fastxFile = args.fa
    fastx = "fasta"
elif args.fq:
    fastxFile = args.fq
    fastx = "fastq"
if fastxFile in ("","-","stdin"):
    fastxFile = sys.stdin
if args.out:
    args.out = open(args.out,"r")
else:
    args.out = sys.stdout

i = 0
for record in SeqIO.parse(fastxFile, fastx):
    i += 1
    if args.only:
        if i in entrynums:
            args.out.write(">"+record.description + "\n" + str(record.seq.reverse_complement()) + "\n")
        else:
            args.out.write(">"+record.description + "\n" + str(record.seq) + "\n")
    elif args.even:
        if i%2 == 0:
           args.out.write(">"+record.description + "\n" + str(record.seq.reverse_complement()) + "\n")
        else:
            args.out.write(">"+record.description + "\n" + str(record.seq) + "\n")
    elif args.odd:
        if i%2 == 1:
           args.out.write(">"+record.description + "\n" + str(record.seq.reverse_complement()) + "\n")
        else:
            args.out.write(">"+record.description + "\n" + str(record.seq) + "\n")
    else: ##ALL
        args.out.write(">"+record.description + "\n" + str(record.seq.reverse_complement()) + "\n")
args.out.close()
