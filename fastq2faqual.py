#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Given a fastq file - generate fasta and/or qual file

    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('--fastq', '-f',
                   type= str,
                   help='''Path to fastq file.''')

parser.add_argument('--fa', action='store_true',
                   help='''Generate fasta file.''')

parser.add_argument('--qual', action='store_true',
                   help='''Generate qual file.''')

parser.add_argument('--out', type=str, default=False,
                   help='''Output prefix. Default: same as fastq''')

args = parser.parse_args()

assert args.fa or args.qual

if not args.out:
    args.out = args.fastq.split("/")[-1].split(".")[0]

if args.fa:
    SeqIO.convert(args.fastq, "fastq", args.out+".fasta", "fasta")
if args.qual:
    SeqIO.convert(args.fastq, "fastq", args.out+".qual", "qual")
