#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Given a fasta/fastq file return 2 column table: seqname seqlength.

    """, formatter_class= argparse.RawTextHelpFormatter)


filetype = parser.add_mutually_exclusive_group()
filetype.add_argument('-fa', type=str, default="",
                   help='''Fasta input. If stdin, either: leave empty''')
filetype.add_argument('-fq', type=str, default="",
                   help='''Fastq input.''')
parser.add_argument("-out",
                    type=int, default=False,
                    help='''Output file prefix. Default: stdout.''')

parser.add_argument("--minlen",
                    type=int, default=0,
                    help='''Use this to only report sequences >= int given. Default: 0.''')

parser.add_argument("--maxlen",
                    type=int, default=30000000000,
                    help='''Use this to report sequences <= int given. Default: 30 billion.''')

args = parser.parse_args()


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
for record in SeqIO.parse(fastxFile, fastx):
    if len(record) >= args.minlen and len(record) <= args.maxlen:
        args.out.write(record.name + "\t" + str(len(record)) + "\n")
args.out.close()
