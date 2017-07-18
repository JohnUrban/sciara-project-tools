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
filetype.add_argument('-c', '--cmdline', type=str, default="",
                   help='''Raw sequence (no fa/fq formatting) from command-line. If given directly, needs to be comma-separated. If coming from stdin, needs to be one per line.''')
##filetype.add_argument('-fq', type=str, default="",
##                   help='''Fastq input.''')
parser.add_argument("-out",
                    type=int, default=False,
                    help='''Output file prefix. Default: stdout.''')
special = parser.add_mutually_exclusive_group()
special.add_argument("-even",
                    action='store_true', default=False,
                    help='''Only revcomp even numbered entries. Return odd-numbered as is.''')
special.add_argument("-odd",
                    action='store_true', default=False,
                    help='''Only revcomp odd numbered entries. Return even-numbered as is.''')
special.add_argument("-only",
                    type=str, default=False,
                    help='''Only revcomp entry numbers supplied as comma-separated list. Return others as is.''')
parser.add_argument("-no_other",
                    action='store_true', default=False,
                    help='''When using -even, -odd, or -only, the default behavior still returns 'other' sequences as is
                    w/o revcomping. Setting this says only return even, odd, or only arguments with revcomps -- nothing else, no other.''')
args = parser.parse_args()

## Process arguments
if args.only:
    entrynums = set([int(e) for e in args.only.strip().split(',')])
    
if args.fa:
    fastxFile = args.fa
    fastx = "fasta"
elif args.fq:
    fastxFile = args.fq
    fastx = "fastq"
elif args.cmdline and args.cmdline in ("","-","stdin"):
    fastxFile = args.cmdline
else:
    fastxFile = None


## seq info coming from stdin?
if fastxFile in ("","-","stdin"):
    fastxFile = sys.stdin

## Output file?
if args.out:
    args.out = open(args.out,"r")
else:
    args.out = sys.stdout


## Define "records"
if args.cmdline and fastxFile == sys.stdin:
    falist = fastxFile
elif args.cmdline:
    falist = args.cmdline.split(",")
elif args.fa or args.fq:
    falist = SeqIO.parse(fastxFile, fastx)


def getrevcomp(record, args):
    if args.fa or args.fq:
        return ">"+ record.description + "\n" + str(record.seq.reverse_complement()) + "\n"
    elif args.cmdline:
        return revcomp(record) + "\n"

def returnidentity(record, args):
    if args.fa or args.fq:
        return ">"+record.description + "\n" + str(record.seq) + "\n"
    elif args.cmdline:
        return record + "\n"


## Execute
i = 0
for record in falist:
    i += 1
    if args.only:
        if i in entrynums:
            args.out.write( getrevcomp(record, args) )
        elif not args.no_other:
            args.out.write( returnidentity(record, args) )
    elif args.even:
        if i%2 == 0:
           args.out.write( getrevcomp(record, args) )
        elif not args.no_other:
            args.out.write( returnidentity(record, args) )
    elif args.odd:
        if i%2 == 1:
           args.out.write( getrevcomp(record, args) )
        elif not args.no_other:
            args.out.write( returnidentity(record, args) )
    else: ##ALL
        args.out.write( getrevcomp(record, args) )
args.out.close()
