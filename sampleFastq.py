#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict
from fastqTools import *


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Sample reads from fastq file with or without replacement
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fastq', '-f',
                   type= str,
                   help='''Path to fastq file of orphan/unpaired reads. Can be gzipped (must end with .gz).''',
                   default= None)

parser.add_argument('--fastq1', '-1',
                   type= str,
                   help='''Path to fastq file for read#1 in each pair. Can be gzipped (must end with .gz).''',
                   default= None)

parser.add_argument('--fastq2', '-2',
                   type= str,
                   help='''Path to fastq file for read#2 in each pair. Can be gzipped (must end with .gz).''',
                   default= None)

parser.add_argument('--outprefix', '-o',
                   type= str,
                   help='''Prefix for output files. If not specified and sampling unpaired reads,
it will just print to stdout. If not specified and sampling paired reads, it will print to downsampled.1.fastq and downsampled.2.fastq.''',
                   default= None)

parser_replacement = parser.add_mutually_exclusive_group()
parser_replacement.add_argument('--with-replacement', '-w',
                   dest="with_replacement", action="store_true",
                   help='''Sample with replacement.''',
                   default=False)
parser_replacement.add_argument('--without-replacement', '-wo',
                   dest="without_replacement", action="store_true",
                   help='''Sample without replacement.''',
                   default=False)
parser_sample = parser.add_mutually_exclusive_group()
parser_sample.add_argument('--num-reads', '-n', dest='num_reads',
                           type=int,
                           help='''Use this option when doing 'sample with replacement'. Provide integer value.''')
parser_sample.add_argument('--proportion', '-p',
                           type=float,
                           help='''Use this option when doing "sample without replacement". Float between 0 and 1''')
args = parser.parse_args()

#filter bad args -- input fastq
if not args.fastq and not (args.fastq1 and args.fastq2):
    print "Specify either only -f, OR specify -1 and -2"
    quit()
if (args.fastq and (args.fastq1 or args.fastq2)):
    print "Specify either only -f, OR specify -1 and -2"
    quit()


#filter bad args -- replacement
if args.num_reads and args.without_replacement:
    print "For now: Use -p with -wo"
    quit()
if args.proportion and args.with_replacement:
    print "For now: Use -n with -w"
    quit()

# filter outprefix, if given, add .fastq
if args.outprefix:
    args.outprefix += ".fastq"

### Execute:
if args.fastq:
    if args.with_replacement:
        downSampleReadsWithReplacement(args.fastq,args.num_reads, outputFile=args.outprefix)
    elif args.without_replacement:
        downSampleReadsWithoutReplacement(args.fastq, args.proportion, outputFile=args.outprefix)
elif args.fastq1 and args.fastq2:
    if not args.outprefix:
        args.outprefix = "downsampled.fastq"
    if args.with_replacement:
        downSamplePairsWithReplacement(args.fastq1, args.fastq2, args.num_reads, outputFile=args.outprefix)
    elif args.without_replacement:
        downSamplePairsWithoutReplacement(args.fastq1, args.fastq2, args.proportion, outputFile=args.outprefix)
   
