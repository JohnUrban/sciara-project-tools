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
                   help='''Path to fastq file.''',
                   default= None, required=True)
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

#filter
if args.num_reads and args.without_replacement:
    print "For now: Use -p with -wo"
    quit()
if args.proportion and args.with_replacement:
    print "For now: Use -n with -w"
    quit()

### Execute: 
if args.with_replacement:
    downSampleReadsWithReplacement(args.fastq,args.num_reads, outputFile=None)
elif args.without_replacement:
    downSampleReadsWithoutReplacement(args.fastq, args.proportion, outputFile=None)
