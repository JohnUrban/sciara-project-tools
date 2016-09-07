#!/usr/bin/env python
import sys
import argparse
from collections import defaultdict
from genomicFileClasses import *

## DEFINE PARSER ARGUMENTS
parser = argparse.ArgumentParser(description="""

    Take in fixed step wig, return only data from chromosomes in names file.
    
    ...


    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('wig', metavar='wig', nargs='+',
                   type= str, 
                   help='''Path(s) to wig file(s).
                        Can handle more than one, though it might not be recommended to use more than one.''')

names = parser.add_mutually_exclusive_group()

names.add_argument('--namesfile', '-n', type=str,
                   help='''Path to file with chrom names -- one name per line''')

names.add_argument('--names', '-c', type=str, help='Enter comma-separated names at command line with this flag')


parser.add_argument('-v', '--verbose', action='store_true')


args = parser.parse_args()




names = set()
if args.namesfile:
    if args.namesfile == "-" or args.namesfile == "stdin":
        nfile = sys.stdin
    else:
        nfile = open(args.namesfile)
    for line in nfile:
        names.add(line.rstrip())
elif args.names:
    for name in args.names.split(","):
        names.add(name)
            

for wig in args.wig:
    w = FixedWig(wig, extractnames=names)



## ALT WAY:
##w = FixedWig(wig)
##w.get_wig_for(chroms=sorted(list(names)))


        
