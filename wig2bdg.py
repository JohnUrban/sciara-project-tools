#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
from genomicFileClasses import *

## DEFINE PARSER ARGUMENTS
parser = argparse.ArgumentParser(description="""

    Take in fixed step wig, return fixed step bdg.
    
    ...
    
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('wig', metavar='wig', nargs='+',
                   type= str, 
                   help='''Path(s) to wig file(s).
                        Can handle more than one, though it might not be recommended to use more than one.''')

parser.add_argument('-e', '--make_span_eq_step', action='store_true',
                    help='''Some times a wig file has a span size that is different from step size (usually <= step).
                        When vizualizing, this creates white space betwen bars in chart.
                        Using this flag will use step size for both step and span when converting to bedGraph.''')
                    
parser.add_argument('-v', '--verbose', action='store_true')


args = parser.parse_args()



for wig in args.wig:
    FixedWig(wig, wig2bdg=True, make_span_eq_step=args.make_span_eq_step)


        
