#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
from genomicFileClasses import *

## DEFINE PARSER ARGUMENTS
parser = argparse.ArgumentParser(description="""

    Take in bdg that may have many adjacent intervals with the same score.
    Return collapsed version that can shrink the file size consierable.
    The adjacent intervals with same score are converted to a single interval.
    
    ...
    
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('bdg', metavar='bdg', nargs='+',
                   type= str, 
                   help='''Path(s) to wig file(s).
                        Can handle more than one, though it might not be recommended to use more than one.''')

parser.add_argument('-v', '--verbose', action='store_true')


args = parser.parse_args()



for bdg in args.bdg:
    stdin = False
    if bdg == "stdin" or bdg == "-":
        stdin = True
    Bdg(bdg, collapse_and_exit=True, stdin=stdin)


        
