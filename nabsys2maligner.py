#!/usr/bin/env python2.7
import sys
import argparse
import numpy as np
import re

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in bnx files, returns single line RMAP format (e.g. maligner input).
    
    RMAP format (columns):
    1. MAP NAME (STRING): A unique identifier for the map.
    2. MAP LENGTH (INT): The total length of the map, in base pairs.
    3. NUMBER OF FRAGMENTS (INT): The number of fragments in the map.
    4. FRAGMENT 1 (INT): The length of the first restriction fragment in the map, in base pairs.
    5. FRAGMENT 2 (INT): The length of the second restriction fragment in the map. ... and so on.
    .. .....
    .. .....
    N. Fragment X (INT): where X = N-4+1; length of Xth restriction fragment in the map. ... and so on.

    Fields 4 onwards provide the ordered listing of restriction fragment lengths in bp units.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-f', "--nabsys",
                   type=str, required=True,
                   help='''BNX file from BioNano Irys''')


args = parser.parse_args()



with open(args.nabsys) as f:
    for line in f:
        line = line.strip().split()
        if line[0] == "FragmentUID": #header
            continue
        else:
            fragid = line[0]
            bplen = line[1]
            pos = np.array([float(e) for e in line[2:] if e != 'x'])
            nfrag = str(len(pos))
            fraglens = [str(int(e)) for e in [pos[0]]+list(pos[1:]-pos[:-1])]
            print ('\t').join([fragid, bplen, nfrag] + fraglens)
            
            
            

        
                
