#!/usr/bin/env python2.7

import sys, argparse
#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

    EXPLAIN.

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('inputs', metavar='inputs', nargs='+',
                   type= str, 
                   help='''...''')






parser.add_argument('--numreadscol', type=int, default=8, help='''1-based. If this is bedpe straight from Sniffles, use ___. If converted it to BEDtools-like bedpe (sometimes I refer to as 'forR.bedpe'), use 8.''')


args = parser.parse_args()
                

if __name__ == "__main__":
    nreadcol = args.numreadscol-1
    with open(args.inputs[0]) as f:
        for line in f:
            if line.startswith('FILLTHISIN'):
                continue
            line = line.strip().split()
            n = int(line[nreadcol])
            for i in range(n):
                print '\t'.join(line[:nreadcol] + ['1'] + line[nreadcol+1:]) 
    
