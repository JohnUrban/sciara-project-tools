#!/usr/bin/env python2.7

import sys, argparse
#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

    EXPLAIN.

   Do this to sniffles bedpe output:
   awk '$1 !~ /^#/ {OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$12,$9,$10,$11}' SV.bedpe > SV.forR.bedpe

   This works on that file directly by default...

   expandSnifflesBedpe.py SV.forR.bedpe > SV.forR-expanded.bedpe

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
    
