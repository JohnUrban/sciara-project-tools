#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION - parse BEDPE files from sniffles to make more PAF-like.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bedpe', '-i', '-f', '-b',
                   type= str, default=False, required=True,
                   help='''Path to Sniffles BedPe.''')

parser.add_argument('--genome', '-g', 
                   type= str, default=False, required=True,
                   help='''Path to BEDtools-like genome file (chr, length).''')

args = parser.parse_args()



genome = {}


with open(args.genome) as f:
    for line in f:
        line = line.strip().split()
        genome[line[0]] = line[1]

getidx = {"chr1":0,"start1":1,"end1":2,"chr2":3,"start2":4,"end2":5,"name":6,"score":7,"strand1":8,"strand2":9,"type":10,"nreads":11,"best_chr1":12,"best_start":13,"best_chr2":14,"best_stop":15,"length":16,"filter":17}

def bedpeToPaf(genome, getidx, line):
    get = lambda x: line[getidx[x]]
    out1 = [get("chr1"), genome[get("chr1")], get("start1"), get("end1"),get("strand1")]
    out2 = [get("chr2"), genome[get("chr2")], get("start2"), get("end2"),get("strand2")]
    out3 = [get("type"), get("nreads")]
    return out1 + out2 + out3

with open(args.bedpe) as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split()
            print '\t'.join( bedpeToPaf(genome, getidx, line) )
        
