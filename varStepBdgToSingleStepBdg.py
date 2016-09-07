#!/usr/local/bin/python

import sys

if len(sys.argv) == 1:
    f = sys.stdin
else:
    f = open(sys.argv[1], 'r')

for line in f:
    line = line.rstrip().split()
    chrom = line[0]
    start = int(line[1])
    end = int(line[2])
    score = float(line[3])
    while end > start:
        sys.stdout.write("\t".join([str(e) for e in [chrom, start, start+1, score]])+"\n")
        start += 1
