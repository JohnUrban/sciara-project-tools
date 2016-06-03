#!/usr/bin/env python
import sys

def gc(x):
    c=0
    for b in x:
        if b in 'GC':
            c+=1
    return 100.0*c/len(x)

seqs = []

for line in open(sys.argv[1]):
    seq = line.strip().split()
    for s in seq:
        seqs.append(s)

for seq in seqs:
    print gc(seq)
