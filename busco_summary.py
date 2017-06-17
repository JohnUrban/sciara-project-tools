#!/usr/bin/env python

## takes in fulltable.tsv output from buscov3

import sys
import numpy as np
from collections import defaultdict

fulltable = open(sys.argv[1])

total = 0
dup = 0
sc = 0
status = defaultdict(int)
scores = defaultdict(list)
lengths = defaultdict(list)

for line in fulltable:
    if line[0] == '#':
        continue
    else:
        line = line.strip().split()
        total += 1
        status[line[1]] += 1
        if line[1] != 'Missing':
            scores[line[0]].append( float(line[5]) )
            lengths[line[0]].append( float(line[6]) )
fulltable.close()

scorelist = []
lenlist = []

for gene in scores.keys():
    scorelist.append( np.mean(scores[gene]) )
    lenlist.append( np.mean(lengths[gene]) )
    if len(scores[gene]) > 1:
        dup += 1
    elif len(scores[gene]) == 1:
        sc += 1  ## also counts fragmented

print 'Complete_Single_copy\t' + str(status['Complete'])
print 'Complete_multi_copy\t' + str(dup)
print 'Complete_total\t' + str(status['Complete']+dup)
print 'Fragmented\t' + str(status['Fragmented'])
print 'Total_complete_and_fragmented\t' + str(status['Complete']+dup+status['Fragmented'])
print 'Missing\t' + str(status['Missing'])
print 'Total\t' + str( status['Complete']+dup+status['Fragmented']+status['Missing']  )
print 'meanscore\t' + str(np.mean(scorelist))
print 'meanlen\t' + str(np.mean(lenlist))

