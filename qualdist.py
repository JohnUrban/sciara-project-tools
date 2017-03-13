#!/bin/usr/env python
from Bio import SeqIO
from collections import defaultdict
import sys
## adapted from poretools qualdist


qual_count = defaultdict(int)
total_nucs = 0

for fq in SeqIO.parse(sys.argv[1], "fastq"):
    for q in fq.letter_annotations['phred_quality']:
        qual_count[q] += 1
        total_nucs += 1
for q in qual_count:
    print '\t'.join(str(s) for s in [chr(q+33), q, qual_count[q], 
        total_nucs, float(qual_count[q]) / float(total_nucs)])
