#!/usr/bin/env python
from Bio import SeqIO
from collections import defaultdict
import sys
import numpy as np


## ONT fastqs have Q values spanning 0-90.
## This is breaking programs that do not expect it to exceed 64 (e.g. ALE).
## Can re-compile those programs with different expectations.... maybe.
## Can also re-vamp fastq to not exceed 64 -- give it cieling.
helpmsg='''
Usage:
enforceQualCeiling.py fastqfile[required] ceiling[default:64]
'''
if len(sys.argv) == 1:
    print helpmsg
    quit()

if sys.argv[1] == "-" or sys.argv[1] == "stdin":
    sys.argv[1] = sys.stdin

if len(sys.argv) == 3:
    ceiling = int(sys.argv[2])
else:
    ceiling = 64



for fq in SeqIO.parse(sys.argv[1], "fastq"):
    fq.letter_annotations['phred_quality'] = np.array(fq.letter_annotations['phred_quality'])
    fq.letter_annotations['phred_quality'][fq.letter_annotations['phred_quality'] > ceiling] = ceiling
    SeqIO.write(fq, sys.stdout, "fastq")
