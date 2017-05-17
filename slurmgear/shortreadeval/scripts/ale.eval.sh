#!/bin/bash

module load ale/20140120
PATH=~/software/frcbam/FRC_align/bin/:/users/jurban/software/sciaratools/sciara-project-tools/:$PATH


###NOTE:
#  -I/--minins <int>  minimum fragment length (0)
#  -X/--maxins <int>  maximum fragment length (500)
#  Should probably have -X as 600-800


## FIRST, ENSURE ALL BASES ARE ACGT -- Ns and Rs (for example) CAUSE ALE TO CRASH
IUPAC-to-ACGT.py -f $REF --samplesize 1000000 --convertN --upper > tmp.fasta

REF=tmp.fasta

ALE ${BAM} $REF ${BASE}.ALE.txt >> $BASE.ale.err

rm tmp.fasta
