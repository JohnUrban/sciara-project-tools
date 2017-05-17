#!/bin/bash

PATH=~/software/sciaratools/sciara-project-tools/alenano/src/:$PATH
PATH=~/software/alenano/src/:/users/jurban/software/sciaratools/sciara-project-tools/:$PATH

## FIRST, ENSURE ALL BASES ARE ACGT -- Ns and Rs (for example) CAUSE ALE TO CRASH
IUPAC-to-ACGT.py -f $REF --samplesize 1000000 --convertN --upper > tmp.fasta

REF=tmp.fasta

ALE ${BAM} ${REF} ${BASE}.ALE.txt >> ${BASE}.ale.err

rm tmp.fasta

if $CLEAN; then bash $SCRIPTS/ale.clean.sh ; fi

