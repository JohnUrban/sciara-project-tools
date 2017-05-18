#!/bin/bash


## FIRST, ENSURE ALL BASES ARE ACGT -- Ns and Rs (for example) CAUSE ALE TO CRASH
IUPAC-to-ACGT.py -f $ASM --samplesize 1000000 --convertN --upper > ${OUTFASTA}

