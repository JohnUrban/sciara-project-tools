#!/bin/bash

source ~/software/bwa/source.sh

if $CONVERT_REF_N_TO_ACGT; then
  ## FIRST, ENSURE ALL BASES ARE ACGT -- Ns and Rs (for example) CAUSE ALE TO CRASH
  IUPAC-to-ACGT.py -f $ASM --samplesize 1000000 --convertN --upper > tmp.fasta
  ASM=tmp.fasta
  bwa index $ASM -p $BASE
  rm tmp.fasta
else 
  bwa index $ASM -p $BASE
fi

