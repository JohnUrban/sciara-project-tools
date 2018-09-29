#!/bin/bash

PATH=~/software/reapr/Reapr_1.0.18/:$PATH

## Need: REF, BASE, P, FACHECK, PERFECTMAP, SMALTMAP, PIPELINE, AGGRESSIVE (t/f)

if $FACHECK; then
 echo facheck
 reapr facheck $REF ${BASE}_renamed
fi

if $PERFECTMAP; then
 echo perfectmap
 reapr perfectmap ${BASE}_renamed.fa $R1 $R2 432 perfect
fi

if $SMALTMAP; then
 echo smaltmap
 reapr smaltmap ${BASE}_renamed.fa $R1 $R2 mapped.bam  -n $P
fi

if $PIPELINE; then
 echo pipeline
 if $AGGRESSIVE ; then A="-break a=1"; else A=""; fi
 reapr pipeline ${BASE}_renamed.fa mapped.bam output_directory perfect ${A}
fi

## EXPERIMENTAL
## This used to be in reapr.clean.sh
## However, the obvious shortcoming there is that it is not executed if CLEAN=false
## I moved it here and adjusted the sbatch command to include "G"
## Not yet tested as of Sep29 2018
cd output_directory/
zcat 03.score.per_base.gz | cut -f 3 | awkMean > per-base-mean-score.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py -G $G > broken_assembly.sizestats.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py -t -G $G > broken_assembly.sizestats.csv

