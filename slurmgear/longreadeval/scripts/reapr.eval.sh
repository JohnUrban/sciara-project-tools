#!/bin/bash

PATH=~/software/reapr/Reapr_1.0.18/:$PATH


## Need: REF, BASE, AGGRESSIVE (t/f), BAM


echo facheck
reapr facheck $REF ${BASE}_renamed

if $AGGRESSIVE ; then A="-break a=1"; else A=""; fi

echo pipeline
reapr pipeline ${BASE}_renamed.fa $BAM output_directory ${A}
## -stats s=1 -score l=1000


## do stuff
cd output_directory/
zcat 03.score.per_base.gz | cut -f 3 | awkMean > per-base-mean-score.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py > broken_assembly.sizestats.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py -t > broken_assembly.sizestats.csv
cd ../
