#!/bin/bash

##NEEDS GENOME SIZE, G

rm *


## THE COMMENTED-OUT LINES ARE NOW PART OF repar.eval.sh as of 29Sep2018 -- not yet tested though.
#cd output_directory/
#zcat 03.score.per_base.gz | cut -f 3 | awkMean > per-base-mean-score.txt
#faSize -detailed 04.break.broken_assembly.fa | asm-stats.py -G $G > broken_assembly.sizestats.txt
#faSize -detailed 04.break.broken_assembly.fa | asm-stats.py -t -G $G > broken_assembly.sizestats.csv
 
rm -r 00.* 01.* 02.* 03.* 04.*

cd ../
