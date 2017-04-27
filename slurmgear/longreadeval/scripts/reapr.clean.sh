#!/bin/bash


rm *

cd output_directory/
zcat 03.score.per_base.gz | cut -f 3 | awkMean > per-base-mean-score.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py > broken_assembly.sizestats.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py -t > broken_assembly.sizestats.csv
 
rm -r 00.* 01.* 02.* 03.* 04.*

cd ../
