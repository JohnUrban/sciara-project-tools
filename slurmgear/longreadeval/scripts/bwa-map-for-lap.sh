#!/bin/bash

##NEED: MTRHEADS, TYPE, BWAIDX, FASTQ
 
source ~/software/bwa/source.sh
bwa mem -t $MTHREADS -a -x $TYPE $BWAIDX $FASTQ | samtools sort -n -T $TYPE.forlap --threads $MTHREADS --output-fmt sam -o $TYPE.lap.sam
