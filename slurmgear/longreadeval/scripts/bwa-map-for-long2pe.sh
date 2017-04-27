#!/bin/bash


source ~/software/bwa/source.sh

##bwa mem -t $MTHREADS -M -x $TYPE $BWAIDX $R1 $R2 | samtools sort -T $TYPE.long2pe --threads $MTHREADS -o $TYPE.long2pe.bam


bwa mem -t $MTHREADS -x $TYPE $BWAIDX $R1 $R2 | samtools sort -T $TYPE.long2pe --threads $MTHREADS -o $TYPE.long2pe.bam

