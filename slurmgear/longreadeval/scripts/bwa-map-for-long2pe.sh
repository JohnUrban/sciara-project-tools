#!/bin/bash


source ~/software/bwa/source.sh

######bwa mem -t $MTHREADS -M -x $TYPE $BWAIDX $R1 $R2 | samtools sort -T $TYPE.long2pe --threads $MTHREADS -o $TYPE.long2pe.bam
echo MTHREADS ${MTHREADS} 
echo TYPE ${TYPE}
echo BWAIDX ${BWAIDX}
echo R1 $R1
echo R2 $R2
echo

bwa mem -t $MTHREADS -x $TYPE $BWAIDX $R1 $R2 | samtools sort -T $TYPE.long2pe --threads $MTHREADS -o $TYPE.long2pe.bam

