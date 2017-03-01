#!/bin/bash

PATH=~/software/frcbam/FRC_align/bin/:/users/jurban/software/sciaratools/sciara-project-tools/:$PATH


 if [ ! -f reads.bam ]; then
   bowtie2 -p $MTHREADS -q --very-sensitive -N 1 -x $BT2 -1 $R1 -2 $R2 2> $BASE.mapreads.err | samtools sort --threads $MTHREADS -o reads.bam
 fi
 
 if [ ! -f reads.bam.bai ]; then 
   samtools index reads.bam
 fi
