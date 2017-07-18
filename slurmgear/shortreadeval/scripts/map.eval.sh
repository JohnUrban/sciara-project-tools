#!/bin/bash

PATH=~/software/frcbam/FRC_align/bin/:/users/jurban/software/sciaratools/sciara-project-tools/:$PATH


 if [ ! -f reads.bam ]; then
   bowtie2 -p $MTHREADS -q --very-sensitive -N 1 -x $BT2 -1 $R1 -2 $R2 2> $BASE.mapreads.err | samtools sort --threads $MTHREADS -o reads.bam
 fi

 #quality check
 count=`samtools quickcheck -vvvvv reads.bam 2>&1 | grep -c "EOF marker is absent"`
 if [ $count -gt 0 ]; then
   rm reads.bam
   bowtie2 -p $MTHREADS -q --very-sensitive -N 1 -x $BT2 -1 $R1 -2 $R2 2> $BASE.mapreads.err | samtools sort --threads $MTHREADS -o reads.bam
 fi

 if [ ! -f reads.bam.bai ]; then 
   samtools index reads.bam
 fi
