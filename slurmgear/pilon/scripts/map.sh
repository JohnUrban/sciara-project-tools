#!/bin/bash


#--export=MEM=${MMEM},BT2=${BT2},P=${MTHREADS},R1=${R1},R2=${R2},PRE=#{}

#N=`echo $MEM-2 | bc`
#D=`echo $P*1.5 | bc`
#M=`echo $N/$D | bc`

#if [ $M -lt 1 ]; then M=1; fi
#M=${M}G
#echo Using $M as M parameter.

#bowtie2 -p $P --very-sensitive -N 1 --minins 0 --maxins 1000 -x $BT2 -1 $R1 -2 $R2 | samtools sort --threads $P -m $M > ${PRE}.bam


bowtie2 -p $P --very-sensitive -N 1 --minins 0 --maxins 1000 -x $BT2 -1 $R1 -2 $R2 | samtools sort --threads $P -o ${PRE}.bam
samtools index ${PRE}.bam
