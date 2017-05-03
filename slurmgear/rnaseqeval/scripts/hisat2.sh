#!/bin/bash

source /users/jurban/data/software/hisat/source.sh

echo P, ${P}
echo READSFOFN, $READSFOFN
echo HIDX, $HIDX
echo PRE, $PRE
echo



## PROCESS READSFOFN TO GET READLISTS
NLINES=`cat $READSFOFN | wc -l`
i=0
while read line; do
  let i++
  if [ $i -eq $NLINES ]; then
    R1LIST+=`echo $line | awk '{print $1}'`;
    R2LIST+=`echo $line | awk '{print $2}'`;
  else
    R1LIST+=`echo $line | awk '{print $1","}'`;
    R2LIST+=`echo $line | awk '{print $2","}'`;
  fi
done < $READSFOFN




## MAP READS
hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS 2> mapreads.err | samtools view -bSh - | samtools sort -T ${PRE} > ${PRE}.bam
samtools index ${PRE}.bam
