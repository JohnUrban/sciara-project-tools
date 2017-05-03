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
#hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS 2> mapreads.err | samtools view -bSh - | samtools sort -T ${PRE} > ${PRE}.bam
#samtools index ${PRE}.bam

## MAP READS
hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS 2> mapreads.err | samtools view -bSh - | samtools sort -n -T ${PRE} > ${PRE}.bam 


# MAP READS AND GET MAPQ STATS 
samtools view -F 4 ${PRE}.bam | awk '{s+=$5}END{print NR,s,s/NR}' > mapq.txt


# GET INFO ON HOW MANY PAIRS MAP TO DIFFERENT CONTIGS
bedtools bamtobed -bedpe -i ${PRE}.bam 2> bamtobed.err | awk '{diff+=$1!=$4; same+=$1==$4; total+=1}END{print same, diff, total}' > pairinfo.txt

# GET INFO ON HOW MANY PAIRS MAP TO DIFFERENT CONTIGS - filtered
# This allows one to get rid of noise from multireads
# However, it also no longer counts pairs where 1 read mapped and 1 did not -- or where 1 read mapq > cutoff and other < cutoff
for i in 2 10 20 30 40; do
  samtools view -bSh -q $i ${PRE}.bam | bedtools bamtobed -bedpe -i - 2> bamtobed.q${i}.err | awk '{diff+=$1!=$4; same+=$1==$4; total+=1}END{print same, diff, total}' > pairinfo.${i}.txt &
done
wait

# CLEAN
if $CLEAN; then rm ${PRE}.bam; fi
