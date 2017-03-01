#!/bin/bash

echo PBBAM, $PBBAM
echo ONTBAM, $ONTBAM
echo PBFQ, $PBFQ
echo ONTFQ, $ONTFQ
echo ASM, $ASM

D=snifflestats
if [ ! -d $D ]; then mkdir $D; fi
cd $D


## Gets number of 2d, temp, and comp alignments 
samtools view -F 4 $ONTBAM | grep -c ^2d > ont.2d.numaln &
samtools view -F 4 $ONTBAM | grep -c ^template > ont.template.numaln &
samtools view -F 4 $ONTBAM | grep -c ^complement > ont.complement.numaln &


##Gets number of uniq 2d, temp, and comp read names that aligned
samtools view -F 4 $ONTBAM | awk '{print $1}' | sort | uniq | grep -c ^2d > ont.2d.numuniqaln &
samtools view -F 4 $ONTBAM | awk '{print $1}' | sort | uniq | grep -c ^template > ont.template.numuniqaln &
samtools view -F 4 $ONTBAM | awk '{print $1}' | sort | uniq | grep -c ^complement > ont.complement.numuniqaln &
wait


##Gets number of uniq read names that aligned or did not align for 2d, temp, comp
samtools view $ONTBAM | awk '{print $1}' | sort | uniq | grep -c ^2d > ont.2d.numuniqentries &
samtools view $ONTBAM | awk '{print $1}' | sort | uniq | grep -c ^template > ont.template.numuniqentries &
samtools view $ONTBAM | awk '{print $1}' | sort | uniq | grep -c ^complement > ont.complement.numuniqentries &
wait

## Sums MAPQ regardless of everything -- 2d, temp, comp
samtools view $ONTBAM | grep ^2d | awk '{s+=$5}END{print s}' > ont.2d.sum.mapq &
samtools view $ONTBAM | grep ^template | awk '{s+=$5}END{print s}' > ont.template.sum.mapq &
samtools view $ONTBAM | grep ^complement | awk '{s+=$5}END{print s}' > ont.complement.sum.mapq &

wait


cd ../
