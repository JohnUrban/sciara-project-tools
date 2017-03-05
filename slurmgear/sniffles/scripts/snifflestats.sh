#!/bin/bash

## MARGINSTATS COMMENTED OUT FOR NOW
##requires PB bam, PB fastq, ONT bam, ONT fastq, ASM fasta

echo PBBAM, $PBBAM
echo ONTBAM, $ONTBAM
echo PBFQ, $PBFQ
echo ONTFQ, $ONTFQ
echo ASM, $ASM

D=snifflestats
if [ ! -d $D ]; then mkdir $D; fi
cd $D


## Gets number of alignments -- however each long read can be split into >1 aln, at least with bwa
samtools view -c -F 4 $PBBAM > pacbio.numaln &
samtools view -c -F 4 $ONTBAM > ont.numaln &


##Gets number of uniq read names that aligned
samtools view -F 4 $PBBAM | awk '{print $1}' | sort | uniq | wc -l > pacbio.numuniqaln &
samtools view -F 4 $ONTBAM | awk '{print $1}' | sort | uniq | wc -l > ont.numuniqaln &
wait


## Gets number of alignments + non-aln -- i.e. number of entries in BAM
samtools view -c $PBBAM > pacbio.numentries &
samtools view -c $ONTBAM > ont.numentries &


##Gets number of uniq read names that aligned or did not align --- num uniq entries -- should be equivalent to number of reads in input Fastq file -- should be same for all comparisons
samtools view $PBBAM | awk '{print $1}' | sort | uniq | wc -l > pacbio.numuniqentries &
samtools view $ONTBAM | awk '{print $1}' | sort | uniq | wc -l > ont.numuniqentries &
wait

## Sums MAPQ regardless of everything
samtools view $PBBAM | awk '{s+=$5}END{print s}' > pacbio.sum.mapq &
samtools view $ONTBAM | awk '{s+=$5}END{print s}' > ont.sum.mapq &
wait


##Gets all stats....
#marginStats --printValuePerReadAlignment --identity $PBBAM $PBFQ $ASM > pacbio.marginstats &
#marginStats --printValuePerReadAlignment --identity $ONTBAM $ONTFQ $ASM > ont.marginstats &

wait

## PCTs and RATIOs and Avg MAPQ
for e in ont pacbio; do
 ## PCTs\
 nualn=`cat ${e}.numuniqaln`
 nuent=`cat ${e}.numuniqentries`
 pctaln=`echo $nualn/$nuent | bc -l`
 echo $pctaln > ${e}.pctaln
 ## RATIO
 naln=`cat ${e}.numaln`
 alnratio=`echo $naln/$nualn | bc -l`
 echo $alnratio > ${e}.alnratio
 ## Avg MAPQ
 mapqsum=`cat ${e}.sum.mapq`
 mapqavg=`echo $mapqsum/$naln | bc -l`
 echo $mapqavg > ${e}.avg.mapq
done

## SVs from Sniffles
## NUM
grep -c -v ^# $PBSNIFF > pbnumsv
grep -c -v ^# $ONTSNIFF > ontnumsv
grep -c -v ^# $COMBSNIFF > combnumsv
## SUM PREDICTED LENGTHS
awk '{print $NF}' $PBSNIFF | grep -v pred | awkSum > pbsumsv
awk '{print $NF}' $ONTSNIFF | grep -v pred | awkSum > ontsumsv
awk '{print $NF}' $PBSNIFF | grep -v pred | awkSum > combnumsv

## PUT ALL IN ONE FILE
for f in *; do
  echo -e $f"\t"`cat $f`
done > all-metrics.txt

cd ../
