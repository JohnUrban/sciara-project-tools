#!/bin/bash

source ~/software/bwa/source.sh

BAM=${TYPE}.bam
bwa mem -t $MTHREADS -M -x $TYPE $BWAIDX $FASTQ | samtools sort -T $TYPE --threads $MTHREADS -o ${BAM}


## stats

## ont2d --> ont
if [ $TYPE == "ont2d" ] && [ -z $PRE ]; then PRE=ont; elif [ -z $PRE ]; then PRE=${TYPE}; fi


## Gets number of alignments -- however each long read can be split into >1 aln, at least with bwa
samtools view -c -F 4 $BAM > ${PRE}.numaln &

##Gets number of uniq read names that aligned
samtools view -F 4 $BAM | awk '{print $1}' | sort | uniq | wc -l > ${PRE}.numuniqaln &

## wait for last 2
wait


## Gets number of alignments + non-aln -- i.e. number of entries in BAM
samtools view -c $BAM > ${PRE}.numentries &

##Gets number of uniq read names that aligned or did not align --- num uniq entries -- should be equivalent to number of reads in input Fastq file -- should be same for all comparisons
samtools view $BAM | awk '{print $1}' | sort | uniq | wc -l > ${PRE}.numuniqentries &

## wait for last 2
wait



## Sums MAPQ regardless of everything
samtools view $BAM | awk '{s+=$5}END{print s}' > ${PRE}.sum.mapq &
wait


## PCTs and RATIOs and Avg MAPQ
## PCTs\
nualn=`cat ${PRE}.numuniqaln`
nuent=`cat ${PRE}.numuniqentries`
pctaln=`echo $nualn/$nuent | bc -l`
echo $pctaln > ${PRE}.pctaln

## RATIO
naln=`cat ${PRE}.numaln`
alnratio=`echo $naln/$nualn | bc -l`
echo $alnratio > ${PRE}.alnratio

## Avg MAPQ
mapqsum=`cat ${PRE}.sum.mapq`
mapqavg=`echo $mapqsum/$naln | bc -l`
echo $mapqavg > ${PRE}.avg.mapq


