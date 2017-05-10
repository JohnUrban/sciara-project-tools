#!/bin/bash

## MARGINSTATS COMMENTED OUT FOR NOW
##requires PB bam, PB fastq, ONT bam, ONT fastq, ASM fasta

echo PBBAM, $PBBAM
echo ONTBAM, $ONTBAM
echo PBFQ, $PBFQ
echo ONTFQ, $ONTFQ
echo ASM, $ASM

D=stats
#snifflestats
if [ ! -d $D ]; then mkdir $D; fi
cd $D


## Gets number of alignments -- however each long read can be split into >1 aln, at least with bwa
if $MAPPB; then samtools view -c -F 4 $PBBAM > pacbio.numaln & fi
if $MAPONT; then samtools view -c -F 4 $ONTBAM > ont.numaln & fi


##Gets number of uniq read names that aligned
if $MAPPB; then samtools view -F 4 $PBBAM | awk '{print $1}' | sort | uniq | wc -l > pacbio.numuniqaln & fi
if $MAPONT; then samtools view -F 4 $ONTBAM | awk '{print $1}' | sort | uniq | wc -l > ont.numuniqaln & fi
wait


## Gets number of alignments + non-aln -- i.e. number of entries in BAM
if $MAPPB; then samtools view -c $PBBAM > pacbio.numentries & fi
if $MAPONT; then samtools view -c $ONTBAM > ont.numentries & fi


##Gets number of uniq read names that aligned or did not align --- num uniq entries -- should be equivalent to number of reads in input Fastq file -- should be same for all comparisons
if $MAPPB; then samtools view $PBBAM | awk '{print $1}' | sort | uniq | wc -l > pacbio.numuniqentries & fi
if $MAPONT; then samtools view $ONTBAM | awk '{print $1}' | sort | uniq | wc -l > ont.numuniqentries & fi
wait

## Sums MAPQ regardless of everything
if $MAPPB; then samtools view $PBBAM | awk '{s+=$5}END{print s}' > pacbio.sum.mapq & fi
if $MAPONT; then samtools view $ONTBAM | awk '{s+=$5}END{print s}' > ont.sum.mapq & fi
wait


##Gets all stats....
#marginStats --printValuePerReadAlignment --identity $PBBAM $PBFQ $ASM > pacbio.marginstats &
#marginStats --printValuePerReadAlignment --identity $ONTBAM $ONTFQ $ASM > ont.marginstats &
#wait

## PCTs and RATIOs and Avg MAPQ
if $MAPPB && $MAPONT; then E="ont pacbio"; 
elif $MAPPB; then E=pacbio; 
elif $MAPONT; then E=ont; 
fi

for e in $E ; do
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
## NUM and SUM PREDICTED LENGTHS
if $SNIFFLESPB; then
 grep -c -v ^# $PBSNIFF > pbnumsv
 grep -v ^# $PBSNIFF | awk -v "s=0" '{s+=$NF}END{print s}' > pbsumsv
fi

if $SNIFFLESONT; then
 grep -c -v ^# $ONTSNIFF > ontnumsv
 grep -v ^# $ONTSNIFF | awk -v "s=0" '{s+=$NF}END{print s}' > ontsumsv
fi

if $SNIFFLESCOMBINED; then
 grep -c -v ^# $COMBSNIFF > combnumsv
 grep -v ^# $COMBSNIFF | awk -v "s=0" '{s+=$NF}END{print s}' > combsumsv
fi


###awk '{print $NF}' $PBSNIFF | grep -v pred | awkSum > pbsumsv
###awk '{print $NF}' $ONTSNIFF | grep -v pred | awkSum > ontsumsv
###awk '{print $NF}' $COMBSNIFF | grep -v pred | awkSum > combsumsv


## PUT ALL IN ONE FILE
if $MAPONT; then
 PRE=ont
 ORDEREDONT="${PRE}.numaln  ${PRE}.numentries ${PRE}.numuniqaln ${PRE}.numuniqentries ${PRE}.sum.mapq ${PRE}.pctaln ${PRE}.alnratio ${PRE}.avg.mapq"
fi
if $MAPPB; then
 PRE=pacbio
 ORDEREDPB="${PRE}.numaln  ${PRE}.numentries ${PRE}.numuniqaln ${PRE}.numuniqentries ${PRE}.sum.mapq ${PRE}.pctaln ${PRE}.alnratio ${PRE}.avg.mapq"
fi

if $SNIFFLESONT; then ORDEREDNUMSV[1]=ontnumsv; ORDEREDSUMSV[1]=ontsumsv; fi  
if $SNIFFLESPB; then ORDEREDNUMSV[2]=pbnumsv; ORDEREDSUMSV[2]=pbsumsv; fi  
if $SNIFFLESCOMBINED; then ORDEREDNUMSV[3]=combnumsv; ORDEREDSUMSV[3]=combsumsv; fi  

ORDEREDFILES="$ORDEREDONT $ORDEREDPB ${ORDEREDNUMSV[@]} ${ORDEREDSUMSV[@]}"

for f in ${ORDEREDFILES}; do
  echo -e $f"\t"`cat $f`
done > all-metrics.txt

cd ../
