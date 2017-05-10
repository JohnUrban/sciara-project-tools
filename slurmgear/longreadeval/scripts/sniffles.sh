#!/bin/bash

echo BAM $BAM
echo BEDPE $BEDPE
echo MINSUPPORT ${MINSUPPORT}
echo OUTPRE $OUTPRE
echo

source ~/software/sniffles/source.sh
sniffles -m $BAM -b $BEDPE.bedpe --min_support $MINSUPPORT


## SVs from Sniffles
## NUM
grep -c -v ^# $$BEDPE.bedpe > ${OUTPRE}numsv

## SUM PREDICTED LENGTHS
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '{s+=$NF}END{print s}' > ${OUTPRE}sumsv

