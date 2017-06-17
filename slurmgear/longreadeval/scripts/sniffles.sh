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
grep -c -v ^# $BEDPE.bedpe > numsv

## SUM PREDICTED LENGTHS
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '{s+=$NF}END{print s}' > sumsv


## TYPES
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="DEL" {s+=1}END{print s}' > numdel
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="DUP" {s+=1}END{print s}' > numdup
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="INS" {s+=1}END{print s}' > numins
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="INV" || $11=="DEL/INV" || $11=="INVDUP" || $11=="INV/INVDUP" {s+=1}END{print s}' > numinv
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="TRA" {s+=1}END{print s}' > numtra


## SV TYPES I ENCOUNTERED
#  13419 DEL
#      1 DEL/INV
#   1268 DUP
#   1037 INS
#    167 INV
#     65 INVDUP
#      5 INV/INVDUP
#   9185 TRA
