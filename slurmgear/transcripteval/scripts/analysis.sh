#!/bin/bash

for i in $(seq $NJOBS); do 
  f=${BLASTDIR}/*.${i}.blastout; 
  awk -v "i=$i" 'OFS="\t" {p+=$3*$4/100.0; b+=$10}END{print i,p,b}' $f; 
done > indiv.txt

awk '{p+=$2; b+=$3}END{print p, b}' indiv.txt > total.txt


