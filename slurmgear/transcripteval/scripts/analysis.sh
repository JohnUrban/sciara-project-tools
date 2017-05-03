#!/bin/bash

# matchlength, bitscore, length

for i in $(seq $NJOBS); do 
  f=${BLASTDIR}/*.${i}.blastout; 
  n=`awk '{print $1}' $f | sort | uniq | wc -l`
  awk -v "i=$i" -v "n=$n" 'OFS="\t" {p+=$3*$4/100.0; b+=$10; l+=$4}END{printf "%d\t%f\t%f\t%f\t%f\n", i,p,b,l,n}' $f; 
done > indiv.txt

awk '{p+=$2; b+=$3; l+=$4; n+=$5}END{printf "%f\t%f\t%f\t%f\n", p, b, l, n}' indiv.txt > total.txt


