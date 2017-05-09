#!/bin/bash

# matchlength, bitscore, length

echo job identity_len sum_bit len num_queries_ge_1_hit | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' > indiv.txt
for i in $(seq $NJOBS); do 
  f=${BLASTDIR}/*.${i}.blastout; 
  n=`awk '{print $1}' $f | sort | uniq | wc -l`
  awk -v "i=$i" -v "n=$n" 'OFS="\t" {p+=$3*$4/100.0; b+=$10; l+=$4}END{printf "%d\t%f\t%f\t%f\t%f\n", i,p,b,l,n}' $f; 
done >> indiv.txt

echo identity_len sum_bit len num_queries_ge_1_hit | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' > total.txt
awk '{p+=$2; b+=$3; l+=$4; n+=$5}END{printf "%f\t%f\t%f\t%f\n", p, b, l, n}' indiv.txt >> total.txt


## clean.txt is a file made by parent script if CLEAN parameter set to true
if [ -f ../clean.txt ]; then 
  rm -r ${BLASTDIR}/
else 
  echo nocleaning..... $PWD
  ls
fi
