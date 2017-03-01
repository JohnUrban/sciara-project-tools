#!/bin/bash


ASMDIR=asms

EVAL=/gpfs/scratch/jurban/male-ilmn/long_read_evals/scripts/eval.ARGS.sh


i=0

for f in ${ASMDIR}/*fasta; do 
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  ref=`readlink -f $f` 
  b=`basename $f .fasta`; 
  echo $b; 
  mkdir $b; 
  cd $b;
  $EVAL $ref $QOS
  cd ../
done 

