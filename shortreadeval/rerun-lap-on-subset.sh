#!/bin/bash


ASMDIR=asms

###EVAL=/gpfs/scratch/jurban/male-ilmn/long_read_evals/scripts/eval.ARGS.sh
EVAL=/gpfs/scratch/jurban/male-ilmn/long_read_evals/scripts/eval.ARGS.laponly.sh 


SUBSET="falcon.seed25.relaxed.pbfilt.ont2d.quiver7x"

i=0
CLEAN=true
for d in $SUBSET; do 
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  f=${ASMDIR}/$d.fasta
  ref=`readlink -f $f` 
  echo $d; 
  cd $d; pwd
  $EVAL $ref $QOS ###$CLEAN
  cd ../
done 

pwd
