#!/bin/bash


ASMDIR=asms

###EVAL=/gpfs/scratch/jurban/male-ilmn/long_read_evals/scripts/eval.ARGS.sh
EVAL=/gpfs/scratch/jurban/male-ilmn/long_read_evals/scripts/eval.ARGS.buscoOnly.sh 

SUBSET="platanus.dbg2olc.pball.ontmol.quiver7x platanus.dbg2olc.pball.quiver7x platanus.dbg2olc.pbfilt.ont2d.quiver7x platanus.dbg2olc.pbfilt.quiver7x"

i=0
CLEAN=true
for d in $SUBSET; do 
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  f=${ASMDIR}/$d.fasta
  ref=`readlink -f $f` 
  echo $d; 
  cd $d; pwd
  if [ -d busco ]; then rm -r busco; mkdir busco; fi
  $EVAL $ref $QOS ###$CLEAN
  cd ../
done 

pwd
