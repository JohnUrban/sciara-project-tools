#!/bin/bash

if [ $# -eq 0 ]; then echo "
Arg1=Logical true/false -- clean as it works?
Arg2=Config file
Arg3=FOFN - file of filenames -- paths to each assembly FASTA to polish - extension for all fasta files should be .fasta
"; exit; fi

CLEAN=$1
CONFIG=$2
FOFN=$3

PIPELINE=/gpfs_home/jurban/software/kineticsTools/scripts/kineticspipeline.sh


i=0
while read f; do
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  ref=`readlink -f $f` 
  b=`basename $f ${SUFFIX}.fasta`; 
  echo $b; 
  if [ ! -d $b ]; then mkdir $b; fi
  cd $b;
  asmout=${b}${NEWSUFFIX}
  $PIPELINE $ref $QOS $CLEAN $CONFIG
  cd ../
done < $FOFN


