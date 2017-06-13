#!/bin/bash

##$RUN $SCRIPTS $CONFIG $CLEAN $ASMFOFN 

if [ $# -eq 0 ]; then echo "
Arg1 = SCRIPTS location
Arg2 = CONFIG location
Arg3 = CLEAN logical 
Arg4 = ASM FOFN
"; exit; fi

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
ASMFOFN=$4

PIPELINE=${SCRIPTS}/buscov3-pipeline.sh


i=0
while read f; do
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  REF=`readlink -f $f` 
  b=`basename $f .fasta`; 
  echo buscov3 $b; 
  if [ ! -d $b ]; then mkdir $b; fi
  cd $b;
  $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF
  cd ../
done < $FOFN


