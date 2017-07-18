#!/bin/bash

## Read in options
FOFN=$1
LR1=$2
LR2=$3
R1=$4
R2=$5
##EvalThese=$6
CONFIG=$6
SCRIPTS=$7


## paths
EVAL=${SCRIPTS}/eval.ARGS.sh
##CONFIG=${SCRIPTS}/configs/${EvalThese}


## Process
i=0
while read f; do
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  ref=`readlink -f $f` 
  b=`basename $f ${SUFFIX}.fasta`; 
  echo $b; 
  if [ ! -d $b ]; then mkdir $b; fi
  cd $b;
  $EVAL $ref $QOS $SCRIPTS $LR1 $LR2 $R1 $R2 $CONFIG
  cd ../
done < $FOFN



