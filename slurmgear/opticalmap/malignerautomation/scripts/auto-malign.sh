#!/bin/bash

if [ $# -eq 0 ]; then echo "
Arg1=Logical true/false -- clean as it works?
Arg2=Config file
Arg3=FOFN - file of filenames -- paths to each assembly FASTA to polish - extension for all fasta files should be .fasta
Arg4=FOFN -- for smoothed bionano maps
Arg5=REC_ENZ (BssSI)
Arg6=REC_SEQ (CACGAG)
Arg7=scripts dir
"; exit; fi

CLEAN=$1
CONFIG=$2
ASMFOFN=$3
MAPSFOFN=$4 
REC_ENZ=$5
REC_SEQ=$6
SCRIPTS=$7

MALPIPE=${SCRIPTS}/maligner-pipeline.sh


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
  $MALPIPE $ref $MAPSFOFN $QOS $CLEAN $CONFIG $REC_ENZ $REC_SEQ $SCRIPTS
  cd ../
done < $ASMFOFN


