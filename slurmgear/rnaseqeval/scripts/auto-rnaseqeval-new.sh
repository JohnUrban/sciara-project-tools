#!/bin/bash


if [ $# -eq 0 ]; then echo "
Arg1=scripts dir
Arg2=Config file
Arg3=Logical true/false -- clean as it works?
Arg4=ASMFOFN - file of filenames -- paths to each assembly FASTA - extension for all fasta files should be .fasta
Arg5=READSFOFN - tab-sep 2 column -- col1=R1 paths ; col2=R2 paths -- each line has corresponding pairs
"; exit; fi


SCRIPTS=$1
CONFIG=$2
CLEAN=$3
ASMFOFN=$4
READSFOFN=$5



PIPELINE=${SCRIPTS}/rnaseqeval-pipeline-new.sh

i=0
while read f; do
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=ccmb-condo; i=0; else QOS=biomed-condo; fi
  REF=`readlink -f $f` 
  if [[ "$f" == *.fa ]]; then B=`basename $f .fa`; 
  elif [[ "$f" == *.fasta ]]; then B=`basename $f .fasta`; fi
  echo $B; 
  if [ ! -d $B ]; then mkdir $B; fi
  cd $B;
  $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF $READSFOFN
  cd ../
done < $ASMFOFN


