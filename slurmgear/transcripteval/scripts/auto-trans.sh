#!/bin/bash

##$RUN $SCRIPTS $CONFIG $CLEAN $ASMFOFN $TRANSQUERYFOFN

if [ $# -eq 0 ]; then echo "
Arg1=scripts dir
Arg2=Config file
Arg3=Logical true/false -- clean as it works?
Arg4=ASMFOFN - file of filenames -- paths to each assembly FASTA - extension for all fasta files should be .fasta
Arg5=TRANSFASTA -- fasta file of transcripts (or other) to align to each assembly
Arg6=NJOBS -- how many jobs the TRANSFASTA should be broken up into.
"; exit; fi

##TRANSQUERYFOFN - file of filenames -- paths to each FASTA to be used as a query file in BLAST.

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
ASMFOFN=$4
TRANSFASTA=$5
NJOBS=$6

PIPELINE=${SCRIPTS}/trans-pipeline.sh

TMP=tmpdir

if [[ "$TRANSFASTA" == *.fa ]]; then PRE=`basename $TRANSFASTA .fa`;
elif [[ "$TRANSFASTA" == *.fasta ]]; then PRE=`basename $TRANSFASTA .fasta`; fi

## BREAK UP QUERY FASTA FILE IN ONE SPOT FOR ALL ASMS TO USE
###if $BREAKUP; then

facount=`grep -c ">" $TRANSFASTA`
count=`echo $facount+$NJOBS | bc`
nlines=`echo $count/$NJOBS | bc`
splitFastA.py -f $TRANSFASTA -n $nlines -o $TMP

QUERYDIR=`readlink -f $TMP`

i=0
while read f; do
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  REF=`readlink -f $f` 
  if [[ "$f" == *.fa ]]; then B=`basename $f .fa`; 
  elif [[ "$f" == *.fasta ]]; then B=`basename $f .fasta`; fi
  echo $B; 
  if [ ! -d $B ]; then mkdir $B; fi
  cd $B;
  $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF $QUERYDIR $PRE $NJOBS
  cd ../
done < $ASMFOFN


