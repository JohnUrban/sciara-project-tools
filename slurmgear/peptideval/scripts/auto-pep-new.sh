#!/bin/bash

##$RUN $SCRIPTS $CONFIG $CLEAN $ASMFOFN $TRANSQUERYFOFN

if [ $# -eq 0 ]; then echo "
Arg1=scripts dir
Arg2=Config file
Arg3=Logical true/false -- clean as it works?
Arg4=ASMFOFN - file of filenames -- paths to each assembly FASTA - extension for all fasta files should be .fasta
Arg5=TRANSFASTA -- fasta file of transcripts (or other) to align to each assembly
Arg6=NJOBS -- how many jobs the TRANSFASTA should be broken up into.
Arg7 = Job Prefix
"; exit; fi

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
ASMFOFN=$4
PEPFASTA=$5
NJOBS=$6
JOBPRE=$7

PIPELINE=${SCRIPTS}/peptide-pipeline-new.sh

TMP=tmpdir

if [[ "$PEPFASTA" == *.fa ]]; then PRE=`basename $PEPFASTA .fa`;
elif [[ "$PEPFASTA" == *.fasta ]]; then PRE=`basename $PEPFASTA .fasta`; fi


facount=`grep -c ">" $PEPFASTA`
count=`echo $facount+$NJOBS | bc`
nlines=`echo $count/$NJOBS | bc`
splitFastA.py -f $PEPFASTA -n $nlines -o $TMP

QUERYDIR=`readlink -f $TMP`

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
  $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF $QUERYDIR $PRE $NJOBS ${JOBPRE}
  cd ../
done < $ASMFOFN


