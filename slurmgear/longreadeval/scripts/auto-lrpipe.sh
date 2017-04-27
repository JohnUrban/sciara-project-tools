#!/bin/bash

if [ $# -eq 0 ]; then echo "
Arg1=Logical true/false -- clean as it works?
Arg2=Config file
Arg3=FOFN - file of filenames -- paths to each assembly FASTA to polish - extension for all fasta files should be .fasta
Arg4=Scripts dir
Arg5=Ont Fastq Loc
Arg6=PacBio Fastq Loc
Arg7=Ont long2pe R1 Loc
Arg8=Ont long2pe R2 Loc
Arg9=Pb long2pe R1 Loc
Arg10=Pb long2pe R2 Loc

"; exit; fi

CLEAN=$1
CONFIG=$2
FOFN=$3
SCRIPTS=$4
ONT=$5
PACBIO=$6
ONT1=$7
ONT2=$8
PACBIO1=$9
PACBIO2=$10

LRPIPE=${SCRIPTS}/longread-pipeline.sh


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
  $LRPIPE $ref $QOS $CLEAN $CONFIG $SCRIPTS $ONT $PACBIO $ONT1 $ONT2 $PACBIO1 $PACBIO2
  cd ../
done < $FOFN


