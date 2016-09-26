#!/bin/bash

function help { echo; echo Usage:; echo digestMultiFasta multi.fa comma-sep-restriction-seqs; echo ; 
  echo This wraps around multiFasta2manyFastas.py and restrictionDigest.py.;
  echo It is a band-aid to deal with multi-fastas.;
  echo TODO: modify restrictionDigest.py to deal directly with it;
  echo;
}

if [ $# -eq 0 ]; then help; exit; fi

FA=$1
R=$2
TEMP=multidigest_$RANDOM

mkdir $TEMP
cd $TEMP
multiFasta2manyFastas.py -f ../$FA
cd ../


for fa in $TEMP/*.fa; do
 b=`basename $fa .fa`
 echo ${b}: $R
 restrictionDigest.py -f $fa -r $R -T -O > $TEMP/${b}.txt
done

