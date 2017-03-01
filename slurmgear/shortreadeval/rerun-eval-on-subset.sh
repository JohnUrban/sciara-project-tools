#!/bin/bash

# specify paths to lap read sample (LR1,LR2) and all reads (R1,R2) -- give dummy answers if will not be using (that will serve as place-holders)
LR1=
LR2=
R1=
R2=
## supply prefixes (prior to .fasta)
SUBSET="platanus.1 platanus.dbg2olc.2 falcon.1 canu.3"



## What programs to use? FILL IN BELOW
ALL=eval.cfg
OnlyAle=eval.aleonly.cfg
OnlyBusco=eval.buscoOnly.cfg
OnlyLap=eval.laponly.cfg
OnlyReapr=eval.reapronly.cfg

## FILL IN WITH CORRECT VARIABLE
EvalThese=$ALL



## May need to adjust the following
ASMDIR=asms
SCRIPTS=/gpfs_home/jurban/software/sciaratools/sciara-project-tools/slurmgear/shortreadeval/scripts/
EVAL=${SCRIPTS}/eval.ARGS.sh
CONFIG=${SCRIPTS}/${EvalThese}



QOS=epscor-condo
for d in $SUBSET; do 
  f=${ASMDIR}/$d.fasta
  ref=`readlink -f $f` 
  echo $d; 
  mkdir $d; 
  cd $d; pwd
  $EVAL $ref $QOS $SCRIPTS $LR1 $LR2 $R1 $R2 $CONFIG
  cd ../
done 
