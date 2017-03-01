#!/bin/bash

# specify paths to lap read sample (LR1,LR2) and all reads (R1,R2)
LR1=
LR2=
R1=
R2=



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


# eval.ARGS.sh
#Arg1=/Path/To/Reference.fasta
#Arg2=QOS
#Arg3=scripts dir path
#Arg4=lap reads 1
#Arg5=lap reads 2
#Arg6=all reads 1
#Arg7=all reads 2
#Arg8=ConfigFile


i=0

for f in ${ASMDIR}/*fasta; do 
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  ref=`readlink -f $f` 
  b=`basename $f .fasta`; 
  echo $b; 
  mkdir $b; 
  cd $b;
  $EVAL $ref $QOS $SCRIPTS $LR1 $LR2 $R1 $R2 $CONFIG
  cd ../
done 

