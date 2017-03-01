#!/bin/bash


ASMDIR=asms

SUBSET="canu.corcov500minrl500.aspb-e025.pball.ontmol.ALTquiver6x" 

EVAL=/gpfs_home/jurban/software/sciaratools/sciara-project-tools/slurmgear/shortreadeval/scripts/eval.ARGS.sh


QOS=epscor-condo
for d in $SUBSET; do 
  f=${ASMDIR}/$d.fasta
  ref=`readlink -f $f` 
  echo $d; 
  mkdir $d; 
  cd $d; pwd
  $EVAL $ref $QOS 
  cd ../
done 

pwd
