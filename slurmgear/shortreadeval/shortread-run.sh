#!/bin/bash

## Input Assemblies FOFN
ASMFOFN=input.fofn

# specify paths to lap read sample (LR1,LR2) and all reads (R1,R2)-- give dummy answers if will not be using (that will serve as place-holders)
LR1=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.1.fastq
LR2=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.2.fastq
R1=~/data/scratch/male-ilmn/data/ilmnraw/R1.fastq
R2=~/data/scratch/male-ilmn/data/ilmnraw/R2.fastq


## What programs to use? FILL IN BELOW
ALL=eval.cfg
OnlyAle=eval.aleonly.cfg
OnlyBusco=eval.buscoOnly.cfg
OnlyLap=eval.laponly.cfg
OnlyReapr=eval.reapronly.cfg
OnlyReaprNoClean=eval.reapronly.noclean.cfg
OnlyReaprNoCleanAggressive=eval.reapronly.noclean.aggressive.cfg

## FILL IN WITH CORRECT VARIABLE
EvalThese=$ALL



## May need to adjust the following
SCRIPTS=/gpfs_home/jurban/software/sciaratools/sciara-project-tools/slurmgear/shortreadeval/scripts/
AUTO=${SCRIPTS}/auto-shortreadeval.sh
EVAL=${SCRIPTS}/eval.ARGS.sh
CONFIG=${SCRIPTS}/configs/${EvalThese}


# eval.ARGS.sh
#Arg1=/Path/To/Reference.fasta
#Arg2=QOS
#Arg3=scripts dir path
#Arg4=lap reads 1
#Arg5=lap reads 2
#Arg6=all reads 1
#Arg7=all reads 2
#Arg8=ConfigFile



################ EXECUTE #####################

bash $AUTO $ASMFOFN $LR1 $LR2 $R1 $R2 $EvalThese $SCRIPTS

################ EXECUTE #####################
## Read in options
#FOFN=$1
#LR1=$2
#LR2=$3
#R1=$4
#R2=$5
#EvalThese=$6  
#SCRIPTS=$7
