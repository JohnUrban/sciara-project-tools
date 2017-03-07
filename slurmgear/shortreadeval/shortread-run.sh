#!/bin/bash

function help {
    echo "
    Usage1: bash $0 FOFN
    Usage2: bash $0  ((FOFN filled in manually inside script -- defaults to input.fofn))

    ...where FOFN has list of all assemblies used in the assembly evaluations in subdirs you are trying to evaluate and summarize.
    (( typically called input.fofn ))
    "
}

## NEED HELP?
if [ $# -eq 1 ]; then
    if [ $1 == "-h" ] || [ $1 == "-help" ] || [ $1 == "--help" ]; then
        help; exit
    fi
fi

## DEFAULT ASMFOFN
ASMFOFN=input.fofn

## OPTIONAL ASMFOFN
if [ $# -eq 1 ]; then ASMFOFN=$1; fi

## IF NEITHER ASMFOFN WORKED -- THEN REPORT ERROR AND EXIT
if [ ! -f $ASMFOFN ]; then echo; echo "    ASMFOFN ERROR: FILE NOT FOUND"; echo "    YOU GAVE:"; echo "    $ASMFOFN"; help; exit; fi


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
