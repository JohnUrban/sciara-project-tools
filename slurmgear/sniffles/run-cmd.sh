#!/bin/bash

# copy this script into your working dir
# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.


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




## LONG READ LOCATIONS
ONT=~/data/scratch/minion2016/fast5fastqs/allReadsFromAllONTlibsCombined.fastq
PACBIO=~/data/scratch/pac_bio_data/filt/all_subreads.fastq

## RUN INFO LOCATIONS
BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/sniffles
SCRIPTS=${BASE}/scripts/
AUTOSNIFF=${SCRIPTS}/auto-sniff.sh
CONFIGS=${BASE}/configs/
CONFIG=${CONFIGS}/sniffles-config-sciara.cfg

## OTHER OPTIONS
CLEAN=false


################ EXECUTE #####################

bash $AUTOSNIFF $CLEAN $CONFIG $ASMFOFN $SCRIPTS $ONT $PACBIO

################ EXECUTE #####################
