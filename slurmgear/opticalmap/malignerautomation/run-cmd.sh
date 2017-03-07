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



CLEAN=false
REC_ENZ=BssSI
REC_SEQ=CACGAG

BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/opticalmap/malignerautomation
SCRIPTS=${BASE}/scripts/
CONFIGS=${BASE}/configs/
FOFNS=${BASE}/fofns/

CONFIG=${CONFIGS}/maligner-config-sciara.cfg
MAPSFOFN=${FOFNS}/bionanomaps.examp2.fofn
RUN=${SCRIPTS}/auto-malign.sh


################ EXECUTE #####################

$RUN $CLEAN $CONFIG $ASMFOFN $MAPSFOFN $REC_ENZ $REC_SEQ $SCRIPTS

################ EXECUTE #####################


##CONFIG=/gpfs/data/sgerbi/jurban/software/maligner/scripts/maligner-config-sciara.cfg
##MAPSFOFN=/gpfs/data/sgerbi/jurban/software/maligner/scripts/bionanomaps.fofn
##RUN=/gpfs/data/sgerbi/jurban/software/maligner/scripts/auto-malign.sh 
##SCRIPTS=/users/jurban/data/software/maligner/scripts
