#!/bin/bash

# copy this script into your working dir
# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.

## FILL THIS IN:
PEPFASTA=
NJOBS=100
## Also do TBLASTX in addition to BLASTN ???

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

BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/peptideval
SCRIPTS=${BASE}/scripts/
CONFIGS=${BASE}/configs/
FOFNS=${BASE}/fofns/

CONFIG=${CONFIGS}/peptide-config-sciara.cfg




RUN=${SCRIPTS}/auto-pep.sh


################ EXECUTE #####################

$RUN $SCRIPTS $CONFIG $CLEAN $ASMFOFN $PEPFASTA $NJOBS 

################ EXECUTE #####################


