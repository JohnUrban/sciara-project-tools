#!/bin/bash

# copy this script into your working dir
# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.

## FILL THIS IN:
## READSFOFN - tab-sep 2 column -- col1=R1 paths ; col2=R2 paths -- each line has corresponding pairs
READSFOFN=

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

BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/rnaseqeval
SCRIPTS=${BASE}/scripts/
CONFIGS=${BASE}/configs/
FOFNS=${BASE}/fofns/

CONFIG=${CONFIGS}/rnaseq-config-sciara.cfg




RUN=${SCRIPTS}/auto-rnaseqeval.sh


################ EXECUTE #####################

$RUN $SCRIPTS $CONFIG $CLEAN $ASMFOFN $READSFOFN

################ EXECUTE #####################


## This pipeline will be designed to use 1 FOFN file
## It can then be re-run with different specified FOFNs if one wants more... such as from related species.
##TRANSQUERYFOFN=${FOFNS}/sciara.fofn

