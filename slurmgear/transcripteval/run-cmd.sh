#!/bin/bash

# copy this script into your working dir
# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.

## FILL THIS IN:
TRANSFASTA=
NJOBS=100
PRE="" ## When running multiple blast analyses - this will add prefix to job name to help distinguish. Need to add underscore if you want it.
## Also do TBLASTX in addition to BLASTN ???
TBLASTX=true

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

BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/transcripteval
SCRIPTS=${BASE}/scripts/
CONFIGS=${BASE}/configs/
FOFNS=${BASE}/fofns/

CONFIG=${CONFIGS}/trans-config-sciara.cfg




RUN=${SCRIPTS}/auto-trans.sh


################ EXECUTE #####################

$RUN $SCRIPTS $CONFIG $CLEAN $ASMFOFN $TRANSFASTA $NJOBS $TBLASTX $PRE

################ EXECUTE #####################


## This pipeline will be designed to use 1 FOFN file
## It can then be re-run with different specified FOFNs if one wants more... such as from related species.
##TRANSQUERYFOFN=${FOFNS}/sciara.fofn

