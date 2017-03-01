#!/bin/bash

# copy this script into your working dir
# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.


BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/sniffles
SCRIPTS=${BASE}/scripts/
CONFIGS=${BASE}/configs/
AUTOSNIFF=${SCRIPTS}/auto-sniff.sh


ASMFOFN=input.fofn
CLEAN=false
CONFIG=${CONFIGS}/sniffles-config-sciara.cfg

################ EXECUTE #####################

bash $AUTOSNIFF $CLEAN $CONFIG $ASMFOFN $SCRIPTS

################ EXECUTE #####################
