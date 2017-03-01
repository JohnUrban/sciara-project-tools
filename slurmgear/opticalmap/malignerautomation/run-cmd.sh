#!/bin/bash

# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.

CLEAN=false
ASMFOFN=input.fofn
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
