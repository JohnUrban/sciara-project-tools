#!/bin/bash
CLEAN=false
CONFIG=/gpfs/data/sgerbi/jurban/software/maligner/scripts/maligner-config-sciara.cfg
ASMFOFN=input.fofn
MAPSFOFN=/gpfs/data/sgerbi/jurban/software/maligner/scripts/bionanomaps.fofn
REC_ENZ=BssSI
REC_SEQ=CACGAG

##RUN=/gpfs/data/sgerbi/jurban/software/maligner/scripts/auto-malign.sh 
RUN=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/opticalmap/malignerautomation/auto-malign.sh

$RUN $CLEAN $CONFIG $ASMFOFN $MAPSFOFN $REC_ENZ $REC_SEQ

