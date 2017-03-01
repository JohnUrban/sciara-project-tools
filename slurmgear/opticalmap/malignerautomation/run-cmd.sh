#!/bin/bash

# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.

CLEAN=false
CONFIG=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/opticalmap/malignerautomation/maligner-config-sciara.cfg
ASMFOFN=input.fofn
MAPSFOFN=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/opticalmap/malignerautomation/bionanomaps.examp2.fofn
REC_ENZ=BssSI
REC_SEQ=CACGAG

RUN=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/opticalmap/malignerautomation/auto-malign.sh



################ EXECUTE #####################

$RUN $CLEAN $CONFIG $ASMFOFN $MAPSFOFN $REC_ENZ $REC_SEQ

################ EXECUTE #####################


##CONFIG=/gpfs/data/sgerbi/jurban/software/maligner/scripts/maligner-config-sciara.cfg
##MAPSFOFN=/gpfs/data/sgerbi/jurban/software/maligner/scripts/bionanomaps.fofn
##RUN=/gpfs/data/sgerbi/jurban/software/maligner/scripts/auto-malign.sh 
