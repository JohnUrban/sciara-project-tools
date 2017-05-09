#!/bin/bash

module load boost/1.55.0  ##for frc
PATH=~/software/frcbam/FRC_align/bin/:/users/jurban/software/sciaratools/sciara-project-tools/:$PATH

echo BAM ${BAM}
echo L2PE_MAXINS ${L2PE_MAXINS}
echo GSIZE ${GSIZE}
echo BASE ${BASE}

FRC --pe-sam $BAM --pe-max-insert $L2PE_MAXINS --genome-size $GSIZE --output ${BASE}.long2pe

if $CLEAN; then bash $SCRIPTS/frc.clean.sh ; fi
