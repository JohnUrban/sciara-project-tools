#!/bin/bash

module load boost/1.55.0  ##for frc
PATH=~/software/frcbam/FRC_align/bin/:/users/jurban/software/sciaratools/sciara-project-tools/:$PATH

FRC --pe-sam $BAM --pe-max-insert 800 --genome-size 292000000 --output ${BASE}.frc
