#!/bin/bash


ASMFOFN=input.fofn

SFX2RM="" ## Leave "" if none
SFX2ADD="" ## Leave "" if none
CLEAN=true

##READS
R1=/gpfs/data/sgerbi/jurban/illumina/sciara/male_hiseq_130721/MSF0007_TGACCA_L004_R1_001.fastq.gz
R2=/gpfs/data/sgerbi/jurban/illumina/sciara/male_hiseq_130721/MSF0007_TGACCA_L004_R2_001.fastq.gz

## Pilon 
##Fix Options (can give comma-sep list): bases, gaps, local, all, none, amb, breaks (can be used with local), novel
FIX=bases
NOSTRAYS=true

##MAP=reads

## Pipeline Options
BUILDBT2=true
MAPREADS=true
FLAG1=true
FLAG2=true
MARK=true
PILON=true

## RUN INFO LOCATIONS
BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/pilon
SCRIPTS=${BASE}/scripts/
CONFIGS=${BASE}/configs/
AUTO=${SCRIPTS}/auto-pilon.sh
PILONPIPE=${SCRIPTS}/pilon-pipeline.sh


## Diff cfg setups
SCIARA=${CONFIGS}/pilon-config-sciara.cfg
##SKIP2MK=${CONFIGS}/pilon-config-sciara-skip-to-mkdup.cfg
##SKIP2MAP=${CONFIGS}/pilon-config-sciara-skip-to-map.cfg

## config -- essentially has time, mem, cpu instructions for slurm
CONFIG=${SCIARA}



################ EXECUTE #####################

bash $AUTO $SFX2RM $SFX2ADD $CLEAN $CONFIG $ASMFOFN $BUILDBT2 $MAPREADS $FLAG1 $FLAG2 $MARK $PILON $FIX $NOSTRAYS $PILONPIPE $R1 $R2

################ EXECUTE #####################

