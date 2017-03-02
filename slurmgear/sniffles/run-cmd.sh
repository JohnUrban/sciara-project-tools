#!/bin/bash

# copy this script into your working dir
# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.

## ASM FOFN
ASMFOFN=input.fofn

## LONG READ LOCATIONS
ONT=~/data/scratch/minion2016/fast5fastqs/allReadsFromAllONTlibsCombined.fastq
PACBIO=~/data/scratch/pac_bio_data/filt/all_subreads.fastq

## RUN INFO LOCATIONS
BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/sniffles
SCRIPTS=${BASE}/scripts/
AUTOSNIFF=${SCRIPTS}/auto-sniff.sh
CONFIGS=${BASE}/configs/
CONFIG=${CONFIGS}/sniffles-config-sciara.cfg

## OTHER OPTIONS
CLEAN=false


################ EXECUTE #####################

bash $AUTOSNIFF $CLEAN $CONFIG $ASMFOFN $SCRIPTS $ONT $PACBIO

################ EXECUTE #####################
