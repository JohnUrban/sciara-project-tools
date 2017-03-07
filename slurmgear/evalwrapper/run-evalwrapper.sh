#!/bin/bash

## Input Assemblies FOFN
ASMFOFN=input.fofn

###################  SHORT READ  ###########################
# specify paths to lap read sample (LR1,LR2) and all reads (R1,R2)-- give dummy answers if will not be using (that will serve as place-holders)
LR1=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.1.fastq
LR2=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.2.fastq
R1=~/data/scratch/male-ilmn/data/ilmnraw/R1.fastq
R2=~/data/scratch/male-ilmn/data/ilmnraw/R2.fastq

## OPTIONS for what programs to use. 
## FILL IN   ==>"EvalThese"<==   BELOW WITH ONE OF THESE
ALL=eval.cfg
OnlyAle=eval.aleonly.cfg
OnlyBusco=eval.buscoOnly.cfg
OnlyLap=eval.laponly.cfg
OnlyReapr=eval.reapronly.cfg
OnlyReaprNoClean=eval.reapronly.noclean.cfg
OnlyReaprNoCleanAggressive=eval.reapronly.noclean.aggressive.cfg

## FILL IN WITH CORRECT VARIABLE
EvalThese=$ALL

## May need to adjust the following
SHORTSCRIPTS=/gpfs_home/jurban/software/sciaratools/sciara-project-tools/slurmgear/shortreadeval/scripts/
SHORTAUTO=${SHORTSCRIPTS}/auto-shortreadeval.sh
SHORTEVAL=${SHORTSCRIPTS}/eval.ARGS.sh
SHORTCONFIG=${SHORTSCRIPTS}/configs/${EvalThese}


############### BIONANO MALIGNER SECTION #################
BIONANOCLEAN=false
REC_ENZ=BssSI
REC_SEQ=CACGAG

BIONANOBASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/opticalmap/malignerautomation
BIONANOSCRIPTS=${BIONANOBASE}/scripts/
BIONANOCONFIGS=${BIONANOBASE}/configs/
BIONANOFOFNS=${BIONANOBASE}/fofns/

BIONANOCONFIG=${BIONANOCONFIGS}/maligner-config-sciara.cfg
MAPSFOFN=${BIONANOFOFNS}/bionanomaps.examp2.fofn
BIONANORUN=${BIONANOSCRIPTS}/auto-malign.sh




############### LONG READ SNIFFLES SECTION #################
## LONG READ LOCATIONS
ONT=~/data/scratch/minion2016/fast5fastqs/allReadsFromAllONTlibsCombined.fastq
PACBIO=~/data/scratch/pac_bio_data/filt/all_subreads.fastq

## RUN INFO LOCATIONS
SNIFFBASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/sniffles
SNIFFSCRIPTS=${SNIFFBASE}/scripts/
AUTOSNIFF=${SNIFFSCRIPTS}/auto-sniff.sh
SNIFFCONFIGS=${SNIFFBASE}/configs/
SNIFFCONFIG=${SNIFFCONFIGS}/sniffles-config-sciara.cfg

## OTHER OPTIONS
SNIFFCLEAN=false





##############################################
##############################################
##############################################
################ EXECUTE #####################

mkdir shortread
cd shortread
bash $SHORTAUTO $ASMFOFN $LR1 $LR2 $R1 $R2 $EvalThese $SHORTSCRIPTS
cd ../

mkdir bionano
cd bionano
$BIONANORUN $BIONANOCLEAN $BIONANOCONFIG $ASMFOFN $MAPSFOFN $REC_ENZ $REC_SEQ $BIONANOSCRIPTS
cd ../

mkdir longread
cd longread 
bash $AUTOSNIFF $SNIFFCLEAN $SNIFFCONFIG $ASMFOFN $SNIFFSCRIPTS $ONT $PACBIO
cd ../

################ EXECUTE #####################
##############################################
##############################################
##############################################





