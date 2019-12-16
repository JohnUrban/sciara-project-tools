#!/bin/bash
###########################
## 1. Make a file called input.fofn with path to all .bax.h5 files
## 2. Fill in REF and ASMOUTPREFIX below; change INPUT if necessary
## 3. Make changes to other parameters if necessary
## 4. Run as  "bash quiverpipeline.parallel.sh"

REF=                ## path to pre-quiver asm to polish
ASMOUTPREFIX=       ## e.g. falcon.pball
INPUT=/gpfs_home/jurban/software/quiver/scripts/input.fofn #input.fofn

module load blasr/2015Oct22-8cc8621

PATH_TO_SCRIPTS=/gpfs_home/jurban/software/quiver/scripts
QUIVER=true

MERGEDCMP=out_all.cmp.h5
MAIN_TMP_DIR=temp
TMP_SORT=${MAIN_TMP_DIR}/cmp_sort
QTHREADS=8
QMEM=100g
OUTGFF=${ASMOUTPREFIX}.quiver.variants.gff
OUTFASTQ=${ASMOUTPREFIX}.quiver.fastq
OUTFASTA=${ASMOUTPREFIX}.quiver.fasta
TIME=48:00:00
QOS=biomed-sb-condo
ALTQOS=${QOS} #ccmb-condo
CMPDIR=cmpfiles
SLURMOUTDIR=slurmout


##############################################################################
## QUIVER
##############################################################################

sbatch -J ${ASMOUTPREFIX}_quiver -o ${SLURMOUTDIR}/quiver.slurm.%A.out --mem=$QMEM --time=$TIME -c $QTHREADS --qos=$QOS --export=THREADS=${QTHREADS},INPUT=${MERGEDCMP},REF=${REF},OUTGFF=${OUTGFF},OUTFASTQ=${OUTFASTQ},OUTFASTA=${OUTFASTA} ${PATH_TO_SCRIPTS}/quiver.sh | awk '{print $4}'
