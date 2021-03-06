#!/bin/bash

## ASSUMES DOING THIS IN TOPDIR ABOVE REF BASE_NAME DIRS 
##	AS OPPOSED TO ASSUMING ALREADY IN REF_BASENAME_DIR AND THAT THIS IS LAST STEP IN KINETICSPIPELINE 
## MAKES/ENTERS REF_BASENAME_DIR
## 	B/C OF THIS, IT ASSUMES THE FIRST 3 STEPS OF KINETICS PIPELINE ALREADY DONE (PBALIGN, MERGE, SORT) AND EXIST IN REF_BASENAME_DIR
##	ALSO B/C OF THIS, IT CAN BE RUN IN THE TOPDIR IN FOR LOOP OVER REFs
## MAKES IPDSUMMARY SUBDIR, ENTERS, ISSUES 1 JOB PER CONTIG

. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27

## ASSUMES CONFIG FILE PRESENT
CFG=`readlink -f ../vars.cfg`

## AUTO CONFIG:	VARS from CONFIG FILE
QOS=`grep -w ^QOS ${CFG} | awk '{print $2}'`
TOPDIR=`grep -w ^TOPDIR ${CFG} | awk '{print $2}'`
REF=`grep -w ^REF ${CFG} | awk '{print $2}'`
REFBASE=`grep -w ^REFBASE ${CFG} | awk '{print $2}'`
INPUT=`grep -w ^INPUT ${CFG} | awk '{print $2}'`
NAMES=`grep -w ^NAMES ${CFG} | awk '{print $2}'`
MAX=`grep -w ^MAX ${CFG} | awk '{print $2}'`
PVALUE=`grep -w ^PVALUE ${CFG} | awk '{print $2}'`
IDENTIFY=`grep -w ^IDENTIFY ${CFG} | awk '{print $2}'`
THREADS=`grep -w ^THREADS ${CFG} | awk '{print $2}'`
TIME=`grep -w ^TIME ${CFG} | awk '{print $2}'`
MEM=`grep -w ^MEM ${CFG} | awk '{print $2}'`

## AUTO CONFIG: detecting where are in DIR hierarchy
contig=$(basename `readlink -f ../`)
debug=$(basename `readlink -f ./`)
KC=${TOPDIR}/${REFBASE}/debug/${contig}/${debug}/kinetics-contig.sh



## LAUNCH SLURM JOB
sbatch --mem=${MEM} --time=${TIME} -c ${THREADS} -J debugging_${REFBASE}_${contig}_ipd --qos=${QOS} \
      --export=PVALUE=${PVALUE},IDENTIFY=${IDENTIFY},THREADS=${THREADS},INPUT=${INPUT},REF=${REF},contigName=${contig} ${KC}

