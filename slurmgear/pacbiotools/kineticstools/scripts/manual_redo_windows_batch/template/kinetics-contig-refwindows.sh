#!/bin/bash

. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27



## ASSUMES CONFIG FILE PRESENT
CFG=`readlink -f vars.cfg`

## AUTO CONFIG: VARS from CONFIG FILE
REF=`grep -w ^REF ${CFG} | awk '{print $2}'`
REFBASE=`grep -w ^REFBASE ${CFG} | awk '{print $2}'`
REFWINDOWS=`grep -w ^REFWINDOWS ${CFG} | awk '{print $2}'`
PROBLEM_REFWINDOWS=`grep -w ^PROBLEM_REFWINDOWS ${CFG} | awk '{print $2}'`
QOS=`grep -w ^QOS ${CFG} | awk '{print $2}'`
PVALUE=`grep -w ^PVALUE ${CFG} | awk '{print $2}'`
IDENTIFY=`grep -w ^IDENTIFY ${CFG} | awk '{print $2}'`
THREADS=`grep -w ^THREADS ${CFG} | awk '{print $2}'`
TIME=`grep -w ^TIME ${CFG} | awk '{print $2}'`
MEM=`grep -w ^MEM ${CFG} | awk '{print $2}'`
TOPDIR=`grep -w ^TOPDIR ${CFG} | awk '{print $2}'`
INPUT=`grep -w ^INPUT ${CFG} | awk '{print $2}'`
KC=`grep -w ^KC ${CFG} | awk '{print $2}'`


if [ ! -f ${REF}.fai ]; then samtools faidx ${REF}; fi

#-- specify integer 1-7; 1=m6A, 2=m4C, 3=m5C_TET, 4=m6A,m4C, 5=m6A,m5C_TET, 6=m4C,m5C_TET, 7=m6A,m4C,m5C_TET
if [ ${IDENTIFY} -eq 1 ]; then IDENTIFY=m6A ;
elif [ ${IDENTIFY} -eq 2 ]; then IDENTIFY=m4C ;
elif [ ${IDENTIFY} -eq 3 ]; then IDENTIFY=m5C_TET ;
elif [ ${IDENTIFY} -eq 4 ]; then IDENTIFY=m6A,m4C ;
elif [ ${IDENTIFY} -eq 5 ]; then IDENTIFY=m6A,m5C_TET ;
elif [ ${IDENTIFY} -eq 6 ]; then IDENTIFY=m4C,m5C_TET ;
elif [ ${IDENTIFY} -eq 7 ]; then IDENTIFY=m6A,m4C,m5C_TET ;
else echo IDENTIFY NEEDS TO BE 1,2,3,4,5,6,7; exit
fi


## REPORT VARS
for var in REF REFBASE REFWINDOWS PROBLEM_REFWINDOWS QOS PVALUE IDENTIFY THREADS TIME MEM TOPDIR INPUT KC; do
  echo ${var} ${!var}
done

## MAKE CMDs BASED ON ARGS/VARS
ID_CMD="ipdSummary ${INPUT} -v -v -v --identify ${IDENTIFY} --methylFraction \
	--gff basemods-identified.gff --csv kinetics-identified.csv --bigwig ipd-identified.bigWig \
	--reference ${REF} --pvalue ${PVALUE} --minCoverage 3 --methylMinCov 10 --identifyMinCov 5 -j ${THREADS} --maxAlignments 1000000 \
	--referenceWindows ${REFWINDOWS}" 

NOID_CMD="ipdSummary ${INPUT} -v -v -v \
	--gff basemods-unidentified.gff --csv kinetics-unidentified.csv --bigwig ipd-unidentified.bigWig \
	--reference ${REF} --pvalue ${PVALUE} --minCoverage 3 --methylMinCov 10 --identifyMinCov 5 -j ${THREADS} --maxAlignments 1000000 \
	--referenceWindows ${PROBLEM_REFWINDOWS}" 



## EXECUTE CMDs IF WINDOWS VARS ARE NON-EMPTY
date
if [ ! -z ${PROBLEM_REFWINDOWS} ]; then ${NOID_CMD} > unidentified.err 2>&1 & fi

if [ ! -z ${REFWINDOWS} ]; then ${ID_CMD} > identified.err 2>&1 & fi

## WAIT FOR BOTH BG JOBS TO BE DONE
wait
date
