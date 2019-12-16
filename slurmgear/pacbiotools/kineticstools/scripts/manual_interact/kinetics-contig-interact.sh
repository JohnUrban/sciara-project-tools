#!/bin/bash

. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27



## ASSUMES CONFIG FILE PRESENT
CFG=`readlink -f vars.cfg`

## AUTO CONFIG: VARS from CONFIG FILE
TOPDIR=`grep -w ^TOPDIR ${CFG} | awk '{print $2}'`
REF=`grep -w ^REF ${CFG} | awk '{print $2}'`
REFBASE=`grep -w ^REFBASE ${CFG} | awk '{print $2}'`
INPUT=`grep -w ^INPUT ${CFG} | awk '{print $2}'`
NAMES=`grep -w ^NAMES ${CFG} | awk '{print $2}'`
MAX=`grep -w ^MAX ${CFG} | awk '{print $2}'`
PVALUE=`grep -w ^PVALUE ${CFG} | awk '{print $2}'`
IDENTIFY=`grep -w ^IDENTIFY ${CFG} | awk '{print $2}'`
THREADS=1

REFWINDOWS=`grep -w ^REFWINDOWS ${CFG} | awk '{print $2}'`

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

for var in MEM TIME THREADS REFBASE contigName QOS PVALUE IDENTIFY THREADS INPUT REF KC REFWINDOWS; do
  echo ${var} ${!var}
done

date
#ipdSummary ${INPUT} -v -v -v --gff basemods.gff --csv kinetics.csv --reference $REF --identify ${IDENTIFY} --methylFraction \
#	--pvalue ${PVALUE} --minCoverage 3 \
#	--methylMinCov 10 --identifyMinCov 5 -j ${THREADS} --maxAlignments 1000000 \
#	--bigwig ipd.bigWig --refContigs ${contigName}


ipdSummary ${INPUT} -v -v -v --gff basemods.gff --csv kinetics.csv --reference $REF --identify ${IDENTIFY} --methylFraction \
	--pvalue ${PVALUE} --minCoverage 3 \
	--methylMinCov 10 --identifyMinCov 5 -j ${THREADS} --maxAlignments 1000000 \
	--bigwig ipd.bigWig --referenceWindows ${REFWINDOWS}


###--outfile basemods
###--gff basemods.gff --csv kinetics.csv
### --ms_csv multisite.csv
###--m5Cgff m5C.scores.gff
date
