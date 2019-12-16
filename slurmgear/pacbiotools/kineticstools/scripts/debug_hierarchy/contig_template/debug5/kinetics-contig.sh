#!/bin/bash

. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27

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

echo PVALUE ${PVALUE}
echo IDENTIFY ${IDENTIFY}
echo THREADS ${THREADS}
echo INPUT ${INPUT}
echo REF ${REF}
echo contigName $contigName

echo "
ipdSummary ${INPUT} --reference $REF --identify ${IDENTIFY} --methylFraction \
	--gff basemods.gff --csv kinetics.csv --pvalue ${PVALUE} --minCoverage 3 \
	--methylMinCov 10 --identifyMinCov 5 -j ${THREADS} --maxAlignments 1000000 \
	--ms_csv multisite.csv --bigwig ipd.bigWig --refContigs $contigName
"

date
ipdSummary ${INPUT} --reference $REF  \
	--gff basemods.gff --csv kinetics.csv --refContigs ${contigName}
date
