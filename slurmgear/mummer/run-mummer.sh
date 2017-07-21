#!/bin/bash

mkdir nuc
mkdir delfil
mkdir coords

echo NUCMER
##default breaklen = 200
QUERY=${READSDIR}/${PRE}.${SLURM_ARRAY_TASK_ID}.fa
${MUMMER}/nucmer -minmatch ${minMatchLength} -prefix nuc/${PREFIX}.${SLURM_ARRAY_TASK_ID} -breaklen ${breakLen} -mincluster ${mincluster} -maxgap ${maxgap} ${REF} ${QUERY}

echo; echo; echo

echo DELTA FILTER
${MUMMER}/delta-filter -i ${identity} -l ${minAlnLen} -q nuc/${PREFIX}.${SLURM_ARRAY_TASK_ID}.delta > delfil/${PREFIX}.${SLURM_ARRAY_TASK_ID}.filtered.delta

echo; echo; echo

echo SHOW COORD
${MUMMER}/show-coords -c -H -T -d -l -q delfil/${PREFIX}.${SLURM_ARRAY_TASK_ID}.filtered.delta | awk 'OFS="\t" {print $14,$1-1,$2,$7,$15,$3-1,$4,$5,$6,$8,$9,$10,$11,$12,$13}' > coords/${PREFIX}.${SLURM_ARRAY_TASK_ID}.filtered.delta.bed

