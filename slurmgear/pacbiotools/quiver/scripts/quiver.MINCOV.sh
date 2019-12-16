#!/bin/bash


date

if [ ! -f ${REF}.fai ]; then samtools faidx ${REF}; fi
quiver --minCoverage ${MINCOV} -j${THREADS} ${INPUT} -r ${REF} -o $OUTGFF -o $OUTFASTQ -o $OUTFASTA --noEvidenceConsensusCall=lowercasereference --verbose



date
