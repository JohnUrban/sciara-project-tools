#!/bin/bash

export PATH=/users/jurban/software/localpy/pbh5tools/bin/:${PATH}
#module load pbh5tools/1.0

PYTHONPATH=""
. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate genomicconsensus


date

if [ ! -f ${REF}.fai ]; then samtools faidx ${REF}; fi
which quiver
quiver --version
quiver -j${THREADS} ${INPUT} -r ${REF} -o $OUTGFF -o $OUTFASTQ -o $OUTFASTA --noEvidenceConsensusCall=lowercasereference --verbose



date
