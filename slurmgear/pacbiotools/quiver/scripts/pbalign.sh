#!/bin/bash

module load blasr/2015Oct22-8cc8621
export PATH=/users/jurban/software/localpy/pbh5tools/bin/:${PATH}
. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27

#module load pbh5tools/1.0

date
pbalign $BAX $REF $OUTPRE.cmp.h5 --forQuiver --tmpDir $TMPDIR --nproc $THREADS
date
