#!/usr/bin/env bash

export PATH=~/data/software/maligner/maligner/bin/:$PATH
export PYTHONPATH=~/data/software/maligner/maligner/lib/:$PYTHONPATH
export PATH=~/data/software/maligner/maligner/build/bin/:$PATH

###################################################################

BASE=`basename $ASM_FASTA .fasta`
BASE2=`basename $SMOOTH_MAPS`


# outputs
RMAP_OUT_PFX=${BASE2}
OUT_PFX=${BASE2}

# setttings
###REC_SEQ=CACGAG # recognition sequence for BssSI == define in above env

MIN_FRAG_SIZE=1000 # bp units
QUERY_MISS_PENALTY=3.0
REF_MISS_PENALTY=3.0
QUERY_MAX_MISSES=5
REF_MAX_MISSES=5
SD_RATE=0.05
MIN_SD=750
MAX_SCORE_PER_INNER_CHUNK=1.0
##MAX_ALIGNMENTS_PER_QUERY=5
MAX_ALIGNMENTS_PER_QUERY=1


###################################################################
# RUN MALIGNER_DP

# Align the smoothed query rmaps file to the smoothed contig maps file with maligner_dp.
# 
# The alignments are saved to $OUT_PFX.aln,
# Log file written to $OUT_PFX.
echo run
maligner_dp \
  -q $QUERY_MISS_PENALTY \
  -r $REF_MISS_PENALTY \
  --query-max-misses $QUERY_MAX_MISSES \
  --ref-max-misses $REF_MAX_MISSES \
  --max-score-per-inner-chunk $MAX_SCORE_PER_INNER_CHUNK \
  --sd-rate $SD_RATE \
  --min-sd $MIN_SD \
  --max-alignments $MAX_ALIGNMENTS_PER_QUERY \
  ${SMOOTH_MAPS} \
  ${ASM_MAP} \
  2>&1 1> ${OUT_PFX}.aln | tee ${OUT_PFX}.log


##  ${RMAP_OUT_PFX}.smoothed.maps \

