#!/usr/bin/env bash

###################################################################

BASE=`basename $ASM_FASTA .fasta`
BASE2=`basename $SMOOTH_MAPS`


# outputs
RMAP_OUT_PFX=${BASE2}
OUT_PFX=${BASE2}

for var in REC_ENZ REC_SEQ MIN_FRAG_SIZE QUERY_MISS_PENALTY REF_MISS_PENALTY QUERY_MAX_MISSES REF_MAX_MISSES SD_RATE MIN_SD MAX_SCORE_PER_INNER_CHUNK MAX_ALIGNMENTS_PER_QUERY BASE BASE2 RMAP_OUT_PFX OUT_PFX; do
  echo ${var} ${!var}
done

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

