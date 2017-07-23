#!/usr/bin/env bash

export PATH=~/data/software/maligner/maligner/bin/:$PATH
export PYTHONPATH=~/data/software/maligner/maligner/lib/:$PYTHONPATH
export PATH=~/data/software/maligner/maligner/build/bin/:$PATH


##Need: ASM_FASTA, BASE, REC_ENZ, REC_SEQ, MIN_FRAG_SIZE, FASTAFOFN, MAPSFOFN

# outputs


# convert the asm fasta file to the Maligner maps format and smooth the maps file by merging consecutive fragments that are less than 1kb
i=0
while read fastaloc; do
  let i++
  if [[ "$fastaloc" == *.fasta ]]; then BASE=`basename ${fastaloc} .fasta`; fi
  if [[ "$fastaloc" == *.fa ]]; then BASE=`basename ${fastaloc} .fa`; fi
  OUT_PFX=fastaloc_${i}.${BASE}.${REC_ENZ}
  make_insilico_map -o $OUT_PFX $fastaloc $REC_SEQ
  smooth_maps_file -m $MIN_FRAG_SIZE ${OUT_PFX}.maps > ${OUT_PFX}.smoothed.maps
done < $FASTAFOFN

## APPEND locations of newly generated smoothed map files to maps.fofn
for smoothmapfile in *.smoothed.maps; do readlink -f $smoothmapfile ; done >> ${MAPSFOFN}

