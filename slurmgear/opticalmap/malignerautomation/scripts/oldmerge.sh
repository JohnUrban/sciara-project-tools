#!/bin/bash

head -n 1 ${ALN}/bionano.smoothed.maps00.aln > all.bionano.smoothed.maps.aln
for f in ${ALN}/bionano.smoothed.maps*aln; do
    tail -n +2 $f >> all.bionano.smoothed.maps.aln
done
tail -n +2 all.bionano.smoothed.maps.aln | wc -l > num_alignments.txt
