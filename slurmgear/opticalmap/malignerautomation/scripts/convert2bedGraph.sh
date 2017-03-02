#!/bin/bash

#need: ASM,BASE,ALL

faSize -detailed $ASM > ${BASE}.genome
tail -n +2 $ALL | awk 'OFS="\t" { if ($10<0) print $2,0,$11; else print $2,$10,$11}' | sortBed -i - | genomeCoverageBed -i - -g ${BASE}.genome -bg > ${ALL}.bedGraph
awk '{s+=($3-$2)*$4}END{print s}' ${ALL}.bedGraph > total_base_cov.txt

## NOTE: total_base_cov was in error --	($3-$2)*4 instead of ($3-$2)*$4
