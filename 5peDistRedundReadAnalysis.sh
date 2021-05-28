#!/bin/bash

## Preceded by:
## for f in 01-bam/bwa/*bam ; do 
##   BASE=$( basename $f .bam );  
##   bedtools genomecov -5 -ibam ${f} -g ${G} > ${OUTDIR}/${BASE}.txt & 
##   bedtools genomecov -strand "+" -5 -ibam ${f} -g ${G} > ${OUTDIR}/${BASE}.pos.txt & 
##   bedtools genomecov -strand "-" -5 -ibam ${f} -g ${G} > ${OUTDIR}/${BASE}.neg.txt & 
## done


F=$1
K=5
if [ $# -gt 1 ]; then K=$2; fi

B=$( basename $F .txt )
D=$( dirname $F )
PRE=${D}/${B}

for C in ${PRE}.pos.txt ${PRE}.neg.txt; do
  grep genome ${C} | head -n 2 ; 
  grep genome ${C} | awk '$2>1 {s+=$5}END{print "genome\t2+\t-\t-\t"s}'
  grep genome ${C} | awk '$2>1 {s+=$2*$3}END{print "genome\t2+\t-\tNreads\t"s}'
  grep genome ${C} | awk '$2>1 {s+=($2-1)*($3)}END{print "genome\t2+\t-\tNreadsLeave1\t"s}'
  echo
  grep genome ${C} | awk -v "K=$K" '$2>K {s+=$5}END{print "genome\t"K+1"+\t-\t-\t"s}'
  grep genome ${C} | awk -v "K=$K" '$2>K {s+=$2*$3}END{print "genome\t"K+1"+\t-\tNreads\t"s}'
  grep genome ${C} | awk -v "K=$K" '$2>K {s+=($2-K)*($3)}END{print "genome\t"K+1"+\t-\tNreadsLeaveK\t"s}'
  echo
  echo
done
