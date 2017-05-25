#!/bin/bash




cd merge/

RM=false
F=allstats.txt
if [ -f $F ]; then X=`awk '$2!=""' $F | wc -l`; if [ $X -ge 9 ]; then RM=true; fi; fi

if $RM; then 
  rm *aln *bedGraph *genome; 
  rm -r ../aln/ ../asm_map/
fi

cd ../


