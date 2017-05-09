#!/bin/bash


rm -r aln/ asm_map/


cd merge/

RM=false
F=allstats.txt
if [ -f $F ]; then X=`cat $F | wc -l`; if [ $X -ge 9 ]; then RM=true; fi; fi

if $RM; then rm *aln *bedGraph *genome; fi

cd ../


