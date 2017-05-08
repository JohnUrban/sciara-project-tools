#!/bin/bash

#logicals
# MKDUPS, RUNPILON
# NOSTRAYS, CHANGES, VCF, TRACKS

#non-logicals 
# JX, PICARDJAR, BAM, PILONJAR, ASM, FIX

VARS="MKDUPS RUNPILON NOSTRAYS CHANGES VCF TRACKS JX PICARDJAR BAM PILONJAR ASM FIX"
for var in $VARS; do echo $VAR ${!VAR}; done; echo

if $MKDUPS; then
 java -Xmx${JX} -jar $PICARDJAR MarkDuplicates INPUT=${BAM} OUTPUT=markdup.bam METRICS_FILE=markdup.metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true 2> mkdups.err
 BAM=markdup.bam
 samtools index $BAM
fi

if $RUNPILON; then
 VCF=false
 if $NOSTRAYS; then nostrays=--nostrays; fi
 if $CHANGES; then changes=--changes; fi
 if $VCF; then vcf=--vcf; fi
 if $TRACKS; then tracks=--tracks; fi
 echo "java -Xmx${JX} -jar $PILONJAR --genome $ASM --output pilon --frags ${BAM} --diploid --fix $FIX $nostrays $changes $vcf $tracks"
 java -Xmx${JX} -jar $PILONJAR --genome $ASM --output pilon --frags ${BAM} --diploid --fix $FIX $nostrays $changes $vcf $tracks 1> pilon.err
 echo confirmed total pct_confirmed | awk 'OFS="\t" {print $1,$2,$3}' > confirmed.txt
 grep Confirmed pilon.err | awk '{s+=$2; t+=$4}END{print s, t, 100.0*s/t}' >> confirmed.txt
fi


if $CLEAN; then
#
 if $MKDUPS; then
  rm markdup.bam
 fi
fi

