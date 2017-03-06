#!/bin/bash

##--export=JX=${JX},JAR=$MKDUPJAR,PRE=reads
## OLD SYNTAX: java -Xmx${JX} -jar $JAR INPUT=${PRE}.bam OUTPUT=${PRE}.markdup.bam METRICS_FILE=${PRE}.metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true

## New Syntax:
java -Xmx${JX} -jar $JAR MarkDuplicates INPUT=${PRE}.bam OUTPUT=${PRE}.markdup.bam METRICS_FILE=${PRE}.metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true
samtools index ${PRE}.markdup.bam
