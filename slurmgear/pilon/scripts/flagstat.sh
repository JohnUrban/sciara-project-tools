#!/bin/bash


##--export=PRE=reads

samtools flagstat ${PRE}.bam > ${PRE}.flagstats.txt
