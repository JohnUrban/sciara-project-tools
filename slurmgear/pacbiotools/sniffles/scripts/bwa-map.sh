#!/bin/bash
source ~/software/bwa/source.sh
####bwa mem -M -x $TYPE $BWAIDXPRE $FASTQ | samtools sort --threads > $TYPE.bam
bwa mem -t $MTHREADS -M -x $TYPE $BWAIDX $FASTQ | samtools sort -T $TYPE --threads $MTHREADS -o $TYPE.bam
