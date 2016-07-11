#!/bin/bash
BT2=$1
P1=$2
P2=$3

if [ $# -eq 0 ]; then
  echo 
  echo Usage: 
  echo get_PCR_product_from_primers.sh    bt2index    primer1.fa   primer2.fa    genome_sequence.fa
  echo
  exit 0
fi



DIR=`dirname $2`
BASE=`basename $2 | stringEdit - .fa`
TARGET=${DIR}/${BASE}-target.fa
bowtie2 --no-discordant --no-mixed -x $1 -1 $2 -2 $3 -f -a --maxins 10000 | samtools view -S -F 4 - | grep -v ^@ | grep ^p1 | awk 'OFS="\t" {print $3,$4-1,$4-1+$9}' 
