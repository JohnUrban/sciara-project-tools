#!/bin/bash
BT2=$1
P1=$2
P2=$3

if [ $# -eq 0 ]; then
  echo 
  echo Usage: 
  echo productSpecificity.sh   bt2index    primer1.fa   primer2.fa    blastn_index    genome_sequence.fa
  echo
  exit 0
fi



echo OBTAINING PCR PRODUCT FROM PRIMER MAPPING
get_PCR_product_from_primers.sh $1 $2 $3 $5
echo

echo BLASTING PCR TARGET

DIR=`dirname $2`
BASE=`basename $2 | stringEdit - .fa`
TARGET=${DIR}/${BASE}-target.fa
BLAST=`blastn -db $4 -query $TARGET -outfmt 6 -culling_limit 1 | tee /dev/stderr | wc -l`
echo

echo "PCR TARGET BLAST RESULTS"
echo The target PCR product matches $BLAST sites in the genome sequence given.
