#!/bin/bash
BT2=$1
P1=$2
P2=$3

if [ $# -eq 0 ]; then
  echo 
  echo Usage: 
  echo testPrimersAndProduct.sh bt2index primer1.fa primer2.fa blastn_index genome_sequence.fa
  echo
  exit 0
fi

primerSpecificity.sh $1 $2 $3
productSpecificity.sh $1 $2 $3 $4 $5
