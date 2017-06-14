#!/bin/bash


for d in augustus_output blast_output hmmer_output single_copy_busco_sequences; do
  rm -r run_*/${d} 
done

rm -r tmp/
