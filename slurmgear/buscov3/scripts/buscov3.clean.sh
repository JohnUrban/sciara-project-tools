#!/bin/bash


for d in augustus augustus_proteins hmmer_output gb gffs single_copy; do
  rm -r run_*/${d} #/*
done

for f in coordinates* full_table* *tblastn missing_buscos_list* training_set*;
  do rm $f
done
