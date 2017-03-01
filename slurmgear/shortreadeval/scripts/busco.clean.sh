#!/bin/bash


for d in augustus augustus_proteins hmmer_output gb gffs single_copy; do
  rm -r run_*/${d} #/*
done
