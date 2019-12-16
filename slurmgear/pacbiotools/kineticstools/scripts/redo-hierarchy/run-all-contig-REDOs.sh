#!/bin/bash

while read contig; do
  echo ${contig}
  cp -r contig_template ${contig}
  cd ${contig}
  bash run-REDO.sh
  cd ../
done < contignames.txt 
