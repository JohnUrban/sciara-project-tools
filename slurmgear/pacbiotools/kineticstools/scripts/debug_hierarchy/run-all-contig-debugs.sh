#!/bin/bash

while read contig; do
  echo ${contig}
  cp -r contig_template ${contig}
  cd ${contig}
  bash run-all-debugs.sh
  cd ../
done < contignames.txt 
