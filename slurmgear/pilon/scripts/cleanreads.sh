#!/bin/bash

  if [ -f reads.bam ] && [ -f reads.markdup.bam ]; then
    echo Both files detected
    c=`samtools view -c reads.markdup.bam`
    if [ $c -gt	0 ]; then 
      echo MarkDup count: $c
      rm reads.bam*
    else
      echo Mark Dup BAM empty... count: $c
    fi
  else
    echo ...Did not detect both files....
  fi
