#!/bin/bash

if [ -f reads.bam ] || [ -f reads.markdup.bam ]; then
    echo Reads detected....
    if [ -f ../pilon/${OUTPRE}.fasta ] && [ -f ../pilon/${OUTPRE}.changes ]; then
        echo Pilon files detected....
        c1=`cat ../pilon/${OUTPRE}.fasta | wc -l`
        c2=`cat ../pilon/${OUTPRE}.changes | wc -l`
        if [ $c1 -gt 0 ] && [ $c2 -gt 0 ]; then
            echo Pilon files are non-empty.... enough evidence to remove reads... removing...
            rm reads*bam*
            cd ../bt2/
            rm *
        else
            echo Pilon files are empty.... keeping reads
        fi
    else
        echo Pilon files not detected.... keeping reads
    fi
else
    echo Reads not detected.... keeping reads
fi
