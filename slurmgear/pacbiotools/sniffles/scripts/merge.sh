#!/bin/bash
source ~/software/sniffles/source.sh


samtools merge --threads $P combined.bam $PBBAM $ONTBAM
