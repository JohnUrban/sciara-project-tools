#!/bin/bash

# copy this script into your working dir
# Make sure variables below work for you.
# Make sure to have the ASMFOFN (usually input.fofn) in working dir.

## ASM FOFN
ASMFOFN=input.fofn

## LONG READ LOCATIONS
ONTFQ=~/data/scratch/minion2016/fast5fastqs/allReadsFromAllONTlibsCombined.fastq
PBFQ=~/data/scratch/pac_bio_data/filt/all_subreads.fastq

## RUN INFO LOCATIONS
BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/sniffles
SCRIPTS=${BASE}/scripts/


i=0
while read ASM; do 
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  b=`basename $ASM .fasta`
  cd $b
  OUT=`readlink -f slurmout/`
  PBBAM=`readlink -f mreads/pacbio.bam`
  ONTBAM=`readlink -f mreads/ont2d.bam`
  COMBBAM=`readlink -f mreads/combined.bam`
  PBSNIFF=`readlink -f sniffles_pb/*_pacbio.bedpe`
  ONTSNIFF=`readlink -f sniffles_ont/*_ont.bedpe`
  COMBSNIFF=`readlink -f sniffles_combined/*_combined.bedpe`
  sbatch -J ${b}_snifflestats -o ${OUT}/snifflestats.slurm.%A.out --mem=32g --time=72:00:00 -c 4 --qos=$QOS --export=ASM=${ASM},PBBAM=${PBBAM},PBFQ=${PBFQ},ONTBAM=${ONTBAM},ONTFQ=${ONTFQ},COMBBAM=${COMBBAM},PBSNIFF=${PBSNIFF},ONTSNIFF=${ONTSNIFF},COMBSNIFF=${COMBSNIFF} ${SCRIPTS}/snifflestats.sh | awk '{print $4}'
  cd ../
done < $ASMFOFN
  
