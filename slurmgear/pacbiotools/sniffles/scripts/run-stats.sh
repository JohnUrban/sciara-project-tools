#!/bin/bash


FOFN=input.fofn


##FOFN=debug.fofn

SCRIPTS=~/software/sniffles/scripts
ONTFQ=/users/jurban/scratch/minion2016/fast5fastqs/allReadsFromAllONTlibsCombined.fastq
PBFQ=/gpfs/scratch/jurban/pac_bio_data/filt/all_subreads.fastq

i=0
while read ASM; do 
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  b=`basename $ASM .fasta`
  cd $b
  OUT=`readlink -f slurmout/`
  PBBAM=`readlink -f mreads/pacbio.bam`
  ONTBAM=`readlink -f mreads/ont2d.bam`
  sbatch -J ${b}_snifflestats -o ${OUT}/snifflestats.slurm.%A.out --mem=32g --time=72:00:00 -c 4 --qos=$QOS --export=ASM=${ASM},PBBAM=${PBBAM},PBFQ=${PBFQ},ONTBAM=${ONTBAM},ONTFQ=${ONTFQ} ${SCRIPTS}/do-stuff.sh | awk '{print $4}'
  cd ../
done < $FOFN
  
