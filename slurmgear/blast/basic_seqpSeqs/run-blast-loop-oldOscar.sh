#!/bin/bash

ASM=

TIGDIR=contigs/
BLASTDIR=blast
SLURMOUT=slurmout
P=8
MEM=30g
TIME=24:00:00
QOS=biomed-sb-condo
##E=1e-25
E=1e-10

mkdir -p $BLASTDIR
mkdir -p $SLURMOUT
mkdir -p $TIGDIR

## Separate assembly into its components - 1 per file here.
multiFasta2manyFastas.py -f ${ASM} -d ${TIGDIR}  ## output is all *.fa

## BLAST each component separately
for fa in ${TIGDIR}/*fa; do
 echo $fa
 b=`basename $fa .fa`
 sbatch --qos=${QOS} --mem=$MEM --time=$TIME -c $P -o ${SLURMOUT}/${b}.%A.out -J ${b}_blast_NT --export=Q=${fa},P=${P},BLASTDIR=${BLASTDIR},E=${E} blastn.sh
done
