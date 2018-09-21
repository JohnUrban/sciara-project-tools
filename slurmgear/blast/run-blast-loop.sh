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

multiFasta2manyFastas.py -f ${ASM} -d ${TIGDIR}

for fa in ${TIGDIR}/*fa ${TIGDIR}/*fasta; do
 echo $fa
 b=`basename $fa .fa`
 b=`basename $b .fasta`
 sbatch --qos=${QOS} --mem=$MEM --time=$TIME -c $P -o ${SLURMOUT}/${b}.%A.out -J ${b}_blast_NT --export=Q=${fa},P=${P},BLASTDIR=${BLASTDIR},E=${E} blast.sh
done
