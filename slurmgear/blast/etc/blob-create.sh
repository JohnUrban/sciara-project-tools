#!/bin/bash

source ~/software/blobtools/source.sh

origfa=~/data/scratch/male-ilmn/long_read_asms/quiver0x/canu.corcov500minrl500.aspb-e02.pball.ont2d.fasta

#canu-to-cov.sh $origfa > canu01.cov


##blobtools create -i $REF --cov $file --nodes $NODES --names $NAMES -o $OUT
NODES=/gpfs/scratch/jurban/taxdb/nodes.dmp
NAMES=/gpfs/scratch/jurban/taxdb/names.dmp
file=/gpfs/data/sgerbi/jurban/scratch/methylation/sciara/identify_bacterial_contigs/canu.pball.ont2d/blob/canu01.cov
REF=/gpfs/data/sgerbi/jurban/scratch/methylation/sciara/identify_bacterial_contigs/canu.pball.ont2d/asm/1.canu.corcov500minrl500.aspb-e02.pball.ont2d.quiverfinal1.pilon2x.fasta
BLASTFILE=/gpfs/data/sgerbi/jurban/scratch/methylation/sciara/identify_bacterial_contigs/canu.pball.ont2d/blob/canu01.cov

OUT=canu01.blast
sbatch -J blobcreate01 -o blobcreate01.%A.slurm.out -c 1 --mem=30g --time=12:00:00 --export=REF=${REF},file=${file},NODES=${NODES},NAMES=${NAMES},OUT=${OUT},BLASTFILE=${BLASTFILE} $SCRIPTS/create.sh

