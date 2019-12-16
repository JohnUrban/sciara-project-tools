#!/bin/bash

module load blast/2.2.30+ 

NT=/gpfs/scratch/shared/blastdb/nt
BASE=`basename $Q .fa`

for var in BASE NT Q P BLASTDIR E; do
  echo ${var} ${!var}
done; echo


date
echo "
blastn -task megablast -query $Q -db $NT \
 -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
 -culling_limit 5 \
 -num_threads $P \
 -evalue ${E} \
 -out ${BLASTDIR}/${BASE}.blastout
"

blastn -task megablast -query $Q -db $NT \
 -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
 -culling_limit 5 \
 -num_threads $P \
 -evalue ${E} \
 -out ${BLASTDIR}/${BASE}.blastout

date
