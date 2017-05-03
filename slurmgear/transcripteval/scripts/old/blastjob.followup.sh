#!/bin/bash
module load blast/2.2.30+ 

echo $i, $fa
#Q=${TEMPDIR}/${PRE}.${i}.fa
NT=/gpfs/scratch/shared/blastdb/nt
#echo $Q

date
echo "
blastn -task megablast -query $fa -db $NT \
 -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
 -culling_limit 5 \
 -num_threads $P \
 -evalue 1e-25 \
 -out ${BLASTDIR}/${PRE}.${i}.blastout
"

blastn -task megablast -query $fa -db $NT \
 -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
 -culling_limit 5 \
 -num_threads $P \
 -evalue 1e-25 \
 -out ${BLASTDIR}/${PRE}.${i}.followup.blastout

date
