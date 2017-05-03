#!/bin/bash
module load blast/2.2.30+ 

echo $SLURM_ARRAY_TASK_ID
Q=${TEMPDIR}/${PRE}.${SLURM_ARRAY_TASK_ID}.fa
NT=/gpfs/scratch/shared/blastdb/nt
echo $Q

date
echo "
blastn -task megablast -query $Q -db $NT \
 -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
 -culling_limit 5 \
 -num_threads $P \
 -evalue 1e-25 \
 -out ${BLASTDIR}/${PRE}.${SLURM_ARRAY_TASK_ID}.blastout
"

blastn -task megablast -query $Q -db $NT \
 -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
 -culling_limit 5 \
 -num_threads $P \
 -evalue 1e-25 \
 -out ${BLASTDIR}/${PRE}.${SLURM_ARRAY_TASK_ID}.blastout

date
