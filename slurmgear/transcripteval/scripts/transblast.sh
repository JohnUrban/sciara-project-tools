#!/bin/bash

## if use "-max_target_seqs 1" you will punish assemblies with more genes split up -- i.e. reward assemblies with more intact genes
## will likely be faster than allowing more....

echo QUERYDIR $QUERYDIR
echo PRE $PRE
echo BLASTDIR $BLASTDIR
echo P $P
echo BDB $BDB
echo TASK $TASK
echo EVAL $EVAL
echo WORDSIZE $WORDSIZE
echo CULL $CULL
echo MAXTARGSEQ ${MAXTARGSEQ}



CMD="blastn -task $TASK -db $BDB -query ${QUERYDIR}/${PRE}.${SLURM_ARRAY_TASK_ID}.fa 
 -evalue $1e-1 -word_size $WORDSIZE -culling_limit $CULL \
 -max_target_seqs $MAXTARGSEQ -num_threads $P \
 -out ${BLASTDIR}/${PRE}.${SLURM_ARRAY_TASK_ID}.blastout \
 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sstrand' "

echo $CMD

$CMD


