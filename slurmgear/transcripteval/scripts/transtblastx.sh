#!/bin/bash

module load blast/2.2.30+ 

## if use "-max_target_seqs 1" you will punish assemblies with more genes split up -- i.e. reward assemblies with more intact genes
## will likely be faster than allowing more....

echo QUERYDIR $QUERYDIR
echo PRE $PRE
echo BLASTDIR $BLASTDIR
echo P $P
echo BDB $BDB
echo EVAL $EVAL
echo WORDSIZE $WORDSIZE
echo CULL $CULL
echo MAXTARGSEQ ${MAXTARGSEQ}
echo JOBNUM $JOBNUM
echo TBLASTEXTRA ${TBLASTEXTRA}
echo; date; echo


CMD="tblastx -strand plus -db $BDB -query ${QUERYDIR}/${PRE}.${JOBNUM}.fa \
 -evalue $EVAL -word_size $WORDSIZE -culling_limit $CULL \
 -max_target_seqs $MAXTARGSEQ -num_threads $P \
 -out ${BLASTDIR}/${PRE}.${JOBNUM}.blastout \
 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sstrand' \
 ${TBLASTEXTRA}"

echo $CMD


tblastx -strand plus -db $BDB -query ${QUERYDIR}/${PRE}.${JOBNUM}.fa \
 -evalue $EVAL -word_size $WORDSIZE -culling_limit $CULL \
 -max_target_seqs $MAXTARGSEQ -num_threads $P \
 -out ${BLASTDIR}/${PRE}.${JOBNUM}.blastout \
 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sstrand' \
 ${TBLASTEXTRA}

echo; date
