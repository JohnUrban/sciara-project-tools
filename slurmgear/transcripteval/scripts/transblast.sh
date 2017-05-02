#!/bin/bash

## if use "-max_target_seqs 1" you will punish assemblies with more genes split up -- i.e. reward assemblies with more intact genes
## will likely be faster than allowing more....

EVAL=1e-1
WORDSIZE=8
CULL=1
MAXTARGSEQ=1

blastn -task blastn -db $db -query $TRANSCRIPTS -evalue $1e-1 -word_size $WORDSIZE -culling_limit $CULL -max_target_seqs $MAXTARGSEQ \
   -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sstrand' > ${PRE}.blastn




