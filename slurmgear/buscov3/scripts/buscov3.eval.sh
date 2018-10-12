#!/bin/bash

## for use with busco.loop.slurm.in or eval script
## Need: FASTA, OUT, CPU

##module load boost/1.55.0
module load boost/1.55 ## to catch error on new system I have both these lines just in case
module load bamtools/2.3.0
module load blast/2.2.30+

## export august config
export AUGUSTUS_CONFIG_PATH=~/software/busco/augustus-3.2.2/config/

## add augustus to path
export PATH=/gpfs_home/jurban/software/busco/augustus-3.2.2/bin:$PATH

## add HMMer to path
export PATH=~/software/busco/hmmer-3.1b2-linux-intel-x86_64/binaries/:$PATH

## add busco v3 to path
export PATH=/users/jurban/data/software/buscov3/busco/scripts:$PATH



for VAR in FASTA OUT CPU LINEAGE MODE REGIONLIMIT PATH; do echo $VAR ${!VAR}; echo; done; echo

run_BUSCO.py --in $FASTA -o $OUT -l $LINEAGE -m $MODE --cpu $CPU --limit $REGIONLIMIT


