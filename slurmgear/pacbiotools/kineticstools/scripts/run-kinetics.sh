#!/bin/bash

. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27

KINPIPE=~/software/kineticsTools/scripts/kineticspipeline.sh 
REF=/gpfs/scratch/jurban/methylation/ecoli/ecoli_k12_mg1655_NC_000913.3.fasta
QOS=epscor-condo
CLEAN=false
CONFIG=/gpfs/scratch/jurban/methylation/ecoli/kineticsTools/pipeline/output07/kineticspipeline.config.cfg
PVALUE=0.01
INPUT=/gpfs_home/jurban/software/kineticsTools/scripts/input.fofn #input.fofn
###INPUT=/gpfs/scratch/jurban/methylation/ecoli/kineticsTools/input.fofn
IDENTIFY=7
##specify integer 1-7; 1=m6A, 2=m4C, 3=m5C_TET, 4=m6A,m4C, 5=m6A,m5C_TET, 6=m4C,m5C_TET, 7=m6A,m4C,m5C_TET

$KINPIPE $REF $QOS $CLEAN $CONFIG $PVALUE $INPUT $IDENTIFY
