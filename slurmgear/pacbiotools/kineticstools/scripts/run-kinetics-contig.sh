#!/bin/bash

## ASSUMES ALREADY IN REF_BASENAME_DIR AND THAT THIS IS LAST STEP IN KINETICSPIPELINE (not yet implemented)
## MAKES IPDSUMMARY SUBDIR, ENTERS, ISSUES 1 JOB PER CONTIG

. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27

## VARS
REF=`readlink -f $1` ## IN FUTURE - WHEN IMPLEMENTED IN KINPIPE, REF WILL JUST BE EXPORTED


## DEFAULTS
INPUT=`readlink -f out_all.cmp.h5`
NAMES=contignames.txt
KC=/gpfs_home/jurban/software/kineticsTools/scripts/kinetics-contig.sh
MAX=9 
PVALUE=0.01
IDENTIFY=7
THREADS=4


## 
mkdir -p ipdsummary && cd ipdsummary

## GET TIG NAMES FROM REF
grep ">" ${REF} | awk '{sub(/>/,""); print}' > ${NAMES}
NAMES=`readlink -f $NAMES`


i=0
while read contig; do
  i=$(( $i+1 ))
  if [ $i -eq $MAX ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  if [ ! -d $contig ]; then mkdir $contig; fi
  cd $contig
  sbatch --mem=30g --time=24:00:00 -c 4 -J ${contig}_ipd --qos=${QOS} \
      --export=PVALUE=${PVALUE},IDENTIFY=${IDENTIFY},THREADS=${THREADS},INPUT=${INPUT},REF=${REF},contigName=${contig} ${KC}
  cd ../
done < $NAMES


## LEAVE REF BASENAME_DIR
cd ../

