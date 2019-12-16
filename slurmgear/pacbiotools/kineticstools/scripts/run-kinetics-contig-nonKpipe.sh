#!/bin/bash

## ASSUMES DOING THIS IN TOPDIR ABOVE REF BASE_NAME DIRS 
##	AS OPPOSED TO ASSUMING ALREADY IN REF_BASENAME_DIR AND THAT THIS IS LAST STEP IN KINETICSPIPELINE 
## MAKES/ENTERS REF_BASENAME_DIR
## 	B/C OF THIS, IT ASSUMES THE FIRST 3 STEPS OF KINETICS PIPELINE ALREADY DONE (PBALIGN, MERGE, SORT) AND EXIST IN REF_BASENAME_DIR
##	ALSO B/C OF THIS, IT CAN BE RUN IN THE TOPDIR IN FOR LOOP OVER REFs
## MAKES IPDSUMMARY SUBDIR, ENTERS, ISSUES 1 JOB PER CONTIG

. /users/jurban/data/software/conda/anaconda3/etc/profile.d/conda.sh
conda activate py27

## VARS
REF=`readlink -f $1` ## IN FUTURE - WHEN IMPLEMENTED IN KINPIPE, REF WILL JUST BE EXPORTED


## DEFAULTS
NAMES=contignames.txt
KC=/gpfs_home/jurban/software/kineticsTools/scripts/kinetics-contig.sh
MAX=9 
PVALUE=0.01
IDENTIFY=7
THREADS=4


## MAKE AND ENTER REF_BASENAME_DIR
DIR=`basename $REF .fasta`
mkdir -p $DIR && cd $DIR

## INPUT SHOULD ALREADY BE IN THIS REF_BASENAME_DIR
INPUT=`readlink -f out_all.cmp.h5`


## MAKE IPDSUMMARY SUBDIR
mkdir -p ipdsummary && cd ipdsummary

## GET TIG NAMES FROM REF
grep ">" ${REF} | awk '{sub(/>/,""); print}' > ${NAMES}
NAMES=`readlink -f $NAMES`


## FOR EACH CONTIG, MAKE SUBSUBDIR, ENTER, ISSUE IPDSUMMARY JOB, LEAVE SUBSUBDIR
REFBASE=`basename $REF .fasta`
i=0
while read contig; do
  i=$(( $i+1 ))
  if [ $i -eq $MAX ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  if [ ! -d $contig ]; then mkdir $contig; fi
  cd $contig
  sbatch --mem=30g --time=24:00:00 -c 4 -J ${REFBASE}_${contig}_ipd --qos=${QOS} \
      --export=PVALUE=${PVALUE},IDENTIFY=${IDENTIFY},THREADS=${THREADS},INPUT=${INPUT},REF=${REF},contigName=${contig} ${KC}
  cd ../
done < $NAMES


## LEAVE IPDSUMARRU SUBDIR AND REF BASENAME_DIR
cd ../../

