#!/bin/bash
#SBATCH -J canu-main-g250m-v1
#SBATCH -c 16
#SBATCH --qos=biomed-sb-condo
#SBATCH --mem=50g
#SBATCH --time=48:00:00

###########################

R=/users/jurban/scratch/pac_bio_data/filtered_subreads.fastq

module load java/8u66 

canu/Linux-amd64/bin/canu \
 -p sciara-v1 -d sciara-v1-t2 \
 genomeSize=250m \
 -pacbio-raw $R \
 "gridOptions=--time 96:00:00" \
 "gridOptionsJobName=g250m-t2" \
 "minReadLength=500" \
 merylMemory=30 \
 oeaMemory=8 \
 cnsMemory=30
