#!/bin/bash
#SBATCH -J canu-opt-loop
#SBATCH -c 4
#SBATCH --qos=biomed-sb-condo
#SBATCH --mem=47g
#SBATCH --time=96:00:00
#SBATCH -o slurm.optloop.out
#SBATCH -e slurm.optloop.err
###########################


sh optimization-loop-2.sh 
