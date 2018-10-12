#!/bin/bash

## for use with busco.loop.slurm.in or eval script
## Need: TARGET, OUT, CPU

module load python/3.5.2 ## In case 3.4.1 below is not on RedHatOscar
module load python/3.4.1 
module load boost/1.55.0
module load boost/1.55 ## need both lines for now to work on 2 different systems
module load bamtools/2.3.0
module load blast/2.2.30+
export PATH=/users/jurban/software/busco/BUSCO_v1.22/:~/software/busco/hmmer-3.1b2-linux-intel-x86_64/binaries/:/gpfs_home/jurban/software/busco/augustus-3.2.2/bin:$PATH
export AUGUSTUS_CONFIG_PATH=~/software/busco/augustus-3.2.2/config/

##LINEAGE=~/software/busco/arthropoda
LINEAGE=~/data/busco/lineages/arthropoda
MODE=genome

echo 1 $TARGET
echo 2 $OUT
echo 3 $CPU
echo 4 $LINEAGE
echo 5 $MODE

BUSCO_v1.22.py -in $TARGET -o $OUT -l $LINEAGE -m $MODE --cpu $CPU 



