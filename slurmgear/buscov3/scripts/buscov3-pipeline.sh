#!/bin/bash
###########################
##$PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF

if [ $# -eq 0 ]; then echo "
Arg1 = SCRIPTS
Arg2 = CONFIG
Arg3 = CLEAN
Arg4 = QOS
Arg5 = REF
"; exit; fi


MAIN=$PWD

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
QOS=$4
ASM=$5
BASE=`basename $ASM .fasta`

source $CONFIG

if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

### PIPELINE

##############################################################################
## BUSCO V3 LOOP
##############################################################################

for DIR in metazoa arthropoda insecta endopterygota diptera; do

 L=`echo $DIR | awk '{print toupper($0)}'`
 LINEAGE=${!L}
 RUNVAR=RUN_${L}
 RUN=${!RUNVAR}

 #for debugging: for VAR in L LINEAGE RUN; do echo $VAR ${!VAR}; done

 if $RUN; then
  if [ ! -d $DIR ]; then mkdir $DIR; fi
  cd $DIR
  DONE=`sbatch -J ${BASE}_buscov3_${DIR} -o ${OUT}/${DIR}.slurm.%A.out --mem=$BMEM --time=$BTIME -c $BTHREADS --qos=$QOS --export=FASTA=${REF},OUT=${BASE},CPU=${BTHREADS},LINEAGE=${LINEAGE} ${SCRIPTS}/buscov3.eval.sh | awk '{print $4}'`
  if $CLEAN; then
   echo ....ADD CLEANING FUNCITONALITY....
  fi
  cd ../
 fi
done


