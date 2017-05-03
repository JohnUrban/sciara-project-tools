#!/bin/bash
###########################


# SLURMOUTDIR, SLURMPRE, QUERYDIR, FOLLOWUPNUM, BMEM, BTIME, BTHREADS, QOS
# QUERYDIR, PRE, BLASTDIR, BDB, TASK, EVAL, WORDSIZE, CULL, MAXTARGSEQ, SCRIPTS

INCOMPLETE_DIR=incomplete_${FOLLOWUPNUM}/

NUM_TO_FIX=0

for i in $(seq $NJOBS); do 
 SLURMFILE=${SLURMOUTDIR}/${SLURMPRE}*_${i}.out
 FA=${QUERYDIR}/${PRE}.${i}.fa
 BLASTOUT=${BLASTDIR}/${PRE}.${i}.blastout
 N=`grep -c slurm $SLURMFILE`
 N2=`grep -c EDT $SLURMFILE` #count how many times date stamped, should be 3
 echo $i $SLURMFILE $BLASTOUT $FA
 echo $N slurm, $N2 EDT
 echo
 if [ $N -ge 1 ] || [ $N2 -lt 3 ]; then
   echo fixing ..............
   echo
   #
   let NUM_TO_FIX++
   mv $BLASTOUT $INCOMPLETE_DIR
   mv $SLURMFILE $INCOMPLETE_DIR

   TASK=blastn
   sbatch -J ${BASE}_blast_trans_${i}_followup_${FOLLOWUPNUM} -o ${SLURMOUTDIR}/blast_trans_followup_${FOLLOWUPNUM}.slurm.%A_%a.out \
     --mem=$BMEM --time=$BTIME -c $BTHREADS --qos=$QOS \
     --export=QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},P=${BTHREADS},BDB=${BDB},TASK=${TASK},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ} \
     ${SCRIPTS}/transblast.sh

 fi
done


## increment FOLLOWUPNUM
let FOLLOWUPNUM++

## IF there's any jobs to follow up on; then issue new followup dependent on all previous jobs finishing
## ELSE issue blastout analysis

if [ $NUM_TO_FIX -gt 0 ]; then
  #sbatch ...
  echo FIX MORE
elif [ $NUM_TO_FIX -eq 0 ]; then
  #sbatch
  echo MOVE ON TO BLASTOUT ANALYSIS
fi


