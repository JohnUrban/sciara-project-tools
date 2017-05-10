#!/bin/bash
###########################


# SLURMOUTDIR, SLURMPRE, QUERYDIR, FOLLOWUPNUM, BMEM, BTIME, BTHREADS, QOS
# QUERYDIR, PRE, BLASTDIR, BDB, TASK, EVAL, WORDSIZE, CULL, MAXTARGSEQ, SCRIPTS

INCOMPLETE_DIR=incomplete_${FOLLOWUPNUM}/


NUM_TO_FIX=0
BLASTDONE=afterany
for i in $(seq $NJOBS); do 
 SLURMFILE=${SLURMOUTDIR}/${SLURMPRE}*_${i}.out
 FA=${QUERYDIR}/${PRE}.${i}.fa
 BLASTOUT=${BLASTDIR}/${PRE}.${i}.blastout
 N=`grep -c slurm $SLURMFILE`
 N2=`grep -c EDT $SLURMFILE` #count how many times date stamped, should be 3
 N3=`grep -c error $SLURMFILE`
 echo $i $SLURMFILE $BLASTOUT $FA
 echo $N slurm, $N2 EDT
 echo
 if [ $N -ge 1 ] || [ $N2 -lt 3 ] || [ $N3 -ge 1  ]; then
   echo fixing ..............
   echo
   #
   if [ ! -d $INCOMPLETE_DIR ]; then mkdir $INCOMPLETE_DIR; fi
   let NUM_TO_FIX++
   mv $BLASTOUT $INCOMPLETE_DIR
   mv $SLURMFILE $INCOMPLETE_DIR

   TASK=blastn
   BDONE=`sbatch -J ${JOBPRE}${BASE}_blast_trans_${i}_followup_${FOLLOWUPNUM} -o ${SLURMOUTDIR}/blast_trans_followup${FOLLOWUPNUM}.slurm.%A_${i}.out \
     --mem=$BMEM --time=$BTIME -c $BTHREADS --qos=$QOS \
     --export=QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},P=${BTHREADS},BDB=${BDB},TASK=${TASK},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ},JOBNUM=${i} \
     ${SCRIPTS}/transblast.sh | awk '{print $4}'`
   
   BLASTDONE+=":"$BDONE

 fi
done


## increment FOLLOWUPNUM
let FOLLOWUPNUM++

## IF there's any jobs to follow up on; then issue new followup dependent on all previous jobs finishing
## ELSE issue blastout analysis

if [ $NUM_TO_FIX -gt 0 ]; then
  #sbatch ...
  echo FIX MORE
  FOLLOWDONE=`sbatch --dependency=${BLASTDONE} -J ${JOBPRE}${BASE}_blast_trans_followup_${FOLLOWUPNUM} \
     -o ${SLURMOUTDIR}/follow_up_${FOLLOWUPNUM}.slurm.%A.out \
     --mem=2g --time=6:00:00 -c 2 --qos=$QOS \
     --export=BASE=${BASE},NJOBS=${NJOBS},SLURMOUTDIR=${SLURMOUTDIR},SLURMPRE=blast_trans,FOLLOWUPNUM=${FOLLOWUPNUM},BMEM=${BMEM},BTIME=${BTIME},BTHREADS=${BTHREADS},QOS=${QOS},SCRIPTS=${SCRIPTS},QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},BDB=${BDB},TASK=${TASK},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ},JOBPRE=${JOBPRE} \
     ${SCRIPTS}/followup.sh | awk '{print $4}'`
elif [ $NUM_TO_FIX -eq 0 ]; then
  #sbatch
  echo MOVE ON TO BLASTOUT ANALYSIS
  D=analysis
  if [ ! -d $D ]; then mkdir $D; fi
  cd $D
  ## there should not be any dependencies by definition at this point
  FOLLOWDONE=`sbatch -J ${JOBPRE}${BASE}_trans_blastout_analysis \
     -o ${SLURMOUTDIR}/analysis.slurm.%A.out \
     --mem=2g --time=6:00:00 -c 2 --qos=$QOS \
     --export=NJOBS=${NJOBS},BLASTDIR=${BLASTDIR} \
     ${SCRIPTS}/analysis.sh | awk '{print $4}'`
fi


