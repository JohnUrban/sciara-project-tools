#!/bin/bash
###########################


# SLURMOUTDIR, SLURMPRE, QUERYDIR, FOLLOWUPNUM, BMEM, BTIME, BTHREADS, QOS
# QUERYDIR, PRE, BLASTDIR, BDB, TASK, EVAL, WORDSIZE, CULL, MAXTARGSEQ, SCRIPTS

INCOMPLETE_DIR=tbx_incomplete_${FOLLOWUPNUM}/


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
 if [ $N -ge 1 ] || [ $N2 -lt 3 ] || [ $N3 -ge 1 ]; then
   echo fixing ..............
   echo
   #
   if [ ! -d $INCOMPLETE_DIR ]; then mkdir $INCOMPLETE_DIR; fi
   let NUM_TO_FIX++
   mv $BLASTOUT $INCOMPLETE_DIR
   mv $SLURMFILE $INCOMPLETE_DIR

   BDONE=`sbatch -J ${JOBPRE}${BASE}_tblastx_${i}_followup_${FOLLOWUPNUM} -o ${SLURMOUTDIR}/tblastx.slurm.followup${FOLLOWUPNUM}.%A_${i}.out \
      --mem=$BMEM --time=$BTIME -c $BTHREADS --qos=$QOS \
      --export=QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},P=${BTHREADS},BDB=${BDB},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ},JOBNUM=${i},TBLASTEXTRA="${TBLASTEXTRA}" \
      ${SCRIPTS}/transtblastx.sh | awk '{print $4}'`

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
  FOLLOWDONE=`sbatch --dependency=${BLASTDONE} -J ${JOBPRE}${BASE}_tblastx_followup_${FOLLOWUPNUM} \
     -o ${SLURMOUTDIR}/tbx_follow_up_${FOLLOWUPNUM}.slurm.%A.out \
     --mem=2g --time=6:00:00 -c 2 --qos=$QOS \
     --export=BASE=${BASE},NJOBS=${NJOBS},SLURMOUTDIR=${SLURMOUTDIR},SLURMPRE=tblastx.slurm,FOLLOWUPNUM=${FOLLOWUPNUM},BMEM=${BMEM},BTIME=${BTIME},BTHREADS=${BTHREADS},QOS=${QOS},SCRIPTS=${SCRIPTS},QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},BDB=${BDB},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ},JOBPRE=${JOBPRE},TBLASTEXTRA="${TBLASTEXTRA}" \
     ${SCRIPTS}/followup-tbx.sh | awk '{print $4}'`
elif [ $NUM_TO_FIX -eq 0 ]; then
  #sbatch
  echo MOVE ON TO BLASTOUT ANALYSIS
  D=analysis-tbx
  if [ ! -d $D ]; then mkdir $D; fi
  cd $D
  ## there should not be any dependencies by definition at this point
  FOLLOWDONE=`sbatch -J ${JOBPRE}${BASE}_tblastx_analysis \
     -o ${SLURMOUTDIR}/tbx_analysis.slurm.%A.out \
     --mem=2g --time=6:00:00 -c 2 --qos=$QOS \
     --export=NJOBS=${NJOBS},BLASTDIR=${BLASTDIR} \
     ${SCRIPTS}/analysis.sh | awk '{print $4}'`
fi


