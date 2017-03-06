#!/bin/bash
###########################

if [ $# -eq 0 ]; then echo "
Arg1=/Path/To/Reference.fasta 
Arg2=QOS
Arg3=scripts dir path
Arg4=lap reads 1
Arg5=lap reads 2
Arg6=all reads 1
Arg7=all reads 2
Arg8=ConfigFile
"; exit; fi

#FILL THESE IN
#SCRIPTS=/users/jurban/scratch/male-ilmn/long_read_evals/scripts/
#LAPR1=/gpfs/scratch/jurban/lap/sample-1.5m/downsampled.1.fastq
#LAPR2=/gpfs/scratch/jurban/lap/sample-1.5m/downsampled.2.fastq
#R1=/gpfs/scratch/jurban/male-ilmn/data/ilmnraw/R1.fastq
#R2=/gpfs/scratch/jurban/male-ilmn/data/ilmnraw/R2.fastq

## READ IN ARGS
##REF=$1
REF=`readlink -f $1`
QOS=$2
SCRIPTS=$3
LAPR1=$4
LAPR2=$5
R1=$6
R2=$7
CONFIG=$8

source $CONFIG

BASE=`basename $REF .fasta`
JOBSFX=""
ALTQOS=$QOS

MAIN=$PWD


SLURMOUTDIR=slurmout
if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR



##############################################################################
## MAKE BOWTIE2 INDEX
##############################################################################
if $BUILDBT2 || [ ! -d bt2 ]; then
  if [ -d bt2 ]; then rm -r bt2; fi
  mkdir bt2
  cd bt2
  BT2DEP=`sbatch -J ${BASE}_buildbt2${JOBSFX} -o ${OUT}/bt2.slurm.%A.out --mem=$B2MEM --time=$B2TIME -c $B2THREADS --qos=$QOS --export=REF=${REF},BASE=${BASE} ${SCRIPTS}/bt2.eval.sh | awk '{print $4}'`
  cd ../
fi
BT2=`readlink -f bt2/`/$BASE

##############################################################################
## LAP
##############################################################################
if $FULLERASELAP; then if [ -d lap ]; then rm -r lap ; fi; fi
if $RUNLAP; then
 if [ ! -d lap ]; then mkdir lap; fi
 cd lap
 if $BUILDBT2; then
  LAPDONE=`sbatch -J ${BASE}_lap${JOBSFX} --dependency=afterok:${BT2DEP} -o ${OUT}/lap.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=P=${LTHREADS},BASE=${BASE},REF=${REF},R1=${LAPR1},R2=${LAPR2},BT2=${BT2} ${SCRIPTS}/lap.eval.sh | awk '{print $4}'`
 else
  LAPDONE=`sbatch -J ${BASE}_lap${JOBSFX} -o ${OUT}/lap.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=P=${LTHREADS},BASE=${BASE},REF=${REF},R1=${LAPR1},R2=${LAPR2},BT2=${BT2} ${SCRIPTS}/lap.eval.sh | awk '{print $4}'`
  ##same but not --dep
 fi
 cd ../
fi

if $CLEANLAP; then
 cd lap
 # if RUNLAP dep on LAPDONE
 if $RUNLAP; then
   sbatch -J ${BASE}_lap_clean${JOBSFX} --dependency=afterok:${LAPDONE} -o ${OUT}/lap.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/lap.clean.sh 
 # else doesnt dep on LAPDONE
 else
   sbatch -J ${BASE}_lap_clean${JOBSFX} -o ${OUT}/lap.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/lap.clean.sh 
 fi
 cd ../
fi
 

##############################################################################
## MAP READS FOR ALE AND FRC ?
##############################################################################
if $MAPREADS; then
 if [ ! -d mreads ]; then mkdir mreads; fi
 cd mreads
 if $BUILDBT2; then
  ALEFRCDEP=`sbatch -J ${BASE}_mapreads${JOBSFX} --dependency=afterok:${BT2DEP} -o ${OUT}/mapreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=BASE=${BASE},BT2=${BT2},MTHREADS=${MTHREADS},R1=${R1},R2=${R2} ${SCRIPTS}/map.eval.sh | awk '{print $4}'`
 else
  ALEFRCDEP=`sbatch -J ${BASE}_mapreads${JOBSFX} -o ${OUT}/mapreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=BASE=${BASE},BT2=${BT2},MTHREADS=${MTHREADS},R1=${R1},R2=${R2} ${SCRIPTS}/map.eval.sh | awk '{print $4}'`
 fi
 cd ../
fi
BAM=`readlink -f $MAIN/mreads/reads.bam`


##############################################################################
## ALE
##############################################################################
if $RUNALE; then
 DIR=ale
 if [ ! -d $DIR ]; then mkdir $DIR; fi
 cd $DIR
 if $MAPREADS; then
  ALEDONE=`sbatch -J ${BASE}_ale${JOBSFX} --dependency=afterok:${ALEFRCDEP} -o ${OUT}/ale.slurm.%A.out --mem=$AMEM --time=$ATIME -c $ATHREADS --qos=$QOS --export=BASE=${BASE},REF=${REF},BAM=${BAM} ${SCRIPTS}/ale.eval.sh | awk '{print $4}'`
 else
  ALEDONE=`sbatch -J ${BASE}_ale${JOBSFX} -o ${OUT}/ale.slurm.%A.out --mem=$AMEM --time=$ATIME -c $ATHREADS --qos=$QOS --export=BASE=${BASE},REF=${REF},BAM=${BAM} ${SCRIPTS}/ale.eval.sh | awk '{print $4}'`
  ##same but not --dep
 fi
 cd ../
fi
if $CLEANALE; then
 cd ale
 # if RUNALE dep on ALEDONE
 if $RUNALE; then
   sbatch -J ${BASE}_ale_clean${JOBSFX} --dependency=afterok:${ALEDONE} -o ${OUT}/ale.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/ale.clean.sh 
 # else doesnt dep on ALEDONE
 else
   sbatch -J ${BASE}_ale_clean${JOBSFX} -o ${OUT}/ale.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/ale.clean.sh 
 fi 
 cd ../
fi


##############################################################################
## FRC
##############################################################################
if $RUNFRC; then
 DIR=frc
 if [ ! -d $DIR ]; then mkdir $DIR; fi
 cd $DIR
 if $MAPREADS; then
  FRCDONE=`sbatch -J ${BASE}_frc${JOBSFX} --dependency=afterok:${ALEFRCDEP} -o ${OUT}/frc.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=BASE=${BASE},BAM=${BAM} ${SCRIPTS}/frc.eval.sh | awk '{print $4}'`
 else
  FRCDONE=`sbatch -J ${BASE}_frc${JOBSFX} -o ${OUT}/frc.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=BASE=${BASE},BAM=${BAM} ${SCRIPTS}/frc.eval.sh | awk '{print $4}'`
  ##same but not --dep
 fi
 cd ../
fi
if $CLEANFRC; then
 cd frc
 # if RUNFRC, dep on FRCDONE
 if $RUNFRC; then
   sbatch -J ${BASE}_frc_clean${JOBSFX} --dependency=afterok:${FRCDONE} -o ${OUT}/frc.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/frc.clean.sh 
 # else doesnt dep on FRCDONE
 else
   sbatch -J ${BASE}_frc_clean${JOBSFX} -o ${OUT}/frc.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/frc.clean.sh 
 fi
 cd ../
fi


##############################################################################
## REAPR
##############################################################################
## Need: REF, BASE, P, FACHECK, PERFECTMAP, SMALTMAP, PIPELINE
if $RUNREAPR; then
 DIR=reapr
 if [ ! -d $DIR ]; then mkdir $DIR; fi
 cd $DIR
 MEM=$RMEM
 TIME=$RTIME
 THREADS=$RTHREADS
 REAPRDONE=`sbatch -J ${BASE}_reapr${JOBSFX} -o ${OUT}/reapr.slurm.%A.out --mem=$MEM --time=$TIME -c $THREADS --qos=$QOS --export=REF=${REF},BASE=${BASE},R1=${R1},R2=${R2},P=${RTHREADS},FACHECK=${FACHECK},PERFECTMAP=${PERFECTMAP},SMALTMAP=${SMALTMAP},PIPELINE=${PIPELINE},AGGRESSIVE=${AGGRESSIVE} ${SCRIPTS}/reapr.eval.sh | awk '{print $4}'`
 cd ../
fi

if $CLEANREAPR; then
 cd reapr
 # if RUNREAPR, dep on REAPRDONE
 if $RUNREAPR; then
   sbatch -J ${BASE}_reapr_clean${JOBSFX} --dependency=afterok:${REAPRDONE} -o ${OUT}/reapr.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/reapr.clean.sh 
 # else doesnt dep on REAPRDONE
 else
   sbatch -J ${BASE}_reapr_clean -o ${OUT}/reapr.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/reapr.clean.sh 
 fi
 cd ../
fi


##############################################################################
## BUSCO
##############################################################################
## Need: TARGET, OUT, CPU
if $RUNBUSCO; then
 DIR=busco
 if [ ! -d $DIR ]; then mkdir $DIR; fi
 cd $DIR
 MEM=$BMEM
 TIME=$BTIME
 THREADS=$BTHREADS
 BUSCODONE=`sbatch -J ${BASE}_busco${JOBSFX} -o ${OUT}/busco.slurm.%A.out --mem=$MEM --time=$TIME -c $THREADS --qos=$QOS --export=TARGET=${REF},OUT=${BASE},CPU=${BTHREADS} ${SCRIPTS}/busco.eval.sh | awk '{print $4}'`
 cd ../
fi

if $CLEANBUSCO; then
 cd busco
 # if RUNBUSCO, dep on BUSCODONE
 if $RUNBUSCO; then
   sbatch -J ${BASE}_busco_clean${JOBSFX} --dependency=afterok:${BUSCODONE} -o ${OUT}/busco.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/busco.clean.sh 
 # else doesnt dep on BUSCODONE
 else
   sbatch -J ${BASE}_busco_clean${JOBSFX} -o ${OUT}/busco.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/busco.clean.sh 
 fi
 cd ../
fi


## CLEAN MAPPED READS
DEPLIST=afterok
if $BUILDBT2; then DEPLIST=${DEPLIST}:$BT2DEP ; fi
if $MAPREADS; then DEPLIST=${DEPLIST}:$ALEFRCDEP ; fi
if $RUNALE; then DEPLIST=${DEPLIST}:$ALEDONE ; fi
if $RUNFRC; then DEPLIST=${DEPLIST}:$FRCDONE ; fi
echo $DEPLIST > deplist.txt

if $CLEANMREADS; then
 if $BUILDBT2 || $MAPREADS || $RUNALE || $RUNFRC; then
 ## has deplist
   sbatch -J ${BASE}_mreads_clean${JOBSFX} --dependency=${DEPLIST} -o ${OUT}/mreads.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/mreads.clean.sh 
 else ##nodeplist
   sbatch -J ${BASE}_mreads_clean${JOBSFX} -o ${OUT}/mreads.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/mreads.clean.sh 
 fi
fi

### cleaning bt2 index -- add lap dependency
if $RUNLAP; then DEPLIST=${DEPLIST}:$LAPDONE ; fi
if $CLEANBT2; then
 if $BUILDBT2 || $MAPREADS || $RUNALE || $RUNFRC || $RUNLAP; then
 ## has deplist
   sbatch -J ${BASE}_bt2_clean${JOBSFX} --dependency=${DEPLIST} -o ${OUT}/bt2.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/bt2.clean.sh 
 else ##nodeplist
   sbatch -J ${BASE}_bt2_clean${JOBSFX} -o ${OUT}/bt2.clean.slurm.%A.out --mem=8g --time=06:00:00 -c 2 --qos=$QOS ${SCRIPTS}/bt2.clean.sh 
 fi
fi
