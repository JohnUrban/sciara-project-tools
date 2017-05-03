#!/bin/bash
###########################

## $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF $QUERYDIR $PRE $NJOBS

if [ $# -eq 0 ]; then echo "
Arg1 = scripts dir
Arg2 = config file
Arg3 = Logical(true/false) should directories be cleaned up...
Arg4 = QOS
Arg5 = /Path/To/Reference.fasta
Arg6 = Query Dir
Arg7 = Prefix to fasta files in query dir
Arg8 = Num Jobs -- equal to num files in query dir
"; exit; fi
##TRANSQUERYFOFN - file of filenames -- paths to each FASTA to be used as a query file in BLAST.

MAIN=$PWD

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
QOS=$4
ASM=$5
QUERYDIR=$6
PRE=$7
NJOBS=$8

BASE=`basename $ASM .fasta`

source $CONFIG

#Dir exist?
if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

if [ ! -d $BLASTOUTDIR ]; then
  mkdir $BLASTOUTDIR
fi
BLASTDIR=`readlink -f $MAIN`/$BLASTOUTDIR



### PIPELINE
CLEAN1DEP=afterok
CLEAN2DEP=afterok


##############################################################################
## MAKE BLAST DB FROM ASM
##############################################################################
###makeblastdb -in $ASM -dbtype $DBTYPE -out $OUT
D=blastdb
if $MAKEBLASTDB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 MAKEDONE=`sbatch -J ${BASE}_makeblastdb -o ${OUT}/makeblastdb.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=ASM=${ASM},DBTYPE=nucl,OUT=asm $SCRIPTS/makeblastdb.sh | awk '{print $4}'`
 cd ../
fi
BDB=$(echo `readlink -f ${MAIN}/${D}`/asm)


##############################################################################
## BLAST
##############################################################################
TASK=blastn
BLASTDONE=`sbatch --dependency=afterok:${MAKEDONE} -a 1-$NJOBS -J ${BASE}_blast_trans -o ${OUT}/blast_trans.slurm.%A_%a.out --mem=$BMEM --time=$BTIME -c $BTHREADS --qos=$QOS \
   --export=QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},P=${BTHREADS},BDB=${BDB},TASK=${TASK},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ} \
   ${SCRIPTS}/transblast.sh | awk '{print $4}'`



##############################################################################
## FOLLOW UP 1
##############################################################################
# SLURMOUTDIR, SLURMPRE, QUERYDIR, FOLLOWUPNUM, BMEM, BTIME, BTHREADS, QOS
# QUERYDIR, PRE, BLASTDIR, BDB, TASK, EVAL, WORDSIZE, CULL, MAXTARGSEQ, SCRIPTS
##echo $BLASTDONE
FOLLOWUPNUM=1
FOLLOWDONE=`sbatch --dependency=afterok:${BLASTDONE} -J ${BASE}_blast_trans_followup_${FOLLOWUPNUM} \
   -o ${OUT}/blast_trans_follow_up_${FOLLOWUPNUM}.slurm.%A.out \
   --mem=2g --time=6:00:00 -c 2 --qos=$QOS \
   --export=SLURMOUTDIR=${OUT},SLURMPRE=blast_trans.slurm,FOLLOWUPNUM=${FOLLOWUPNUM},BMEM=${BMEM},BTIME=${BTIME},BTHREADS=${BTHREADS},QOS=${QOS},SCRIPTS=${SCRIPTS},QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},BDB=${BDB},TASK=${TASK},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ} \
   ${SCRIPTS}/followup.sh | awk '{print $4}'`
