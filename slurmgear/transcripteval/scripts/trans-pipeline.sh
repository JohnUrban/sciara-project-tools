#!/bin/bash
###########################

## $PIPELINE $SCRIPTS $CONFIG $CLEAN $TRANSQUERYFOFN $REF $QOS

if [ $# -eq 0 ]; then echo "
Arg1 = scripts dir
Arg2 = config file
Arg3 = Logical(true/false) should directories be cleaned up...
Arg4 = TRANSQUERYFOFN - file of filenames -- paths to each FASTA to be used as a query file in BLAST.
Arg5 = /Path/To/Reference.fasta
Arg6 = QOS
"; exit; fi

MAIN=$PWD

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
FOFN=$4
ASM=$5
QOS=$6

BASE=`basename $ASM .fasta`

source $CONFIG

#Dir exist?
if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

if [ ! -d $BLASTOUTDIR ]; then
  mkdir $BLASTOUTDIR
fi
BLASTDIR=`readlink -f $MAIN`/$BLASTDIR



### PIPELINE
CLEAN1DEP=afterok
CLEAN2DEP=afterok

##############################################################################
## BLAST
##############################################################################


sbatch -a 1-$NJOBS -J ${BASE}_blast_trans -o ${OUT}/blast_trans.slurm.%A_%a.out --mem=$BMEM --time=$BTIME -c $BTHREADS --qos=$QOS --export=TEMPDIR=${TEMPDIR},BLASTDIR=${BLASTDIR},PRE=${PRE},P=${P} blastjob.sh


D=mreads
if $MAPPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $BUILDBWA; then
  MAPPBDONE=`sbatch -J ${BASE}_mapPBreads --dependency=afterok:${IDXDEP} -o ${OUT}/mapPBreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=pacbio,BWAIDX=${BWAIDX},FASTQ=${PACBIO}  ${SCRIPTS}/bwa-map.sh | awk '{print $4}'`
 else
  MAPPBDONE=`sbatch -J ${BASE}_mapPBreads -o ${OUT}/mapPBreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=pacbio,BWAIDX=${BWAIDX},FASTQ=${PACBIO}  ${SCRIPTS}/bwa-map.sh | awk '{print $4}'`
 fi
 cd ../
 CLEANPBDEP=${CLEANPBDEP}:${MAPPBDONE}
 COMBINEDEP=${COMBINEDEP}:${MAPPBDONE}
fi
PBBAM=`readlink -f ${MAIN}/${D}/pacbio.bam`
