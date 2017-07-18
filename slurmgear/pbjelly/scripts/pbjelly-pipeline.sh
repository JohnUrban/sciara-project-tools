#!/bin/bash
###########################

## $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF $R1 $R2

if [ $# -eq 0 ]; then echo "
Arg1 = scripts dir
Arg2 = config file
Arg3 = Logical(true/false) should directories be cleaned up...
Arg4 = QOS
Arg5 = /Path/To/Reference.fasta
"; exit; fi

MAIN=$PWD

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
QOS=$4
ASM=$5


BASE=`basename $ASM .fasta`

source $CONFIG

#Dir exist?
if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR



## make Protocol.xml
READS_BASEDIR=`abspath.py ${READS} | awk '{print $1}'`
READS_BASENAME=`basename ${READS}`
${SCRIPTS}/get-pbjelly-protocol.py -r ${ASM} -o ${PWD} -b ${READS_BASEDIR} -f ${READS_BASENAME} --minMatch ${minMatch} --sdpTupleSize ${sdpTupleSize} --minPctIdentity ${minPctIdentity} --bestn ${bestn} --nCandidates ${nCandidates} --maxScore ${maxScore} --nproc ${nproc} > Protocol.xml

##############################################################################
## RUN PBJELLY PIPELINE
##############################################################################
MAKEDONE=`sbatch -J ${BASE}_make_hisat_idx -o ${OUT}/make_hisat_idx.slurm.%A.out --mem=$IMEM --time=$ITIME -c $ITHREADS --qos=$QOS \
  --export=G=${ASM},PRE=asm ${SCRIPTS}/hisatbuild.sh | awk '{print $4}'`


