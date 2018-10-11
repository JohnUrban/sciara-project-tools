#!/bin/bash
###########################

## $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF $R1 $R2

if [ $# -eq 0 ]; then echo "
Arg1 = scripts dir
Arg2 = config file
Arg3 = Logical(true/false) should directories be cleaned up...
Arg4 = QOS
Arg5 = /Path/To/Reference.fasta
Arg6=READSFOFN
"; exit; fi

MAIN=$PWD

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
QOS=$4
ASM=$5
READSFOFN=$6



BASE=`basename $ASM .fasta`

source $CONFIG

#Dir exist?
if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR




##############################################################################
## MAKE HISAT2 INDEX
##############################################################################
D=ht2
if $MAKEIDX; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 MAKEDONE=`sbatch -J ${BASE}_make_hisat_idx -o ${OUT}/make_hisat_idx.slurm.%A.out --mem=$IMEM --time=$ITIME -c $ITHREADS --account=${QOS} \
   --export=ALL,G=${ASM},PRE=asm ${SCRIPTS}/hisatbuild.sh | awk '{print $4}'`
 cd ../
fi
HIDX=$(echo `readlink -f ${MAIN}/${D}`/asm)


##############################################################################
## HISAT2
##############################################################################
##hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS 2> mapreads.err | samtools view -bSh - | samtools sort -T ${PRE} > ${PRE}.bam

STRANDEDNESS=RF
D=mreads
if $MAPREADS; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 MAPDONE=`sbatch --dependency=afterok:${MAKEDONE} -J ${BASE}_hisat2 -o ${OUT}/hisat.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --account=${QOS} \
   --export=ALL,P=${MTHREADS},HIDX=${HIDX},READSFOFN=${READSFOFN},STRANDEDNESS=${STRANDEDNESS},PRE=rnaseq,CLEAN=${CLEANRNASEQREADS} \
   ${SCRIPTS}/hisat2.sh | awk '{print $4}'`
 cd ../
fi
