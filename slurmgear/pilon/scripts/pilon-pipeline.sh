#!/bin/bash
###########################

if [ $# -eq 0 ]; then echo "
Arg1=/Path/To/Reference.fasta 
Arg2=AsmOutPrefix
Arg3=QOS
Arg4=Logical(true/false) should directories be cleaned up when Pilon is done?
ARG5=config file
Arg6=BUILDBT2 logical
Arg7=MAPREADS logical
Arg8=FLAG1 logical
Arg9=FLAG2 logical
Arg10=MARK mkdup logical
Arg11=PILON logical
Arg12=FIX bases,all,local etc
Arg13=NOSTRAYS logical
Arg14=R1
Arg15=R2
Arg16=scripts loc
Note: nothing will be appended to AsmOutPrefix.
So if you want '.pilon1x.fasta' (for example), then AsmOutPrefix needs to include '.pilon1x'

Exmaple:
$PILONPIPE $ref $asmout $QOS $CLEAN $CONFIG $BUILDBT2 $MAPREADS $FLAG1 $FLAG2 $MARK $PILON $FIX $NOSTRAYS $R1 $R2 $SCRIPTS
"; exit; fi

MAIN=$PWD

ASM=$1
BASE=`basename $ASM .fasta`
OUTPRE=$2 
QOS=$3
CLEAN=$4
CONFIG=$5
BUILDBT2=$6
MAPREADS=$7
FLAG1=$8
FLAG2=$9
MARK=${10}
PILON=${11}
FIX=${12}
NOSTRAYS=${13}
R1=${14}
R2=${15}
SCRIPTS=${16}

## DEBUG CODE
#for v in BASE OUTPRE QOS CLEAN CONFIG BUILDBT2 MAPREADS FLAG1 FLAG2 MARK PILON FIX NOSTRAYS R1 R2 SCRIPTS; do
# echo -e $v "\t" ${!v}
#done
#exit

source $CONFIG

##DEBUG LINES
###echo 1 $ASM 2 $BASE 3 $OUTPRE $QOS $CLEAN $CONFIG
###echo $SCRIPTS $BUILDBT2 $MAPREADS $FLAG1 $MARK $FLAG2 $PILON
###exit

if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

### PIPELINE
CLEAN1DEP=afterok
CLEAN2DEP=afterok
##############################################################################
## MAKE BOWTIE2 INDEX
##############################################################################
if $BUILDBT2 || [ ! -d bt2 ]; then
  if [ -d bt2 ]; then rm -r bt2; fi
  mkdir bt2
  cd bt2
  BT2DEP=`sbatch -J ${BASE}_buildbt2 -o ${OUT}/bt2.slurm.%A.out --mem=$B2MEM --time=$B2TIME -c $B2THREADS --qos=$QOS --export=ASM=${ASM},BASE=${BASE} ${SCRIPTS}/build.sh | awk '{print $4}'`
  cd ../
fi
BT2=`readlink -f bt2/`/$BASE


##############################################################################
## MAP READS 
##############################################################################
if $MAPREADS; then
 if [ ! -d mreads ]; then mkdir mreads; fi
 cd mreads
 if $BUILDBT2; then
  MAPDONE=`sbatch -J ${BASE}_mapreads --dependency=afterok:${BT2DEP} -o ${OUT}/mapreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MEM=${MMEM},BT2=${BT2},P=${MTHREADS},R1=${R1},R2=${R2},PRE=reads  ${SCRIPTS}/map.sh | awk '{print $4}'`
 else
  MAPDONE=`sbatch -J ${BASE}_mapreads -o ${OUT}/mapreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MEM=${MMEM},BT2=${BT2},P=${MTHREADS},R1=${R1},R2=${R2},PRE=reads  ${SCRIPTS}/map.sh | awk '{print $4}'`
 fi
 cd ../
fi
BAM=`readlink -f $MAIN/mreads/reads.bam`



##############################################################################
## FLAG
##############################################################################
if $FLAG1; then
 cd mreads
 if $MAPREADS; then
  F1DONE=`sbatch -J ${BASE}_flag1 --dependency=afterok:${MAPDONE} -o ${OUT}/flag1.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=PRE=reads ${SCRIPTS}/flagstat.sh | awk '{print $4}'`
 else
  F1DONE=`sbatch -J ${BASE}_flag1 -o ${OUT}/flag1.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=PRE=reads ${SCRIPTS}/flagstat.sh | awk '{print $4}'`
  ##same but not --dep
 fi
 cd ../
  CLEAN1DEP=${CLEAN1DEP}:${F1DONE}
fi


##############################################################################
## MARK DUPS
##############################################################################
if $MARK; then
 cd mreads
 if $MAPREADS; then
  MKDONE=`sbatch -J ${BASE}_mkdup --dependency=afterok:${MAPDONE} -o ${OUT}/mkdup.slurm.%A.out --mem=$MDMEM --time=$MDTIME -c $MDTHREADS --qos=$QOS --export=JX=${PICJX},JAR=$PICARDJAR,PRE=reads ${SCRIPTS}/markdups.sh | awk '{print $4}'`
 else
  MKDONE=`sbatch -J ${BASE}_mkdup -o ${OUT}/mkdup.slurm.%A.out --mem=$MDMEM --time=$MDTIME -c $MDTHREADS --qos=$QOS --export=JX=${PICJX},JAR=$PICARDJAR,PRE=reads ${SCRIPTS}/markdups.sh | awk '{print $4}'`
  ##same but not --dep
 fi
 cd ../
  CLEAN1DEP=${CLEAN1DEP}:${MKDONE}
fi
BAM=`readlink -f $MAIN/mreads/reads.markdup.bam`

if $CLEAN && $MARK; then
  cd mreads
  CRDONE=`sbatch -J ${BASE}_cleanreads --dependency=${CLEAN1DEP} -o ${OUT}/cleanreads.slurm.%A.out --mem=2g --time=2:00:00 -c 2 --qos=$QOS ${SCRIPTS}/cleanreads.sh | awk '{print $4}'`
  ##cleanreads.sh
  cd ../
fi

##############################################################################
## FLAG2
##############################################################################
if $FLAG2; then
 cd mreads
 if $MARK; then
  F2DONE=`sbatch -J ${BASE}_flag2 --dependency=afterok:${MKDONE} -o ${OUT}/flag2.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=PRE=reads.markdup ${SCRIPTS}/flagstat.sh | awk '{print $4}'`
 else
  F2DONE=`sbatch -J ${BASE}_flag2 -o ${OUT}/flag2.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=PRE=reads.markdup ${SCRIPTS}/flagstat.sh | awk '{print $4}'`
  ##same but not --dep
 fi
 cd ../
fi

##############################################################################
## PILON
##############################################################################
if $PILON; then
 DIR=pilon
 if [ ! -d $DIR ]; then mkdir $DIR; fi
 cd $DIR
 if $MARK; then
  PDONE=`sbatch -J ${BASE}_pilon --dependency=afterok:${MKDONE} -o ${OUT}/pilon.slurm.%A.out --mem=$PMEM --time=$PTIME -c $PTHREADS --qos=$QOS --export=NOSTRAYS=${NOSTRAYS},JX=${PILJX},PILONJAR=${PILONJAR},ASM=${ASM},READS=${BAM},PRE=${OUTPRE} ${SCRIPTS}/pilon.sh | awk '{print $4}'`
 else
  PDONE=`sbatch -J ${BASE}_pilon -o ${OUT}/pilon.slurm.%A.out --mem=$PMEM --time=$PTIME -c $PTHREADS --qos=$QOS --export=NOSTRAYS=${NOSTRAYS},JX=${PILJX},PILONJAR=${PILONJAR},ASM=${ASM},READS=${BAM},PRE=${OUTPRE} ${SCRIPTS}/pilon.sh | awk '{print $4}'`
  ##same but not --dep
 fi
 cd ../
fi



## CLEAN MAPPED READS

if $CLEAN && $PILON; then
  cd mreads
  CRDONE=`sbatch -J ${BASE}_cleanreads2 --dependency=afterok:${PDONE} -o ${OUT}/cleanreads2.slurm.%A.out --mem=2g --time=2:00:00 -c 2 --qos=$QOS --export=OUTPRE=${OUTPRE} ${SCRIPTS}/cleanreads2.sh | awk '{print $4}'`
  ##cleanreads2.sh
  cd ../
fi

