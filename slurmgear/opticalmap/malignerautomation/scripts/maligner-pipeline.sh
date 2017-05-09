#!/bin/bash
###########################

if [ $# -eq 0 ]; then echo "
Arg1=/Path/To/Reference.fasta
Arg2=maps.fofn -- paths	to all bionano maps files
Arg3=QOS
Arg4=Logical(true/false) should directories be cleaned up...                
ARG5=config file
Arg6=Rec_ENz == eg BssSI
Arg7=Rec_Seq == eg CACGAG
Arg8=scripts dir
"; exit; fi

MAIN=$PWD

ASM=$1
BASE=`basename $ASM .fasta`
FOFN=$2
QOS=$3
CLEAN=$4
CONFIG=$5
REC_ENZ=$6
REC_SEQ=$7
SCRIPTS=$8

source $CONFIG

if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

### PIPELINE
CLEAN1DEP=afterok
CLEAN2DEP=afterok

##############################################################################
## FASTA ASM TO SMOOTHED MAPS
##############################################################################
D=asm_map
if $CONVERTASM; then
  if [ -d $D ]; then rm -r $D; fi
  mkdir $D
  cd $D
  CONVDEP=`sbatch -J ${BASE}_convertasm -o ${OUT}/convertasm.slurm.%A.out --mem=8g --time=2:00:00 -c 2 --qos=$QOS --export=ASM_FASTA=${ASM},BASE=${BASE},REC_ENZ=${REC_ENZ},REC_SEQ=${REC_SEQ} ${SCRIPTS}/fa2map.sh | awk '{print $4}'`
  cd ../
fi
ASM_MAP=`readlink -f asm_map/`/${BASE}.${REC_ENZ}.smoothed.maps


##############################################################################
## MAP PRE-CONVERTED/PRE-SMOOTHED BIONANO MAPS
##############################################################################
if $MAPBIONANO; then
 D=aln
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 DEPENDS=""
 if $CONVERTASM; then
   DEPENDS=--dependency=afterok:${CONVDEP}
 fi
 while read mapfilename; do
   b=`basename $mapfilename`
   CLEAN1DEP=$CLEAN1DEP:`sbatch -J ${BASE}_${b} $DEPENDS -o ${OUT}/${b}.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=ASM_MAP=${ASM_MAP},SMOOTH_MAPS=${mapfilename},REC_SEQ=${REC_SEQ} ${SCRIPTS}/maligner_dp.sh | awk '{print $4}'`
 done < $FOFN 
 cd ../
fi
ALN=`readlink -f aln/`


##############################################################################
## MERGE
##############################################################################
if $MERGE; then
 D=merge
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 DEPENDS=""
 if $MAPBIONANO; then
   DEPENDS=--dependency=${CLEAN1DEP}
 fi
 MERGEDEP=`sbatch -J ${BASE}_merge $DEPENDS -o ${OUT}/merge.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=ALN=${ALN} ${SCRIPTS}/merge.sh | awk '{print $4}'`
 CLEAN1DEP=${CLEAN1DEP}:${MERGEDEP}
 cd ../
fi

ALL=`readlink -f merge/all.bionano.smoothed.maps.aln`


##############################################################################
## SCORE
##############################################################################
if $MERGE; then
 D=merge
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 DEPENDS=""
 if $MERGE; then
   DEPENDS=--dependency=afterok:${MERGEDEP}
 fi
 CLEAN1DEP=$CLEAN1DEP:`sbatch -J ${BASE}_score $DEPENDS -o ${OUT}/score.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=ALL=${ALL} ${SCRIPTS}/score.sh | awk '{print $4}'`
 cd ../
fi

##############################################################################
## BEDGRAPH
##############################################################################
if $MERGE; then
 D=merge
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 DEPENDS=""
 if $MERGE; then
   DEPENDS=--dependency=afterok:${MERGEDEP}
 fi
 CLEAN1DEP=$CLEAN1DEP:`sbatch -J ${BASE}_score $DEPENDS -o ${OUT}/score.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=ALL=${ALL},ASM=${ASM},BASE=${BASE} ${SCRIPTS}/convert2bedGraph.sh | awk '{print $4}'`
 cd ../
fi





##############################################################################
## BEDGRAPH
##############################################################################
if $CLEAN; then
  CLEANDONE=`sbatch -J ${BASE}_bionano_clean --dependency=$CLEANDEP1 -o ${OUT}/clean.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=${QOS} ${SCRIPTS}/clean.sh | awk '{print $4}'`
fi
  
