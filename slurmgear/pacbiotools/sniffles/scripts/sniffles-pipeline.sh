#!/bin/bash
###########################

if [ $# -eq 0 ]; then echo "
Arg1=/Path/To/Reference.fasta 
Arg2=QOS
Arg3=Logical(true/false) should directories be cleaned up when Pilon is done?
ARG4=config file
"; exit; fi


MAIN=$PWD

ASM=$1
BASE=`basename $ASM .fasta`
QOS=$2
CLEAN=$3
CONFIG=$4

source $CONFIG
source ~/software/bwa/source.sh 
source ~/software/sniffles/source.sh

if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

### PIPELINE
CLEANONTDEP=afterok
CLEANPBDEP=afterok
COMBINEDEP=afterok
##############################################################################
## MAKE BWA INDEX
##############################################################################
D=bwa
if $BUILDBWA || [ ! -d $D ]; then
  if [ -d $D ]; then rm -r $D; fi
  mkdir $D
  cd $D
  IDXDEP=`sbatch -J ${BASE}_buildbwa -o ${OUT}/bwaidx.slurm.%A.out --mem=$BIMEM --time=$BITIME -c $BITHREADS --qos=$QOS --export=ASM=${ASM},BASE=${BASE} ${SCRIPTS}/bwa-idx.sh | awk '{print $4}'`
  cd ../
fi
BWAIDX=`readlink -f bwa/`/$BASE


##############################################################################
## MAP PACBIO READS 
##############################################################################
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

##############################################################################
## MAP ONT READS 
##############################################################################
D=mreads
if $MAPONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $BUILDBWA; then
  MAPONTDONE=`sbatch -J ${BASE}_mapONTreads --dependency=afterok:${IDXDEP} -o ${OUT}/mapONTreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=ont2d,BWAIDX=${BWAIDX},FASTQ=${ONT}  ${SCRIPTS}/bwa-map.sh | awk '{print $4}'`
 else
  MAPONTDONE=`sbatch -J ${BASE}_mapONTreads --dependency=afterok:${IDXDEP} -o ${OUT}/mapONTreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=ont2d,BWAIDX=${BWAIDX},FASTQ=${ONT}  ${SCRIPTS}/bwa-map.sh | awk '{print $4}'`
 fi
 cd ../
 CLEANONTDEP=${CLEANONTDEP}:${MAPONTDONE}
 COMBINEDEP=${COMBINEDEP}:${MAPONTDONE}
fi
ONTBAM=`readlink -f ${MAIN}/${D}/ont2d.bam` ##not necessarily only 2d -- just using the bwa type as an easy prefix



#bwa index $ASM -p $BASE
#bwa mem -t $MTHREADS -M -x $TYPE $BWAIDX $FASTQ | samtools sort --threads $MTHREADS -o $TYPE.bam
#sniffles -m $BAM -b $BEDPE
##############################################################################
## SNIFFLES PACBIO
##############################################################################
D=sniffles_pb
if $SNIFFLESPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPPB; then
  SNIFFLESPBDONE=`sbatch -J ${BASE}_sniffles_pb --dependency=afterok:${MAPPBDONE} -o ${OUT}/snifflesPB.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$PBBAM,BEDPE=${BASE}_pacbio ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 else
  SNIFFLESPBDONE=`sbatch -J ${BASE}_sniffles_pb -o ${OUT}/snifflesPB.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$PBBAM,BEDPE=${BASE}_pacbio ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 fi
 cd ../
 CLEANPBDEP=${CLEANPBDEP}:${SNIFFLESPBDONE}
fi

##############################################################################
## SNIFFLES ONT
##############################################################################
D=sniffles_ont
if $SNIFFLESONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPONT; then
  SNIFFLESONTDONE=`sbatch -J ${BASE}_sniffles_ont --dependency=afterok:${MAPONTDONE} -o ${OUT}/snifflesONT.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$ONTBAM,BEDPE=${BASE}_ont ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 else
  SNIFFLESONTDONE=`sbatch -J ${BASE}_sniffles_ont -o ${OUT}/snifflesONT.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$ONTBAM,BEDPE=${BASE}_ont ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 fi
 cd ../
 CLEANONTDEP=${CLEANONTDEP}:${SNIFFLESONTDONE}
fi


##############################################################################
## SNIFFLES COMBINED
##############################################################################
D=mreads
if $SNIFFLESCOMBINED; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPONT || $MAPPB ; then
    COMBINEDONE=`sbatch -J ${BASE}_merge_reads --dependency=${COMBINEDEP} -o ${OUT}/merge.slurm.%A.out --mem=$CMEM --time=$CTIME -c $CTHREADS --qos=$QOS --export=P=${CTHREADS},PBBAM=${PBBAM},ONTBAM=${ONTBAM} ${SCRIPTS}/merge.sh | awk '{print $4}'`
 else
    COMBINEDONE=`sbatch -J ${BASE}_merge_reads -o ${OUT}/merge.slurm.%A.out --mem=$CMEM --time=$CTIME -c $CTHREADS --qos=$QOS --export=P=${CTHREADS},PBBAM=${PBBAM},ONTBAM=${ONTBAM} ${SCRIPTS}/merge.sh | awk '{print $4}'`
 fi
 cd ../
 COMBBAM=`readlink -f ${MAIN}/${D}/combined.bam` 
 D=sniffles_combined
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 SNIFFLESCOMBDONE=`sbatch -J ${BASE}_sniffles_combined --dependency=afterok:${COMBINEDONE} -o ${OUT}/snifflesCOMBINED.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=${COMBBAM},BEDPE=${BASE}_combined ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 cd ../
 CLEANONTDEP=${CLEANONTDEP}:${SNIFFLESCOMBDONE}
 CLEANPBDEP=${CLEANPBDEP}:${SNIFFLESCOMBDONE}
fi





##############################################################################
## CLEAN UP
##############################################################################
if $CLEAN; then
 if $SNIFFLESPB || $SNIFFLESCOMBINED; then
  GATEFILE="sniffles_pb/${BASE}_pacbio.bedpe"
  DELFILE=$PBBAM
  CLEANPBDONE=`sbatch -J ${BASE}_clean_pb --dependency=${CLEANPBDEP} -o ${OUT}/cleanPB.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`
 fi
 if $SNIFFLESONT || $SNIFFLESCOMBINED; then
  GATEFILE="sniffles_ont/${BASE}_ont.bedpe"
  DELFILE=$ONTBAM
  CLEANPBDONE=`sbatch -J ${BASE}_clean_ont --dependency=${CLEANONTDEP} -o ${OUT}/cleanONT.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`
 fi
 if $SNIFFLESCOMBINED; then
  GATEFILE="sniffles_combined/${BASE}_combined.bedpe"
  DELFILE=$COMBBAM
  CLEANCOMBINEDONE=`sbatch -J ${BASE}_clean_combined --dependency=afterok:${SNIFFLESCOMBDONE} -o ${OUT}/cleanCombined.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'` 
 fi
fi

