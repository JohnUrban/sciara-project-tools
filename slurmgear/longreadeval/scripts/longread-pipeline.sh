#!/bin/bash
###########################

if [ $# -eq 0 ]; then echo "
Arg1=/Path/To/Reference.fasta 
Arg2=QOS
Arg3=Logical(true/false) should directories be cleaned up when Pilon is done?
ARG4=config file
Arg5=Scripts dir
Arg6=Ont fastq loc
Arg7=PacBio fastq loc
Arg8=Ont long2pe R1 loc
Arg9=Ont long2pe R2 loc
Arg10=Pb long2pe R1 loc
Arg11=Pb long2pe R2 loc

"; exit; fi


MAIN=$PWD

ASM=$1
BASE=`basename $ASM .fasta`
QOS=$2
CLEAN=$3
CONFIG=$4
SCRIPTS=$5
ONT=$6
PACBIO=$7
ONT1=$8
ONT2=$9
PACBIO1=${10}
PACBIO2=${11}


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
  SNIFFLESPBDONE=`sbatch -J ${BASE}_sniffles_pb --dependency=afterok:${MAPPBDONE} -o ${OUT}/snifflesPB.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$PBBAM,BEDPE=${BASE}_pacbio,MINSUPPORT=${MINSUPPORT_PACBIO} ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 else
  SNIFFLESPBDONE=`sbatch -J ${BASE}_sniffles_pb -o ${OUT}/snifflesPB.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$PBBAM,BEDPE=${BASE}_pacbio,MINSUPPORT=${MINSUPPORT_PACBIO} ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
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
  SNIFFLESONTDONE=`sbatch -J ${BASE}_sniffles_ont --dependency=afterok:${MAPONTDONE} -o ${OUT}/snifflesONT.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$ONTBAM,BEDPE=${BASE}_ont,MINSUPPORT=${MINSUPPORT_ONT} ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 else
  SNIFFLESONTDONE=`sbatch -J ${BASE}_sniffles_ont -o ${OUT}/snifflesONT.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=$ONTBAM,BEDPE=${BASE}_ont,MINSUPPORT=${MINSUPPORT_ONT} ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
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
 SNIFFLESCOMBDONE=`sbatch -J ${BASE}_sniffles_combined --dependency=afterok:${COMBINEDONE} -o ${OUT}/snifflesCOMBINED.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=BAM=${COMBBAM},BEDPE=${BASE}_combined,MINSUPPORT=${MINSUPPORT_PACBIO} ${SCRIPTS}/sniffles.sh | awk '{print $4}'`
 cd ../
 CLEANONTDEP=${CLEANONTDEP}:${SNIFFLESCOMBDONE}
 CLEANPBDEP=${CLEANPBDEP}:${SNIFFLESCOMBDONE}
fi



##############################################################################
## ALE SECTION ::: ALE CAN USE SAME READS USED FOR SNIFFLES
##############################################################################

##############################################################################
## ALE PACBIO
##############################################################################
D=ale_pb
if $ALEPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPPB; then
  ALEPBDONE=`sbatch -J ${BASE}_ale_pb --dependency=afterok:${MAPPBDONE} -o ${OUT}/alePB.slurm.%A.out --mem=$AMEM --time=$ATIME -c $ATHREADS --qos=$QOS \
    --export=BAM=$PBBAM,BASE=${BASE},REF=${ASM},CLEAN=${CLEAN},SCRIPTS=${SCRIPTS} ${SCRIPTS}/ale.eval.sh | awk '{print $4}'`
 else
  ALEPBDONE=`sbatch -J ${BASE}_ale_pb -o ${OUT}/alePB.slurm.%A.out --mem=$AMEM --time=$ATIME -c $ATHREADS --qos=$QOS \
    --export=BAM=$PBBAM,BASE=${BASE},REF=${ASM},CLEAN=${CLEAN},SCRIPTS=${SCRIPTS} ${SCRIPTS}/ale.eval.sh | awk '{print $4}'`
 fi
 cd ../
 CLEANPBDEP=${CLEANPBDEP}:${ALEPBDONE}
fi

##############################################################################
## ALE ONT
##############################################################################
D=ale_ont
if $ALEONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPONT; then
  ALEONTDONE=`sbatch -J ${BASE}_ale_ont --dependency=afterok:${MAPONTDONE} -o ${OUT}/aleONT.slurm.%A.out --mem=$AMEM --time=$ATIME -c $ATHREADS --qos=$QOS \
    --export=BAM=$ONTBAM,BASE=${BASE},REF=${ASM},CLEAN=${CLEAN},SCRIPTS=${SCRIPTS} ${SCRIPTS}/ale.eval.sh | awk '{print $4}'`
 else
  ALEONTDONE=`sbatch -J ${BASE}_ale_ont -o ${OUT}/aleONT.slurm.%A.out --mem=$AMEM --time=$ATIME -c $ATHREADS --qos=$QOS \
    --export=BAM=$ONTBAM,BASE=${BASE},REF=${ASM},CLEAN=${CLEAN},SCRIPTS=${SCRIPTS} ${SCRIPTS}/ale.eval.sh | awk '{print $4}'`
 fi
 cd ../
 CLEANONTDEP=${CLEANONTDEP}:${ALEONTDONE}
fi



##############################################################################
## STATS
##NOTE: Will need to come re-visit the entire stats section to allow independence of pb and ont...
##############################################################################

STATDEP=${SNIFFLESPBDONE}:${SNIFFLESONTDONE}:${SNIFFLESCOMBDONE}
PBSNIFF=`readlink -f sniffles_pb/*_pacbio.bedpe`
ONTSNIFF=`readlink -f sniffles_ont/*_ont.bedpe`
COMBSNIFF=`readlink -f sniffles_combined/*_combined.bedpe`
STATSDONE=`sbatch -J ${BASE}_snifflestats --dependency=afterok:${STATDEP} -o ${OUT}/snifflestats.slurm.%A.out --mem=32g --time=72:00:00 -c 4 --qos=$QOS --export=ASM=${ASM},PBBAM=${PBBAM},PBFQ=${PACBIO},ONTBAM=${ONTBAM},ONTFQ=${ONT},COMBBAM=${COMBBAM},PBSNIFF=${PBSNIFF},ONTSNIFF=${ONTSNIFF},COMBSNIFF=${COMBSNIFF} ${SCRIPTS}/longreadstats.sh | awk '{print $4}'`

CLEANONTDEP=${CLEANONTDEP}:${STATSDONE}
CLEANPBDEP=${CLEANPBDEP}:${STATSDONE}

##############################################################################
## CLEAN UP -- THIS CLEAN UP ONLY FOR BAM USED FOR SNIFFLES AND ALE
##############################################################################
if $CLEAN; then
 DELFILE=$PBBAM
 #####if $SNIFFLESPB || $SNiFFLESCOMBINED; then
 if $SNIFFLESPB && $ALEPB; then
  GATEFILE1="sniffles_pb/${BASE}_pacbio.bedpe"
  GATEFILE2="ale_pb/${BASE}.ALE.txt"
  CLEANPBDONE=`sbatch -J ${BASE}_clean_pb --dependency=${CLEANPBDEP} -o ${OUT}/cleanPB.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE1=${GATEFILE1},GATEFILE2=${GATEFILE2},DELFILE=${DELFILE} ${SCRIPTS}/clean.2.sh | awk '{print $4}'`
 elif $SNIFFLESPB || $ALEPB; then
  if $SNIFFLESPB; then GATEFILE="sniffles_pb/${BASE}_pacbio.bedpe";
    else GATEFILE="ale_pb/${BASE}.ALE.txt"; fi
  CLEANPBDONE=`sbatch -J ${BASE}_clean_pb --dependency=${CLEANPBDEP} -o ${OUT}/cleanPB.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`
 fi

 #####if $SNIFFLESONT || $SNIFFLESCOMBINED; then
 DELFILE=$ONTBAM
 if $SNIFFLESONT && $ALEONT; then
  GATEFILE1="sniffles_ont/${BASE}_ont.bedpe"
  GATEFILE2="ale_ont/${BASE}.ALE.txt"
  CLEANPBDONE=`sbatch -J ${BASE}_clean_ont --dependency=${CLEANONTDEP} -o ${OUT}/cleanONT.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE1=${GATEFILE1},GATEFILE2=${GATEFILE2},DELFILE=${DELFILE} ${SCRIPTS}/clean.2.sh | awk '{print $4}'`
 elif $SNIFFLESONT || $ALEONT; then
  if $SNIFFLESONT; then GATEFILE="sniffles_ont/${BASE}_ont.bedpe";
    else GATEFILE="ale_ont/${BASE}.ALE.txt"; fi
  CLEANPBDONE=`sbatch -J ${BASE}_clean_ont --dependency=${CLEANONTDEP} -o ${OUT}/cleanONT.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`
 fi

 if $SNIFFLESCOMBINED; then
  GATEFILE="sniffles_combined/${BASE}_combined.bedpe"
  DELFILE=$COMBBAM
  CLEANCOMBINEDONE=`sbatch -J ${BASE}_clean_combined --dependency=afterok:${SNIFFLESCOMBDONE}:${STATSDONE} -o ${OUT}/cleanCombined.slurm.%A.out --mem=2g --time=1:00:00 -c 1 \
    --qos=$QOS --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'` 
 fi
fi




##############################################################################
## LAP SECTION ::: 
##############################################################################

##############################################################################
## MAP FOR LAP PACBIO
##############################################################################
D=mreads
if $MAPLAPPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $BUILDBWA; then
  MAPLAPPBDONE=`sbatch -J ${BASE}_mapPBreads_lap --dependency=afterok:${IDXDEP} -o ${OUT}/lap_mapPBreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=pacbio,BWAIDX=${BWAIDX},FASTQ=${PACBIO}  ${SCRIPTS}/bwa-map-for-lap.sh | awk '{print $4}'`
 else
  MAPLAPPBDONE=`sbatch -J ${BASE}_mapPBreads_lap -o ${OUT}/lap_mapPBreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=pacbio,BWAIDX=${BWAIDX},FASTQ=${PACBIO}  ${SCRIPTS}/bwa-map-for-lap.sh | awk '{print $4}'`
 fi
 cd ../
 LAPCLEANPBDEP=${MAPLAPPBDONE}
fi
PBSAMLAP=`readlink -f ${MAIN}/${D}/pacbio.lap.sam`


##############################################################################
## MAP FOR LAP ONT
##############################################################################
D=mreads
TYPE=ont2d
if $MAPLAPONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $BUILDBWA; then
  MAPLAPONTDONE=`sbatch -J ${BASE}_mapONTreads_lap --dependency=afterok:${IDXDEP} -o ${OUT}/lap_mapONTreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=${TYPE},BWAIDX=${BWAIDX},FASTQ=${ONT}  ${SCRIPTS}/bwa-map-for-lap.sh | awk '{print $4}'`
 else
  MAPLAPONTDONE=`sbatch -J ${BASE}_mapONTreads_lap -o ${OUT}/lap_mapONTreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=${TYPE},BWAIDX=${BWAIDX},FASTQ=${ONT}  ${SCRIPTS}/bwa-map-for-lap.sh | awk '{print $4}'`
 fi
 cd ../
 LAPCLEANONTDEP=${MAPLAPONTDONE}
fi
ONTSAMLAP=`readlink -f ${MAIN}/${D}/ont2d.lap.sam`


##############################################################################
## LAP PACBIO
##############################################################################
D=lap_pb
TYPE=pacbio
if $LAPPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPLAPPB; then
  LAPPBDONE=`sbatch -J ${BASE}_lap_pb --dependency=afterok:${MAPLAPPBDONE} -o ${OUT}/lap_pb.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${PBSAMLAP},REF=${ASM},MISMATCH=${PBMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.longread.sh | awk '{print $4}'`
 else
  LAPPBDONE=`sbatch -J ${BASE}_lap_pb -o ${OUT}/lap_pb.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${PBSAMLAP},REF=${ASM},MISMATCH=${PBMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.longread.sh | awk '{print $4}'`
 fi
 cd ../
 LAPCLEANPBDEP=${LAPCLEANPBDEP}:${LAPPBDONE}
fi



##############################################################################
## LAP ONT
##############################################################################
D=lap_ont
TYPE=ont2d
if $LAPONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPLAPONT; then
  LAPONTDONE=`sbatch -J ${BASE}_lap_ont --dependency=afterok:${MAPLAPONTDONE} -o ${OUT}/lap_ont.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${ONTSAMLAP},REF=${ASM},MISMATCH=${ONTMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.longread.sh | awk '{print $4}'`
 else
  LAPONTDONE=`sbatch -J ${BASE}_lap_ont -o ${OUT}/lap_ont.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${ONTSAMLAP},REF=${ASM},MISMATCH=${ONTMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.longread.sh | awk '{print $4}'`
 fi
 cd ../
 LAPCLEANONTDEP=${LAPCLEANONTDEP}:${LAPONTDONE}
fi


##############################################################################
## LAP CLEAN 
##############################################################################
if $CLEAN; then
 if $MAPLAPPB ; then
  GATEFILE="lap_pb/${BASE}.lapscore"
  DELFILE=$PBSAMLAP
  CLEANPBDONE=`sbatch -J ${BASE}_clean_pb_lap --dependency=afterok:${LAPCLEANPBDEP} -o ${OUT}/cleanPB_lap.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`
 fi
 if $MAPLAPONT ; then
  GATEFILE="lap_ont/${BASE}.lapscore"
  DELFILE=$ONTSAMLAP
  CLEANONTDONE=`sbatch -J ${BASE}_clean_ont_lap --dependency=afterok:${LAPCLEANONTDEP} -o ${OUT}/cleanONT_lap.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`
 fi
fi







##############################################################################
## LONG2PE SECTION ::: 
##############################################################################

##############################################################################
## MAP FOR LONG2PE PACBIO
##############################################################################
D=mreads
TYPE=pacbio
if $MAPPB_PE; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $BUILDBWA; then
  MAPL2PEPBDONE=`sbatch -J ${BASE}_mapPBreads_long2pe --dependency=afterok:${IDXDEP} -o ${OUT}/long2pe_mapPBreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=${TYPE},BWAIDX=${BWAIDX},R1=${PACBIO1},R2=${PACBIO2}  ${SCRIPTS}/bwa-map-for-long2pe.sh | awk '{print $4}'`
 else
  MAPL2PEPBDONE=`sbatch -J ${BASE}_mapPBreads_long2pe -o ${OUT}/long2pe_mapPBreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=${TYPE},BWAIDX=${BWAIDX},R1=${PACBIO1},R2=${PACBIO2}  ${SCRIPTS}/bwa-map-for-long2pe.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANPBDEP=${MAPL2PEPBDONE}
fi
PBBAML2PE=`readlink -f ${MAIN}/${D}/pacbio.long2pe.bam`

##############################################################################
## MAP FOR LONG2PE ONT
##############################################################################
D=mreads
TYPE=ont2d
if $MAPONT_PE; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $BUILDBWA; then
  MAPL2PEONTDONE=`sbatch -J ${BASE}_mapONTreads_long2pe --dependency=afterok:${IDXDEP} -o ${OUT}/long2pe_mapONTreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=${TYPE},BWAIDX=${BWAIDX},R1=${ONT1},R2=${ONT2}  ${SCRIPTS}/bwa-map-for-long2pe.sh | awk '{print $4}'`
 else
  MAPL2PEONTDONE=`sbatch -J ${BASE}_mapONTreads_long2pe -o ${OUT}/long2pe_mapONTreads.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=MTHREADS=${MTHREADS},TYPE=${TYPE},BWAIDX=${BWAIDX},R1=${ONT1},R2=${ONT2}  ${SCRIPTS}/bwa-map-for-long2pe.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANONTDEP=${MAPL2PEONTDONE}
fi
ONTBAML2PE=`readlink -f ${MAIN}/${D}/ont2d.long2pe.bam`






##############################################################################
## FRC PACBIO
##############################################################################
D=frc_pb
if $FRCPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPPB_PE; then
  FRCPBDONE=`sbatch -J ${BASE}_frc_pb --dependency=afterok:${MAPL2PEPBDONE} -o ${OUT}/frc_pb.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=BAM=${PBBAML2PE},L2PE_MAXINS=${L2PE_MAXINS},GSIZE=${GSIZE},BASE=${BASE} ${SCRIPTS}/frc.eval.sh | awk '{print $4}'`
 else
  FRCPBDONE=`sbatch -J ${BASE}_frc_pb -o ${OUT}/frc_pb.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=BAM=${PBBAML2PE},L2PE_MAXINS=${L2PE_MAXINS},GSIZE=${GSIZE},BASE=${BASE} ${SCRIPTS}/frc.eval.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANPBDEP=${L2PECLEANPBDEP}:${FRCPBDONE}
fi




##############################################################################
## FRC ONT
##############################################################################
D=frc_ont
if $FRCONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPONT_PE; then
  FRCONTDONE=`sbatch -J ${BASE}_frc_ont --dependency=afterok:${MAPL2PEONTDONE} -o ${OUT}/frc_ont.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=BAM=${ONTBAML2PE},L2PE_MAXINS=${L2PE_MAXINS},GSIZE=${GSIZE},BASE=${BASE} ${SCRIPTS}/frc.eval.sh | awk '{print $4}'`
 else
  FRCONTDONE=`sbatch -J ${BASE}_frc_ont -o ${OUT}/frc_ont.slurm.%A.out --mem=$FMEM --time=$FTIME -c $FTHREADS --qos=$QOS --export=BAM=${ONTBAML2PE},L2PE_MAXINS=${L2PE_MAXINS},GSIZE=${GSIZE},BASE=${BASE} ${SCRIPTS}/frc.eval.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANONTDEP=${L2PECLEANONTDEP}:${FRCONTDONE}
fi






##############################################################################
## REAPR PACBIO
##############################################################################
D=reapr_pb
if $REAPRPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPPB_PE; then
  REAPRPBDONE=`sbatch -J ${BASE}_reapr_pb --dependency=afterok:${MAPL2PEPBDONE} -o ${OUT}/reapr_pb.slurm.%A.out --mem=$RMEM --time=$RTIME -c $RTHREADS --qos=$QOS --export=REF=${ASM},BASE=${BASE},AGGRESSIVE=${AGGRESSIVE_LONG},BAM=${PBBAML2PE},WINLEN=${WINLEN_PACBIO},MININNER=${MININNER_PACBIO},FCDWINLEN=${FCDWINLEN_PACBIO},MAXFCDSAMPLE=${MAXFCDSAMPLE} ${SCRIPTS}/reapr.eval.sh | awk '{print $4}'`
 else
  REAPRPBDONE=`sbatch -J ${BASE}_reapr_pb -o ${OUT}/reapr_pb.slurm.%A.out --mem=$RMEM --time=$RTIME -c $RTHREADS --qos=$QOS --export=REF=${ASM},BASE=${BASE},AGGRESSIVE=${AGGRESSIVE_LONG},BAM=${PBBAML2PE},WINLEN=${WINLEN_PACBIO},MININNER=${MININNER_PACBIO},FCDWINLEN=${FCDWINLEN_PACBIO},MAXFCDSAMPLE=${MAXFCDSAMPLE} ${SCRIPTS}/reapr.eval.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANPBDEP=${L2PECLEANPBDEP}:${REAPRPBDONE}
fi




##############################################################################
## REAPR ONT
##############################################################################
D=reapr_ont
if $REAPRONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPONT_PE; then
  REAPRONTDONE=`sbatch -J ${BASE}_reapr_ont --dependency=afterok:${MAPL2PEONTDONE} -o ${OUT}/reapr_ont.slurm.%A.out --mem=$RMEM --time=$RTIME -c $RTHREADS --qos=$QOS --export=REF=${ASM},BASE=${BASE},AGGRESSIVE=${AGGRESSIVE_LONG},BAM=${ONTBAML2PE},WINLEN=${WINLEN_ONT},MININNER=${MININNER_ONT},FCDWINLEN=${FCDWINLEN_ONT},MAXFCDSAMPLE=${MAXFCDSAMPLE} ${SCRIPTS}/reapr.eval.sh | awk '{print $4}'`
 else
  REAPRONTDONE=`sbatch -J ${BASE}_reapr_ont -o ${OUT}/reapr_ont.slurm.%A.out --mem=$RMEM --time=$RTIME -c $RTHREADS --qos=$QOS --export=REF=${ASM},BASE=${BASE},AGGRESSIVE=${AGGRESSIVE_LONG},BAM=${ONTBAML2PE},WINLEN=${WINLEN_ONT},MININNER=${MININNER_ONT},FCDWINLEN=${FCDWINLEN_ONT},MAXFCDSAMPLE=${MAXFCDSAMPLE} ${SCRIPTS}/reapr.eval.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANONTDEP=${L2PECLEANONTDEP}:${REAPRONTDONE}
fi





##############################################################################
## LAP Long2pe section -- NOTE: for now (maybe always), this just uses the "long2pe" reads like "se" reads - same as LAP above.
## One reason to do so is that the way LAP is calculated goes to 0 within a few kb for long reads...
## The reads used for long2pe essentially break long reads up into 1kb sections, which allows many more scores >0 to be used to discriminate assemblies
## Shortcoming here is that not all mappings, not even more than 1 (I dont think) are reported here -- so that aspect does not satisfy LAP assumptions.
##############################################################################

##############################################################################
## LAP PACBIO Long2pe
##############################################################################
D=lap_pb_long2pe
if $LAPPB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPPB_PE; then
  LAPPBL2PEDONE=`sbatch -J ${BASE}_lap_pb_long2pe --dependency=afterok:${MAPL2PEPBDONE} -o ${OUT}/lap_pb_long2pe.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${PBBAML2PE},REF=${ASM},MISMATCH=${PBMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.long2pe.sh | awk '{print $4}'`
 else
  LAPPBL2PEDONE=`sbatch -J ${BASE}_lap_pb_long2pe -o ${OUT}/lap_pb_long2pe.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${PBBAML2PE},REF=${ASM},MISMATCH=${PBMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.long2pe.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANPBDEP=${L2PECLEANPBDEP}:${LAPPBL2PEDONE}
fi




##############################################################################
## LAP ONT Long2pe
##############################################################################
D=lap_ont_long2pe
if $LAPONT; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 if $MAPONT_PE; then
  LAPONTL2PEDONE=`sbatch -J ${BASE}_lap_ont_long2pe --dependency=afterok:${MAPL2PEONTDONE} -o ${OUT}/lap_ont_long2pe.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${ONTBAML2PE},REF=${ASM},MISMATCH=${ONTMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.long2pe.sh | awk '{print $4}'`
 else
  LAPONTL2PEDONE=`sbatch -J ${BASE}_lap_ont_long2pe -o ${OUT}/lap_ont_long2pe.slurm.%A.out --mem=$LMEM --time=$LTIME -c $LTHREADS --qos=$QOS --export=BASE=${BASE},SAM=${ONTBAML2PE},REF=${ASM},MISMATCH=${ONTMISMATCHRATE},P=${LTHREADS} ${SCRIPTS}/lap.long2pe.sh | awk '{print $4}'`
 fi
 cd ../
 L2PECLEANONTDEP=${L2PECLEANONTDEP}:${LAPONTL2PEDONE}
fi




##############################################################################
## LONG2PE CLEAN 
##############################################################################
if $CLEAN; then

 DELFILE=$PBBAML2PE
 if $FRCPB && $REAPRPB ; then
  GATEFILE1="frc_pb/${BASE}.long2peFeatures.gff"
  GATEFILE2="reapr_pb/output_directory/05.summary.report.txt"
  CLEANPBL2PE=`sbatch -J ${BASE}_clean_pb_long2pe --dependency=afterok:${L2PECLEANPBDEP} -o ${OUT}/cleanPB_l2pe.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE1=${GATEFILE1},GATEFILE2=${GATEFILE2},DELFILE=${DELFILE} ${SCRIPTS}/clean.2.sh | awk '{print $4}'`
 elif $FRCPB || $REAPRPB ; then
  if $FRCPB; then GATEFILE="frc_pb/${BASE}.long2peFeatures.gff"; 
    else GATEFILE="reapr_pb/output_directory/05.summary.report.txt"; fi
  CLEANPBL2PE=`sbatch -J ${BASE}_clean_pb_long2pe --dependency=afterok:${L2PECLEANPBDEP} -o ${OUT}/cleanPB_l2pe.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`  
 fi

 DELFILE=$ONTBAML2PE
 if $FRCONT && $REAPRONT ; then
  GATEFILE1="frc_ont/${BASE}.long2peFeatures.gff"
  GATEFILE2="reapr_ont/output_directory/05.summary.report.txt"
  DELFILE=$ONTBAML2PE
  CLEANONTL2PE=`sbatch -J ${BASE}_clean_ont_long2pe1 --dependency=afterok:${L2PECLEANONTDEP} -o ${OUT}/cleanONT_l2pe.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE1=${GATEFILE1},GATEFILE2=${GATEFILE2},DELFILE=${DELFILE} ${SCRIPTS}/clean.2.sh | awk '{print $4}'`
 elif $FRCONT || $REAPRONT ; then
  if $FRCONT; then GATEFILE="frc_ont/${BASE}.long2peFeatures.gff"; 
    else GATEFILE="reapr_ont/output_directory/05.summary.report.txt"; fi
  CLEANONTL2PE=`sbatch -J ${BASE}_clean_ont_long2pe --dependency=afterok:${L2PECLEANONTDEP} -o ${OUT}/cleanONT_l2pe.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
    --export=GATEFILE=${GATEFILE},DELFILE=${DELFILE} ${SCRIPTS}/clean.sh | awk '{print $4}'`  
 fi

fi
