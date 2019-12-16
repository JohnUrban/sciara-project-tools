#!/bin/bash
###########################

if [ $# -eq 0 ]; then echo "
Arg1=/Path/To/Reference.fasta 
Arg2=QOS
Arg3=Logical(true/false) should directories be cleaned up when quiver is done?
Arg4=config file
Arg5=PVALUE (0.01 to use default)
Arg6=INPUT FOFN of BAX.h5 files
Arg7=IDENTIFY -- specify integer 1-7; 1=m6A, 2=m4C, 3=m5C_TET, 4=m6A,m4C, 5=m6A,m5C_TET, 6=m4C,m5C_TET, 7=m6A,m4C,m5C_TET
"; exit; fi

set > del.txt
sort del.txt -o del.sort.txt
rm del.txt
module load blasr/2015Oct22-8cc8621
## 1. Make a file called input.fofn with path to all .bax.h5 files
## 2. Fill in REF and ASMOUTPREFIX below; change INPUT if necessary
## 3. Make changes to other parameters if necessary
## 4. Run as  "bash quiverpipeline.parallel.sh"


REF=$1               ## path to pre-quiver asm to polish
QOS=$2
ALTQOS=$2
QQOS=$2
CLEAN=$3
CONFIG=$4
PVALUE=$5
INPUT=$6
IDENTIFY="$7"

source $CONFIG


ASMOUTPREFIX=`basename $REF .fasta`
OUTGFF=${ASMOUTPREFIX}.variants.gff
OUTFASTQ=${ASMOUTPREFIX}.fastq
OUTFASTA=${ASMOUTPREFIX}.fasta

if [ ! -d $MAIN_TMP_DIR ]; then mkdir $MAIN_TMP_DIR; fi
if [ ! -d $CMPDIR ]; then mkdir $CMPDIR; fi
if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi

##############################################################################
## ADD BUBBLES AND/OR RENAME (OPTIONAL)
##############################################################################
if $EXTRACT_CANU_UNASM; then
  filterCanuFasta.py --minreads 2 $CANU_UNASM > $EXTRACTED_CANU_UNASM
  BUBBLES_ETC="${BUBBLES_ETC} $EXTRACTED_CANU_UNASM"
fi
  

if $ADD_BUBBLES_ETC; then
 cat $REF > $COMBINED_CONTIGS_BUBBLES_ETC
 for F in $BUBBLES_ETC; do
  if [ ${F: -6} == ".fasta" ]; 
   then base=`basename $F .fasta`; 
  elif [ ${F: -3} == ".fa" ]; 
   then base=`basename $F .fa`; 
  else base=$F; 
  fi 
  fasta_name_changer.py -f $F --front remove_${base}_ >> $COMBINED_CONTIGS_BUBBLES_ETC
 done
 REF=$COMBINED_CONTIGS_BUBBLES_ETC
fi

if $RENAME; then
 fasta_name_changer.py -f $REF -k 1 > $RENAMED_CONTIGS
 REF=$RENAMED_CONTIGS
fi



##############################################################################
## FINAL VARIABLES
##############################################################################
#printenv > final_variables.txt
set > del.txt
sort del.txt -o del2.sort.txt
rm del.txt
comm -1 -3 del.sort.txt del2.sort.txt > final_variables.txt
rm del.sort.txt del2.sort.txt

##############################################################################
## PBALIGN JOBS
##############################################################################
if $PBALIGN; then
 DEPLIST=''
 i=0
 while read file; do
  i=$(( $i+1 ))
  b=`basename $file .bax.h5`
  TMPDIR=${MAIN_TMP_DIR}/${b}
  if [ ! -d $TMPDIR ]; then mkdir $TMPDIR; fi
  DEP=`sbatch -J ${ASMOUTPREFIX}_pbalign_${b} -o ${SLURMOUTDIR}/${b}.slurm.%A.out --mem=$AMEM --time=$ATIME -c $ATHREADS --qos=$QOS --export=THREADS=${ATHREADS},REF=$REF,OUTPRE=${CMPDIR}/${b},BAX=${file},TMPDIR=${TMPDIR} ${PATH_TO_SCRIPTS}/pbalign.sh | awk '{print $4}'`
  DEPLIST=$DEPLIST":"$DEP
  PREV=$QOS
  QOS=$ALTQOS 
  ALTQOS=$PREV
 done < $INPUT
fi


##############################################################################
## MERGE
##############################################################################

if $MERGE; then
 if $PBALIGN; then
  SORTDEP=`sbatch -J ${ASMOUTPREFIX}_cmp_merge --dependency=afterok${DEPLIST} -o ${SLURMOUTDIR}/cmp_merge.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=CMPDIR=${CMPDIR},MERGEDCMP=${MERGEDCMP}  ${PATH_TO_SCRIPTS}/cmp_merge.sh | awk '{print $4}'`
 else
  SORTDEP=`sbatch -J ${ASMOUTPREFIX}_cmp_merge -o ${SLURMOUTDIR}/cmp_merge.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --qos=$QOS --export=CMPDIR=${CMPDIR},MERGEDCMP=${MERGEDCMP}  ${PATH_TO_SCRIPTS}/cmp_merge.sh | awk '{print $4}'`
 fi
fi 


##############################################################################
## SORT
##############################################################################

if $SORT; then
 if [ ! -d $TMP_SORT ]; then mkdir $TMP_SORT; fi
 if $MERGE; then
  KINDEP=`sbatch -J ${ASMOUTPREFIX}_cmp_sort --dependency=afterok:${SORTDEP} -o ${SLURMOUTDIR}/cmp_sort.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=TMP_SORT=${TMP_SORT},IN_CMP=${MERGEDCMP} ${PATH_TO_SCRIPTS}/cmp_sort.sh | awk '{print $4}'`
 else
  KINDEP=`sbatch -J ${ASMOUTPREFIX}_cmp_sort -o ${SLURMOUTDIR}/cmp_sort.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=TMP_SORT=${TMP_SORT},IN_CMP=${MERGEDCMP} ${PATH_TO_SCRIPTS}/cmp_sort.sh | awk '{print $4}'`
 fi
fi

##############################################################################
## KINETICS TOOLS
##############################################################################
if $KINETICS; then
 if $SORT; then
  CLEANDEP=`sbatch -J ${ASMOUTPREFIX}_kinetics_ipdSummary --dependency=afterok:${KINDEP} -o ${SLURMOUTDIR}/kinetics.slurm.%A.out --mem=$QMEM --time=$QTIME -c $QTHREADS --qos=$QQOS --export=PVALUE=${PVALUE},IDENTIFY=${IDENTIFY},THREADS=${QTHREADS},INPUT=${MERGEDCMP},REF=${REF} ${PATH_TO_SCRIPTS}/kinetics.sh | awk '{print $4}'`
 else
  CLEANDEP=`sbatch -J ${ASMOUTPREFIX}_kinetics_ipdSummary -o ${SLURMOUTDIR}/kinetics.slurm.%A.out --mem=$QMEM --time=$QTIME -c $QTHREADS --qos=$QQOS --export=PVALUE=${PVALUE},IDENTIFY=${IDENTIFY},THREADS=${QTHREADS},INPUT=${MERGEDCMP},REF=${REF} ${PATH_TO_SCRIPTS}/kinetics.sh | awk '{print $4}'`
 fi
fi


##############################################################################
## CLEANUP -- assumes you would only require automated sbatch-dependency clean-up after successful kineticsTools.
##         -- otherwise, just create ad hoc script to clean up dirs manually...
##############################################################################
if $CLEAN && $KINETICS; then
  sbatch -J ${ASMOUTPREFIX}_cleanup --dependency=afterok:${CLEANDEP} -o ${SLURMOUTDIR}/cleanup.slurm.%A.out --mem=4g --time=2:00:00 -c 2 --qos=$QOS ${PATH_TO_SCRIPTS}/cleanup.sh | awk '{print $4}'
fi

