#!/bin/bash
###########################
module load blasr/2015Oct22-8cc8621
set > del.txt
sort del.txt -o del.sort.txt
rm del.txt
## 1. Make a file called input.fofn with path to all .bax.h5 files
## 2. Fill in REF and ASMOUTPREFIX below; change INPUT if necessary
## 3. Make changes to other parameters if necessary
## 4. Run as  "bash quiverpipeline.parallel.sh"

## CANU
## PARAMS IN THIS FILE BEST FOR CANU
## ONLY NEED TO PUT /PATH/TO/CANU/ASM ... up to the ".contigs.fasta" -- it will then fill in rest below for .bubbles.fasta etc
GENERAL_PRE=""
ASMOUTPREFIX=  ## e.g. canu.default.pball

INPUT=/gpfs_home/jurban/software/quiver/scripts/input.fofn #input.fofn
REF=${GENERAL_PRE}.contigs.fasta
BUBBLES_ETC="${GENERAL_PRE}.bubbles.fasta" ##e.g. canu.bubbles.fasta or falcon.a_ctg.fasta; Only add canu.unasm.fa file if it has been pre-filtered or you want to use all contents - else put it in CANU_UNASM variable below
CANU_UNASM="${GENERAL_PRE}.unassembled.fasta"  ##to filter for unasm seqs with >= 2 reads support (can adjust minreads below) -- the resulting file is automatically added to BUBBLES_ETC


PATH_TO_SCRIPTS=/gpfs_home/jurban/software/quiver/scripts
EXTRACT_CANU_UNASM=true
ADD_BUBBLES_ETC=true
RENAME=true
PBALIGN=true
MERGE=true
SORT=true
QUIVER=true

MERGEDCMP=out_all.cmp.h5
MAIN_TMP_DIR=temp
TMP_SORT=${MAIN_TMP_DIR}/cmp_sort
ATHREADS=8
AMEM=24g
ATIME=06:00:00
MTHREADS=8
MMEM=32g
MTIME=12:00:00
STHREADS=8
SMEM=16g
STIME=12:00:00
QTHREADS=8
QMEM=64g
QTIME=48:00:00
OUTGFF=${ASMOUTPREFIX}.quiver.variants.gff
OUTFASTQ=${ASMOUTPREFIX}.quiver.fastq
OUTFASTA=${ASMOUTPREFIX}.quiver.fasta
TIME=48:00:00
QOS=biomed-sb-condo
ALTQOS=${QOS} #ccmb-condo
QQOS=biomed-sb-condo
CMPDIR=cmpfiles
SLURMOUTDIR=slurmout
COMBINED_CONTIGS_BUBBLES_ETC=input_contigs_bubbles_etc.fasta
RENAMED_CONTIGS=input_sequences_renamed.fasta
CANU_MINREADS=2 ## Feel free to adjust this to be higher (any lower than 2 and it will take all seqs - i.e. no filtering necessary, ust add to BUBBLES_ETC)
EXTRACTED_CANU_UNASM=canu_unasm_seqs_with_ge${CANU_MINREADS}_reads.fasta


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
  QUIVDEP=`sbatch -J ${ASMOUTPREFIX}_cmp_sort --dependency=afterok:${SORTDEP} -o ${SLURMOUTDIR}/cmp_sort.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=TMP_SORT=${TMP_SORT},IN_CMP=${MERGEDCMP} ${PATH_TO_SCRIPTS}/cmp_sort.sh | awk '{print $4}'`
 else
  QUIVDEP=`sbatch -J ${ASMOUTPREFIX}_cmp_sort -o ${SLURMOUTDIR}/cmp_sort.slurm.%A.out --mem=$SMEM --time=$STIME -c $STHREADS --qos=$QOS --export=TMP_SORT=${TMP_SORT},IN_CMP=${MERGEDCMP} ${PATH_TO_SCRIPTS}/cmp_sort.sh | awk '{print $4}'`
 fi
fi

##############################################################################
## QUIVER
##############################################################################
if $QUIVER; then
 if $SORT; then
  sbatch -J ${ASMOUTPREFIX}_quiver --dependency=afterok:${QUIVDEP} -o ${SLURMOUTDIR}/quiver.slurm.%A.out --mem=$QMEM --time=$QTIME -c $QTHREADS --qos=$QQOS --export=THREADS=${QTHREADS},INPUT=${MERGEDCMP},REF=${REF},OUTGFF=${OUTGFF},OUTFASTQ=${OUTFASTQ},OUTFASTA=${OUTFASTA} ${PATH_TO_SCRIPTS}/quiver.sh | awk '{print $4}'
 else
  sbatch -J ${ASMOUTPREFIX}_quiver -o ${SLURMOUTDIR}/quiver.slurm.%A.out --mem=$QMEM --time=$QTIME -c $QTHREADS --qos=$QQOS --export=THREADS=${QTHREADS},INPUT=${MERGEDCMP},REF=${REF},OUTGFF=${OUTGFF},OUTFASTQ=${OUTFASTQ},OUTFASTA=${OUTFASTA} ${PATH_TO_SCRIPTS}/quiver.sh | awk '{print $4}'
 fi
fi
