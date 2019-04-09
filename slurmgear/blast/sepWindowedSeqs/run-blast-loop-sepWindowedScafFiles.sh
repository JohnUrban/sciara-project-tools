#!/bin/bash

ASM=/gpfs/scratch/jmurban/final_bacterial_ids/finalscafs/asms/01_canu_hybridscaffold_T8.quiver5.pbj3.quiver3.pbj3.quiver3.pilon12.fasta
NICKNAME=canu-blastn-wc
MINJOBS=1000
GSIZE=292000000
WSIZE=1000

NWIN=`echo ${GSIZE}/${WSIZE} | bc`
MAXLINESPERFILE=`echo ${NWIN}/${MINJOBS} | bc`

echo Given a minimum of ${MINJOBS}, there are
echo ${NWIN} windows...
echo a maximum of ${MAXLINESPERFILE} lines per file....

BLASTDIR=blast
SLURMOUT=slurmout
TIGDIR=contigs
P=8
MEM=30g
TIME=48:00:00
QOS=biomed-condo
E=1e-10

mkdir -p $BLASTDIR
mkdir -p $SLURMOUT
mkdir -p $TIGDIR

echo Separating contigs....
multiFasta2manyFastas.py -f ${ASM} -d ${TIGDIR}


mkdir -p windows
TMP=tmpdir
mkdir -p ${TMP}
QUERYDIR=`readlink -f $TMP`

echo Making windows for each contig.....
for fa in ${TIGDIR}/*fa; do
  ASM=`readlink -f ${fa}`
  B=`basename ${ASM} .fa`
  G=windows/${B}.size
  W=windows/${B}.bed
  WF=windows/${B}.fa
  GAP=windows/${B}_gap.bed
  echo Sizing... ${B}; faSize -detailed ${ASM} > ${G}
  echo Windowing... ${B}; bedtools makewindows -w ${WSIZE} -s ${WSIZE} -g ${G} > ${W}
  echo Gapping... ${B}
  scf-N-to-BED.py ${ASM} > ${GAP}
  echo Ungapping... ${B}
  intersectBed -f 0.75 -v -a ${W} -b ${GAP} > ${W}.1 ## to be considered an overlap at least 75% of the window needs to overlap, otherwise include it
  mv ${W}.1 ${W}
  echo Extracting... ${B}; fastaFromBed -fi ${ASM} -bed ${W} > ${WF}
  echo Calculting.... ${B}
  echo Splitting.... ${B}; splitFastA.py -f ${WF} -n ${MAXLINESPERFILE} -o $TMP
done

echo Looping....
for fa in ${QUERYDIR}/*fa; do
 b=`basename $fa .fa`
 echo Launching... ${b}
 sbatch --account=${QOS} --mem=$MEM --time=$TIME -c $P -o ${SLURMOUT}/${b}.%A.out -J ${b}_blast_NT_${NICKNAME} --export=Q=${fa},P=${P},BLASTDIR=${BLASTDIR},E=${E} blastn-sensitive.sh
done
