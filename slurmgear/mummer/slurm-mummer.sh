#!/bin/bash

## Removed args
##        -c with argument = path to CONFIG file

function help {
    echo "
        Usage: ${0} -r:n:s:m:a:q:x:M:T:C:l:p:b:c:g:i:L:h
        Currently, only required arg is -r READS.
        -r with argument = abs path to fasta/fastq file of query fasta or path from HOME; do not give relative path from pwd unless it is in pwd or subdir.
        -n with argument = NJOBS -- how many jobs to split query file into. (Default: 1)
        -s with argument = path to SCRIPTS DIR (Default: same directory as this script).
        -m with argument = path to MUMMER DIR (Default: ~/software/quickmerge/quickmerge/MUMmer3.23/).
        -a with argument = path to ASM FOFN (Default: input.fofn)
        -q with argument = Primary QOS for sbatch. (Default: epscor-condo)
        -x with argument = Secondary QOS for sbatch. (Default: biomed-sb-condo)
        -I with argument = Higher numbers skew this toward using primary QOS more than secondary. Setting to 2 would be even split. (Default: 9)
        -M with argument = how much memory to tell sbatch. (Default: 16g)
        -T with argument = how much time to tell sbatch. (Default: 24:00:00)
        -C with argument = how many cpus/threads to tell sbatch. (Default: 1)
        -l with argument = minMatchLength for nucmer (Default: 20)
        -p with argument = PREFIX for nucmer (Default: out)
        -b with argument = breakLen for nucmer (Default: 200)
        -c with argument = mincluster for nucmer (Default: 65)
        -g with argument = maxgap for nucmer (Default: 90)
        -i with argument = identity for delta-filter (Default: 0)
        -L with argument = minAlnLen for delta-filter (Default: 0)
        -h help - returns this message; also returns this when no arguments given
"
}

if [ $# -eq 0 ]; then help; exit; fi

## DEFAULTS
READSDIR=readsdir
CLEAN=false
ASMFOFN=input.fofn
HELP=false
JTHREADS=1
JMEM=16g
JTIME=24:00:00
IMAX=9
QOS1=epscor-condo
QOS2=biomed-sb-condo
SCRIPTS=`abspath.py ${0} --split | awk '{print $1}'`
MUMMER=~/software/quickmerge/quickmerge/MUMmer3.23/
NJOBS=1
minMatchLength=20
PREFIX=out
breakLen=200
mincluster=65
maxgap=90
identity=0
minAlnLen=0


## Currently these defaults do not have corresponding options for changing
SLURMOUTDIR=slurmout

EXIT=false

while getopts "r:n:s:m:a:q:x:I:M:T:C:l:p:b:c:g:i:L:h" arg; do
    case $arg in
        r) READSFASTA=$OPTARG;;
        n) NJOBS=$OPTARG;;
        s) SCRIPTS=$OPTARG;;
        s) MUMMER=$OPTARG;;
        a) ASMFOFN=$OPTARG;;
        q) QOS1=$OPTARG;;
        x) QOS2=$OPTARG;;
        I) IMAX=$OPTARG;;
        M) JMEM=$OPTARG;;
        T) JTIME=$OPTARG;;
        C) JTHREADS=$OPTARG;;
        l) minMatchLength=$OPTARG;;
        p) PREFIX=$OPTARG;;
        b) breakLen=$OPTARG;;
        c) mincluster=$OPTARG;;
        g) maxgap=$OPTARG;;
        i) identity=$OPTARG;;
        L) minAlnLen=$OPTARG;;
        h) HELP=true;;
        *) help; exit;;
    esac
done


if ${HELP}; then help; exit; fi


READSDIR=`readlink -f $TMP`


## SPLIT UP READS FASTA INTO MULTIPLE FILES
if [[ "$READSFASTA" == *.fa ]]; then PRE=`basename $READSFASTA .fa`;
elif [[ "$READSFASTA" == *.fasta ]]; then PRE=`basename $READSFASTA .fasta`; fi
facount=`grep -c ">" $READSFASTA`
count=`echo $facount+$NJOBS | bc`
nlines=`echo $count/$NJOBS | bc`
splitFastA.py -f $READSFASTA -n $nlines -o $READSDIR

## SUBMIT BATCH JOBS
i=0
while read REF; do
  i=$(( $i+1 ))
  if [ $i -eq $IMAX ]; then QOS=${QOS2}; i=0; else QOS=${QOS1}; fi
  if [[ "$REF" == *.fasta ]]; then B=`basename $REF .fasta`; fi
  if [[ "$REF" == *.fa ]]; then B=`basename $REF .fa`; fi
  echo $B; 
  if [ ! -d $B ]; then mkdir $B; fi
  cd $B;
    MAIN=$PWD
    if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
    OUT=${MAIN}/${SLURMOUTDIR}
    cp ${REF} input_assembly.fasta
    ASM=`readlink -f ${REF}`
    EXPORTS=`echo PROTOCOL=READSDIR=${READSDIR},PRE=${PRE},NJOBS=${NJOBS},minMatchLength=${minMatchLength},PREFIX=${PREFIX},breakLen=${breakLen},mincluster=${mincluster},maxgap=${maxgap},identity=${identity},minAlnLen=${minAlnLen}`
    MUMMERDONE=`sbatch --dependency=afterok:${MAKEDONE} -a 1-$NJOBS -J ${BASE}_mummer -o ${OUT}/mummer.slurm.%A_%a.out --mem=$MEM --time=$TIME -c $THREADS --qos=$QOS \
      --export=${EXPORTS} ${SCRIPTS}/run-mummer.sh | awk '{print $4}'`
  cd ../
done < $ASMFOFN



