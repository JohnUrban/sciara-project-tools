#!/bin/bash

## Removed args
##        -c with argument = path to CONFIG file

function help {
    echo "
        Usage: ${0} -r:s:a:m:q:x:M:T:C:SGQCh
        Currently, only required arg is -r READS.
        -r with argument = path to fastq file of long reads
        -s with argument = path to SCRIPTS DIR
        -a with argument = path to ASM FOFN (Default: input.fofn)
        -m with argument = minPctIdentity for BLASR; (Default: 75)
        -q with argument = Primary QOS for sbatch. (Default: epscor-condo)
        -x with argument = Secondary QOS for sbatch. (Default: biomed-sb-condo)
        -I with argument = Higher numbers skew this toward using primary QOS more than secondary. Setting to 2 would be even split. (Default: 9)
        -M with argument = how much memory to tell sbatch. (Default: 60g)
        -T with argument = how much time to tell sbatch. (Default: 24:00:00)
        -C with argument = how many cpus/threads to tell sbatch. (Default: 48)
        -S Span only. Only use reads that span entire gap. (Default: False)
        -G Fill in gaps only. Do not try to scaffold. (Default: False)
        -Q Tells it to make fake quals from assembly fastas. Usually this is needed unless the assembly is fastq.
        -C Tells it to clean up when done. 
        -h help - returns this message; also returns this when no arguments given
"
}

if [ $# -eq 0 ]; then help; exit; fi

## DEFAULTS
CLEAN=false
ASMFOFN=input.fofn
MAKEFAKEQUALS=false
SPANONLY=false
GAPSONLY=false
HELP=false
JTHREADS=48
JMEM=60g
JTIME=24:00:00
IMAX=9
QOS1=epscor-condo
QOS2=biomed-sb-condo
SCRIPTS=`abspath.py ${0} --split | awk '{print $1}'`
minPctIdentity=75
MAKEFAKEQUALS=false
GAPSONLY=false
SPANONLY=false


## Currently these defaults do not have corresponding options for changing
minMatch=8
sdpTupleSize=8
bestn=1
nCandidates=10
maxScore=-500
nproc=48
SLURMOUTDIR=slurmout


#### OPTIONS AND COMMANDLINE ARGS
##        c) CONFIG=$OPTARG;;

while getopts "r:s:a:m:q:x:M:T:C:SGQCh" arg; do
    case $arg in
        r) READS=$OPTARG;;
        s) SCRIPTS=$OPTARG;;
        a) ASMFOFN=$OPTARG;;
        m) minPctIdentity=$OPTARG;;
        q) QOS1=$OPTARG;;
        x) QOS2=$OPTARG;;
        I) IMAX=$OPTARG;;
        M) JMEM=$OPTARG;;
        T) JTIME=$OPTARG;;
        C) JTHREADS=$OPTARG;;
        S) SPANONLY=true;;
        G) GAPSONLY=true;;
        Q) MAKEFAKEQUALS=true;;
        C) CLEAN=true;;
        h) HELP=true;;
        *) help; exit;;
    esac
done

if ${HELP}; then help; exit; fi

##source $CONFIG
##PIPELINE=${SCRIPTS}/pbjelly-pipeline.sh


i=0
while read REF; do
  i=$(( $i+1 ))
  if [ $i -eq $IMAX ]; then QOS=${QOS2}; i=0; else QOS=${QOS1}; fi
  ASM=`readlink -f $REF` 
  if [[ "$ASM" == *.fa ]]; then B=`basename $ASM .fa`; 
  elif [[ "$ASM" == *.fasta ]]; then B=`basename $ASM .fasta`; fi
  echo $B; 
  if [ ! -d $B ]; then mkdir $B; fi
  cd $B;
    ##$PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF
    MAIN=$PWD
    if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
    OUT=${MAIN}/${SLURMOUTDIR}
    ## make Protocol.xml
    READS_BASEDIR=`abspath.py ${READS} --split | awk '{print $1}'`
    READS_BASENAME=`basename ${READS}`
    ${SCRIPTS}/get-pbjelly-protocol.py -r ${ASM} -o ${PWD} -b ${READS_BASEDIR} -f ${READS_BASENAME} --minMatch ${minMatch} --sdpTupleSize ${sdpTupleSize} --minPctIdentity ${minPctIdentity} --bestn ${bestn} --nCandidates ${nCandidates} --maxScore ${maxScore} --nproc ${nproc} > Protocol.xml
    PROTOCOL=`abspath.py Protocol.xml`
    EXPORTS='PROTOCOL=${PROTOCOL},MAKEFAKEQUALS=${MAKEFAKEQUALS},GAPSONLY=${GAPSONLY},SPANONLY=${SPANONLY},THREADS=${JTHREADS}'
    PBJDONE=`sbatch -J ${B}_pbjelly -o ${OUT}/pbjelly.slurm.%A.out --mem=$JMEM --time=$JTIME -c $JTHREADS --qos=$QOS \
      --export=${EXPORTS} ${SCRIPTS}/ApplyJelly.sh | awk '{print $4}'`
  cd ../
done < $ASMFOFN



