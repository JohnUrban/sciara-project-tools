#!/bin/bash

## Removed args
##        -c with argument = path to CONFIG file

function help {
    echo "
        Usage: ${0} -r:s:a:p:m:t:n:q:x:M:T:C:SGQCh23456
        Currently, only required arg is -r READS.
        -r with argument = abs path to fastq file of long reads or path from HOME; do not give relative path from pwd unless it is in pwd or subdir.
        -s with argument = path to SCRIPTS DIR
        -a with argument = path to ASM FOFN (Default: input.fofn)
        -p with argument = minPctIdentity for BLASR; (Default: 75)
        -m with argument = minMatch for BLASR; (Default: 8)
        -t with argument = sdpTupleSize for BLASR; (Default: 8)
        -n with argument = nCandidates for BLASR; (Default: 10)
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
        -2 start from step 2: mapping
        -3 start from step 3: support
        -4 start from step 4: extraction
        -5 start from step 5: assembly
        -6 start from step 6: output
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
minMatch=8
sdpTupleSize=8
nCandidates=10
MAKEFAKEQUALS=false
GAPSONLY=false
SPANONLY=false
RUNSETUP=true
RUNMAPPING=true
RUNSUPPORT=true
RUNSUPPORT=true
RUNEXTRACTION=true
RUNASSEMBLY=true
RUNOUTPUT=true

## Currently these defaults do not have corresponding options for changing
bestn=1
maxScore=-500
nproc=48
SLURMOUTDIR=slurmout

EXIT=false

#### OPTIONS AND COMMANDLINE ARGS
##        c) CONFIG=$OPTARG;;

while getopts "r:s:a:p:m:t:n:q:x:M:T:C:SGQCh23456" arg; do
    case $arg in
        r) READS=$OPTARG;;
        s) SCRIPTS=$OPTARG;;
        a) ASMFOFN=$OPTARG;;
        p) minPctIdentity=$OPTARG;;
        m) minMatch=$OPTARG;;
        t) sdpTupleSize=$OPTARG;;
        n) nCandidates=$OPTARG;;
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
        2) RUNSETUP=false;;
        3) RUNSETUP=false; RUNMAPPING=false;;
        4) RUNSETUP=false; RUNMAPPING=false; RUNSUPPORT=false;;
        5) RUNSETUP=false; RUNMAPPING=false; RUNSUPPORT=false; RUNEXTRACTION=false;;
        6) RUNSETUP=false; RUNMAPPING=false; RUNSUPPORT=false; RUNEXTRACTION=false; RUNASSEMBLY=false;;
        *) help; exit;;
    esac
done


if ${HELP}; then help; exit; fi

##source $CONFIG
##PIPELINE=${SCRIPTS}/pbjelly-pipeline.sh


## CHECK FILENAMES TO MAKE SURE THEY ARE NAMED CORRECTLY
while read REF; do 
  if [[ "$REF" != *.fasta ]]; then echo "All files need to have .fasta extension for PBJelly. Not .fa, etc."; exit; fi
done < $ASMFOFN



## SUBMIT BATCH JOBS
READS_BASEDIR=`abspath.py ${READS} --split | awk '{print $1}'`
READS_BASENAME=`basename ${READS}`
i=0
while read REF; do
  i=$(( $i+1 ))
  if [ $i -eq $IMAX ]; then QOS=${QOS2}; i=0; else QOS=${QOS1}; fi
  if [[ "$REF" == *.fasta ]]; then B=`basename $REF .fasta`; fi
  echo $B; 
  if [ ! -d $B ]; then mkdir $B; fi
  cd $B;
    ##$PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF
    MAIN=$PWD
    if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
    OUT=${MAIN}/${SLURMOUTDIR}
    ## make Protocol.xml
    cp ${REF} input_assembly.fasta
    ASM=`readlink -f input_assembly.fasta`
    ${SCRIPTS}/get-pbjelly-protocol.py -r ${ASM} -o ${PWD} -b ${READS_BASEDIR} -f ${READS_BASENAME} --minMatch ${minMatch} --sdpTupleSize ${sdpTupleSize} --minPctIdentity ${minPctIdentity} --bestn ${bestn} --nCandidates ${nCandidates} --maxScore ${maxScore} --nproc ${nproc} > Protocol.xml
    PROTOCOL=`abspath.py Protocol.xml`
    EXPORTS=`echo PROTOCOL=${PROTOCOL},MAKEFAKEQUALS=${MAKEFAKEQUALS},GAPSONLY=${GAPSONLY},SPANONLY=${SPANONLY},THREADS=${JTHREADS},ASM=${ASM},RUNSETUP=${RUNSETUP},RUNMAPPING=${RUNMAPPING},RUNSUPPORT=${RUNSUPPORT},RUNEXTRACTION=${RUNEXTRACTION},RUNASSEMBLY=${RUNASSEMBLY},RUNOUTPUT=${RUNOUTPUT}`
    PBJDONE=`sbatch -J ${B}_pbjelly -o ${OUT}/pbjelly.slurm.%A.out --mem=$JMEM --time=$JTIME -c $JTHREADS --qos=$QOS \
      --export=${EXPORTS} ${SCRIPTS}/ApplyJelly.sh | awk '{print $4}'`
  cd ../
done < $ASMFOFN



