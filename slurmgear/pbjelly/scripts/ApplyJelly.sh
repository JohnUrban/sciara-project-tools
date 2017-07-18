#!/bin/bash
module load blasr

## ECHO OUT VARIABLES
for var in PROTOCOL MAKEFAKEQUALS GAPSONLY SPANONLY THREADS; do echo -e ${var} ${!var}; done

## DEFINE FUNCTIONS
function start {
echo Start ${1} >> JellyApply.log
date >> JellyApply.log
}

function end {
echo End ${1} >> JellyApply.log
date >> JellyApply.log
echo >> JellyApply.log
}

function run {
 #arg1 = step name
 #arg2 = path to protocol.xml
 #arg3 = quotes around any args passed to step -- eg. "-p 2 -x 5"
 step=${1}
 protocol=${2}
 args=${3}
 start ${1}
 Jelly.py ${step} ${protocol} ${args} 2>>JellyApply.log 
 end ${1}
}


touch JellyApply.log
if ${MAKEFAKEQUALS}; then
  BASE=`basename ${ASM} .fasta` 
  fakeQuals.py ${ASM} ${BASE}.qual
fi

## if using fastq, may want the following as a solution
## fastq2faqual.py --fastq /path/to/polished_assembly.fastq --fa --qual --out ${outpre}


## SOME ARGUMENTS
SETUPARGS=""
MAPARGS=""
SUPPORTARGS=""
EXTRACTARGS=""
ASSEMBLYARGS='-x "--nproc=${THREADS}"'
OUTPUTARGS=""

TMPSUPP=""
if ${GAPSONLY}; then TMPSUPP+="--capturedOnly "; fi
if ${SPANONLY}; then TMPSUPP+="--spanOnly "; fi
if ${GAPSONLY} || ${SPANONLY}; then SUPPORTARGS+='-x "${TMPSUPP}"'




## RUNNING
run setup ${PROTOCOL} ${SETUPARGS}
run mapping ${PROTOCOL} ${MAPARGS} 
run support ${PROTOCOL} ${SUPPORTARGS}
run extraction ${PROTOCOL} ${EXTRACTARGS}
run assembly ${PROTOCOL} ${ASSEMBLYARGS}
run output ${PROTOCOL} ${OUTPUTARGS}

