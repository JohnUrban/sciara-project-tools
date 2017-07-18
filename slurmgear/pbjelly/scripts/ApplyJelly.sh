#!/bin/bash

echo inside script now $0

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
  fakeQuals.py $2.fasta $2.qual
fi

## if using fastq, may want the following as a solution
## fastq2faqual.py --fastq /path/to/polished_assembly.fastq --fa --qual --out ${outpre}

SETUPARGS=""
MAPARGS=""
SUPPORTARGS=""
EXTRACTARGS=""
ASSEMBLYARGS="--nproc=${THREADS}"
OUTPUTARGS=""

if ${GAPSONLY}; then SUPPORTARGS+="--capturedOnly "; fi
if ${SPANONLY}; then SUPPORTARGS+="--spanOnly "; fi

run setup ${PROTOCOL} ${SETUPARGS}
run mapping ${PROTOCOL} ${MAPARGS} 
run support ${PROTOCOL} ${SUPPORTARGS}
run extraction ${PROTOCOL} ${EXTRACTARGS}
run assembly ${PROTOCOL} ${ASSEMBLYARGS}
run output ${PROTOCOL} ${OUTPUTARGS}

