## Highest prob of needing changing
ASMFASTA=""
KINETICSCRIPTS=~/software/kineticsTools/scripts


## Optionally change these
QOS=epscor-condo
NAMES=contignames.txt
MAX=9
PVALUE=0.01
IDENTIFY=7
THREADS=4
TIME=24:00:00
MEM=30g
TOPTHREADS=1

## ABSPATHING
KINETICSCRIPTS=`readlink -f ${KINETICSCRIPTS}`
ASMFASTA=`readlink -f ${ASMFASTA}`
TOPKC=${KINETICSCRIPTS}/kinetics-contig.sh

## VARS THAT MAKE ASSUMPSOINS
TOPDIR=`readlink -f ../../..`
REF=`readlink -f ${ASMFASTA}`
REFBASE=`basename ${REF} .fasta`
INPUT=${TOPDIR}/${REFBASE}/out_all.cmp.h5 ## INPUT SHOULD ALREADY BE IN THIS REF_BASENAME_DIR


for var in QOS TOPDIR REF REFBASE INPUT NAMES MAX PVALUE IDENTIFY THREADS TIME MEM TOPKC TOPTHREADS; do
  echo ${var} ${!var}
done > vars.cfg
