## Highest prob of needing changing
ASMFASTA=""
KINETICSCRIPTS=~/software/kineticsTools/scripts


## Optionally change these
#QOS=epscor-condo
QOS=biomed-condo
NAMES=contignames.txt
MAX=9
PVALUE=0.01
IDENTIFY=7
THREADS=1
TIME=168:00:00
MEM=30g

## ABSPATHING
KINETICSCRIPTS=`readlink -f ${KINETICSCRIPTS}`
ASMFASTA=`readlink -f ${ASMFASTA}`
KC=${KINETICSCRIPTS}/kinetics-contig.sh

## VARS THAT MAKE ASSUMPSOINS
TOPDIR=`readlink -f ../../..`
REF=`readlink -f ${ASMFASTA}`
REFBASE=`basename ${REF} .fasta`
INPUT=${TOPDIR}/${REFBASE}/out_all.cmp.h5 ## INPUT SHOULD ALREADY BE IN THIS REF_BASENAME_DIR


for var in QOS TOPDIR REF REFBASE INPUT NAMES MAX PVALUE IDENTIFY THREADS TIME MEM KC; do
  echo ${var} ${!var}
done > vars.cfg
