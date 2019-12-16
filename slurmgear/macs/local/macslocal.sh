#!/bin/bash

function help {
    echo "
        Usage: ${0} -m:f:e:r:d:s:a:q:x:I:M:T:C:1:2:3:4:5:6:7:8:9:0hS
        -t with argument = treatment files
        -c with argument = control files
        -e with argument = --extsize for macs2 pileup (Default: 350)
        -S with argument = scaling factor used with macs2 bdgcmp
        -p with argument = pseudocount to use with macs2 bdgcmp. The pseudocount used for calculating logLR, logFE or FE. The count will be applied after normalization of sequencing depth. DEFAULT: 0.0, no pseudocount is applied.
        -m with argument = method for macs2 bdgcmp. Choose from: ppois,qpois,subtract,logFE,FE,logLR. Default = qpois.
        -B Use this to indicate --both-direction true for macs2 pileup (Default: False)
        -d with argument = NOT NEEDED NOW --- path to MACS2 (above bin)
        -s with argument = path to SCRIPTS DIR
        -q with argument = Primary QOS for sbatch. (Default: epscor-condo)
        -x with argument = Secondary QOS for sbatch. (Default: biomed-sb-condo)
        -I with argument = Higher numbers skew this toward using primary QOS more than secondary. Setting to 2 would be even split. (Default: 9)
        -M with argument = how much memory to tell sbatch. (Default: 8g)
        -T with argument = how much time to tell sbatch. (Default: 12:00:00)
        -C with argument = how many cpus/threads to tell sbatch. (Default: 2)
        -h help - returns this message; also returns this when no arguments given
        -S SUBMIT JOBS TO SLURM.

       In general, provide abs paths or paths from HOME rather than relative path from pwd unless it is in pwd or subdir.
"
}

##############################################################################
## TRIGGER HELP IF NO ARGS
##############################################################################
if [ $# -eq 0 ]; then help; exit; fi



##############################################################################
## DEFAULTS
##############################################################################
## DEFAULTS SPECIFIC TO MALIGNER SLURM FILE
JTHREADS=2
JMEM=8g
JTIME=12:00:00
SUBMITTOSLURM=false


## DEFAULTS - TYPICALLY USED IN ALL SLURM_X.sh FILES
CLEAN=false
HELP=false
IMAX=9
QOS1=epscor-condo
QOS2=biomed-sb-condo
SCRIPTS=`abspath.py ${0} --split | awk '{print $1}'`
SLURMOUTDIR=slurmout
EXIT=false

##############################################################################
## GET OPTS
##############################################################################
while getopts "m:f:e:r:d:s:a:q:x:I:M:T:C:1:2:3:4:5:6:7:8:9:0hS" arg; do
    case $arg in
        m) MAPSFOFN=$OPTARG;;
        h) HELP=true;;
        S) SUBMITTOSLURM=true;;
        *) help; exit;;
    esac
done


##############################################################################
## TRIGGER HELP IF OPTED FOR
##############################################################################
if ${HELP}; then help; exit; fi


##### PIPELINE FUNCTIONS
##############################################################################
## SUBMIT TO SLURM???
##############################################################################
if $SUBMITTOSLURM; then 
  source ${SCRIPTS}/slurm-macslocal-functions.sh; 
else 
  source ${SCRIPTS}/bash-macslocal-functions.sh; 
fi

##############################################################################
## RUNN PIPELINE
##############################################################################

## LOOP
i=0
while read ASM; do
  i=$(( $i+1 ))
  if [ $i -eq $IMAX ]; then QOS=${QOS2}; i=0; else QOS=${QOS1}; fi
  if [[ "$ASM" == *.fasta ]]; then BASE=`basename $ASM .fasta`; fi
  if [[ "$ASM" == *.fa ]]; then BASE=`basename $ASM .fa`; fi
  echo $BASE; 
  if [ ! -d $BASE ]; then mkdir $BASE; fi
  cd $BASE;
    MAIN=$PWD
    if $SUBMITTOSLURM && [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
    OUT=${MAIN}/${SLURMOUTDIR}
    ASM=`readlink -f ${ASM}`
    ### PIPELINE
    CLEAN1DEP=afterok
    echo convert_asm; convert_asm
    echo map_align; map_align
    echo merge_maps; merge_maps
    echo score_map_alns; score_map_alns
    echo map_alns_to_bdg; map_alns_to_bdg
    echo clean_up_map_aln; clean_up_map_aln
  cd ../
done < $ASMFOFN




