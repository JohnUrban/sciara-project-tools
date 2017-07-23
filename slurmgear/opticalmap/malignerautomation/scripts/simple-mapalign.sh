#!/bin/bash

##############################################################################
## FUNCTION: HELP
##############################################################################

function help {
    echo "
        Usage: ${0} -R:Q:e:r:d:1:2:3:4:5:6:7:8:9:h
        -R with argument = maps for reference
        -Q with argument = maps for query/queries
        -e with argument = REC_ENZYME (Default: BssSI)
        -r with argument = REC_SEQ (Default: CACGAG)
        -d with argument = path to MALIGNER dir (above bin/) (Default: ~/data/software/maligner/maligner)
        -s with argument = path to SCRIPTS DIR
        -1 with argument = MIN_FRAG_SIZE for maligner (Default: 1000)
        -2 with argument = QUERY_MISS_PENALTY for maligner (Default: 3.0)
        -3 with argument = REF_MISS_PENALTY for maligner (Default: 3.0)
        -4 with argument = QUERY_MAX_MISSES for maligner (Default: 5)
        -5 with argument = REF_MAX_MISSES for maligner (Default: 5)
        -6 with argument = SD_RATE for maligner (Default: 0.05)
        -7 with argument = MIN_SD for maligner (Default: 750)
        -8 with argument = MAX_SCORE_PER_INNER_CHUNK= for maligner (Default: 1.0)
        -9 with argument = MAX_ALIGNMENTS_PER_QUERY for maligner (Default: 1)
        -h help - returns this message; also returns this when no arguments given

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
REC_ENZ=BssSI
REC_SEQ=CACGAG
MIN_FRAG_SIZE=1000 # bp units
QUERY_MISS_PENALTY=3.0
REF_MISS_PENALTY=3.0
QUERY_MAX_MISSES=5
REF_MAX_MISSES=5
SD_RATE=0.05
MIN_SD=750
MAX_SCORE_PER_INNER_CHUNK=1.0
MAX_ALIGNMENTS_PER_QUERY=1

## DEFAULTS - TYPICALLY USED IN ALL SLURM_X.sh FILES
HELP=false
SCRIPTS=`abspath.py ${0} --split | awk '{print $1}'`

##############################################################################
## GET OPTS
##############################################################################
while getopts "R:Q:e:r:d:1:2:3:4:5:6:7:8:9:h" arg; do
    case $arg in
        R) ASM_MAP=$OPTARG;;
        Q) SMOOTH_MAPS=$OPTARG;;
        e) REC_ENZ=$OPTARG;;
        r) REC_SEQ=$OPTARG;;
        d) MALIGNER=$OPTARG;;
        s) SCRIPTS=$OPTARG;;
        1) MIN_FRAG_SIZE=$OPTARG;;
        2) QUERY_MISS_PENALTY=$OPTARG;;
        3) REF_MISS_PENALTY=$OPTARG;;
        4) QUERY_MAX_MISSES=$OPTARG;;
        5) REF_MAX_MISSES=$OPTARG;;
        6) SD_RATE=$OPTARG;;
        7) MIN_SD=$OPTARG;;
        8) MAX_SCORE_PER_INNER_CHUNK=$OPTARG;;
        9) MAX_ALIGNMENTS_PER_QUERY=$OPTARG;;
        h) HELP=true;;
        *) help; exit;;
    esac
done


##############################################################################
## TRIGGER HELP IF OPTED FOR
##############################################################################
if ${HELP}; then help; exit; fi


##############################################################################
## EXPORT MALIGNER PATHS
##############################################################################
## EXPORT MALIGNER INTO ENV
export PATH=${MALIGNER}/bin/:${MALIGNER}/build/bin/:$PATH
export PYTHONPATH=${MALIGNER}/lib/:$PYTHONPATH


##### PIPELINE FUNCTIONS
##source ${SCRIPTS}/bash-mapalign-functions.sh; 

##############################################################################
## RUNN PIPELINE
##############################################################################
## FIRST CONVERT QUERIES IF NEED BE

if [[ "$ASM" == *.fasta ]]; then BASE=`basename $ASM .fasta`; fi
if [[ "$ASM" == *.fa ]]; then BASE=`basename $ASM .fa`; fi
MAIN=$PWD
ASM=`readlink -f ${ASM}`
BASE2=`basename $SMOOTH_MAPS`
RMAP_OUT_PFX=${BASE2}
OUT_PFX=${BASE2}

       maligner_dp \
         -q $QUERY_MISS_PENALTY \
         -r $REF_MISS_PENALTY \
         --query-max-misses $QUERY_MAX_MISSES \
         --ref-max-misses $REF_MAX_MISSES \
         --max-score-per-inner-chunk $MAX_SCORE_PER_INNER_CHUNK \
         --sd-rate $SD_RATE \
         --min-sd $MIN_SD \
         --max-alignments $MAX_ALIGNMENTS_PER_QUERY \
         ${SMOOTH_MAPS} \
         ${ASM_MAP} \
         1> ${OUT_PFX}.aln 2> ${OUT_PFX}.log
         ###2>&1 1> ${OUT_PFX}.aln | tee ${OUT_PFX}.log





