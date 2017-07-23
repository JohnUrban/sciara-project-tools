#!/bin/bash

##############################################################################
## FUNCTION: HELP
##############################################################################

function help {
    echo "
        Usage: ${0} -m:f:e:r:d:s:a:q:x:I:M:T:C:1:2:3:4:5:6:7:8:9:0h
        -m with argument = maps.fofn file (has paths to all maps to be aligned to assembly) (Default: maps.fofn)
        -f with argument = fasta.fofn file -- alternative to -m/maps.fofn. Contains paths to fastas to be aligned to assembly after being converted to maps. If both maps.fofn and fasta.fofn files available, they will be combined after convrsion of fastas. (Default: fasta.fofn).
        -e with argument = REC_ENZYME (Default: BssSI)
        -r with argument = REC_SEQ (Default: CACGAG)
        -d with argument = path to MALIGNER dir (above bin/) (Default: ~/data/software/maligner/maligner)
        -s with argument = path to SCRIPTS DIR
        -a with argument = path to ASM FOFN (Default: input.fofn)
        -q with argument = Primary QOS for sbatch. (Default: epscor-condo)
        -x with argument = Secondary QOS for sbatch. (Default: biomed-sb-condo)
        -I with argument = Higher numbers skew this toward using primary QOS more than secondary. Setting to 2 would be even split. (Default: 9)
        -M with argument = how much memory to tell sbatch. (Default: 8g)
        -T with argument = how much time to tell sbatch. (Default: 12:00:00)
        -C with argument = how many cpus/threads to tell sbatch. (Default: 2)
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
        -0 Clean when done.

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
##CONFIG=$2
MAPSFOFN=maps.fofn
FASTAFOFN=fasta.fofn
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
JTHREADS=2
JMEM=8g
JTIME=12:00:00
MALIGNER=~/data/software/maligner/maligner

## DO NOT HAVE OPTIONS TO TURN OFF YET
CONVERTASM=true
MAPBIONANO=true
MERGE=true
## NOT YET USED -- score and bdg functions use MERGE for now....
SCORE=true
COVBDG=true

## DEFAULTS - TYPICALLY USED IN ALL SLURM_X.sh FILES
CLEAN=false
ASMFOFN=input.fofn
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
while getopts "m:f:e:r:d:s:a:q:x:I:M:T:C:1:2:3:4:5:6:7:8:9:0h" arg; do
    case $arg in
        m) MAPSFOFN=$OPTARG;;
        f) FASTAFOFN=$OPTARG;;
        e) REC_ENZ=$OPTARG;;
        r) REC_SEQ=$OPTARG;;
        d) MALIGNER=$OPTARG;;
        s) SCRIPTS=$OPTARG;;
        a) ASMFOFN=$OPTARG;;
        q) QOS1=$OPTARG;;
        x) QOS2=$OPTARG;;
        I) IMAX=$OPTARG;;
        M) JMEM=$OPTARG;;
        T) JTIME=$OPTARG;;
        C) JTHREADS=$OPTARG;;
        1) MIN_FRAG_SIZE=$OPTARG;;
        2) QUERY_MISS_PENALTY=$OPTARG;;
        3) REF_MISS_PENALTY=$OPTARG;;
        4) QUERY_MAX_MISSES=$OPTARG;;
        5) REF_MAX_MISSES=$OPTARG;;
        6) SD_RATE=$OPTARG;;
        7) MIN_SD=$OPTARG;;
        8) MAX_SCORE_PER_INNER_CHUNK=$OPTARG;;
        9) MAX_ALIGNMENTS_PER_QUERY=$OPTARG;;
        0) CLEAN=true;;
        h) HELP=true;;
        *) help; exit;;
    esac
done


##############################################################################
## TRIGGER HELP IF OPTED FOR
##############################################################################
if ${HELP}; then help; exit; fi


##############################################################################
## PROCESS ARGS WHERE NECESSARY
##############################################################################
HASMAPSFOFN=false
HASFASTAFOFN=false
if [ -f $MAPSFOFN ]; then HASMAPSFOFN=true; MAPSFOFN=`readlink -f ${MAPSFOFN}`; fi
if [ -f $FASTAFOFN ]; then HASFASTAFOFN=true; FASTAFOFN=`readlink -f ${FASTAFOFN}`; fi
if [ $HASMAPSFOFN  == false ] && [ $HASFASTAFOFN  == false ]; then echo "Could not find MAPSFOFN nor FASTAFOFN file(s). Exiting..."; exit; fi


##############################################################################
## EXPORT MALIGNER PATHS
##############################################################################
## EXPORT MALIGNER INTO ENV
export PATH=${MALIGNER}/bin/:${MALIGNER}/build/bin/:$PATH
export PYTHONPATH=${MALIGNER}/lib/:$PYTHONPATH



##### PIPELINE FUNCTIONS
##############################################################################
## FUNCTION:   FASTA ASM TO SMOOTHED MAPS
##############################################################################
function convert_asm {
    D=asm_map
    if $CONVERTASM; then
      if [ -d $D ]; then rm -r $D; fi
      mkdir $D
      cd $D
      CONVDEP=`sbatch -J ${BASE}_convertasm -o ${OUT}/convertasm.slurm.%A.out --mem=8g --time=2:00:00 -c 2 --qos=$QOS \
         --export=ASM_FASTA=${ASM},BASE=${BASE},REC_ENZ=${REC_ENZ},REC_SEQ=${REC_SEQ} ${SCRIPTS}/fa2map.sh | awk '{print $4}'`
      cd ../
    fi
    export CONVDEP
    export ASM_MAP=`readlink -f asm_map/`/${BASE}.${REC_ENZ}.smoothed.maps
}

##############################################################################
## MAP PRE-CONVERTED/PRE-SMOOTHED BIONANO MAPS
##############################################################################
function map_align {
    if $MAPBIONANO; then
     D=aln
     if [ ! -d $D ]; then mkdir $D; fi
     cd $D
     DEPENDS=""
     if $CONVERTASM; then
       DEPENDS=--dependency=afterok:${CONVDEP}
     fi
     while read mapfilename; do
       mapfilebase=`basename $mapfilename`
       CLEAN1DEP+=:`sbatch -J ${BASE}_${mapfilebase} $DEPENDS -o ${OUT}/${mapfilebase}.slurm.%A.out --mem=$JMEM --time=$JTIME -c $JTHREADS --qos=$QOS \
          --export=ASM_MAP=${ASM_MAP},SMOOTH_MAPS=${mapfilename},REC_SEQ=${REC_SEQ},PATH=${PATH} ${SCRIPTS}/maligner_dp.for_slurm_mapalign.sh | awk '{print $4}'`
     done < $MAPSFOFN 
     cd ../
    fi
    export CLEAN1DEP
    export ALN=`readlink -f aln/`
}

##############################################################################
## MERGE
##############################################################################
function merge_maps {
    if $MERGE; then
     D=merge
     if [ ! -d $D ]; then mkdir $D; fi
     cd $D
     DEPENDS=""
     if $MAPBIONANO; then
       DEPENDS=--dependency=${CLEAN1DEP}
     fi
     MERGEDEP=`sbatch -J ${BASE}_merge $DEPENDS -o ${OUT}/merge.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
        --export=ALN=${ALN} ${SCRIPTS}/merge.sh | awk '{print $4}'`
     CLEAN1DEP+=:${MERGEDEP}
     cd ../
    fi
    export CLEAN1DEP
    export ALL=`readlink -f merge/all.bionano.smoothed.maps.aln`
}

##############################################################################
## SCORE
##############################################################################
function score_map_alns {
    if $MERGE; then
     D=merge
     if [ ! -d $D ]; then mkdir $D; fi
     cd $D
     DEPENDS=""
     if $MERGE; then
       DEPENDS=--dependency=afterok:${MERGEDEP}
     fi
     SCOREDONE=`sbatch -J ${BASE}_score $DEPENDS -o ${OUT}/score.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
       --export=ALL=${ALL} ${SCRIPTS}/score.sh | awk '{print $4}'`
     CLEAN1DEP+=:$SCOREDONE
     cd ../
    fi
    export CLEAN1DEP
}

##############################################################################
## BEDGRAPH
##############################################################################
function map_alns_to_bdg {
    if $MERGE; then
    D=merge
    if [ ! -d $D ]; then mkdir $D; fi
    cd $D
    DEPENDS=""
    if $MERGE; then
      DEPENDS=--dependency=afterok:${SCOREDONE}
    fi
    CLEAN1DEP=$CLEAN1DEP:`sbatch -J ${BASE}_bdg $DEPENDS -o ${OUT}/bedgraph.slurm.%A.out --mem=2g --time=1:00:00 -c 1 --qos=$QOS \
      --export=ALL=${ALL},ASM=${ASM},BASE=${BASE} ${SCRIPTS}/convert2bedGraph.sh | awk '{print $4}'`
    cd ../
    fi
    export CLEAN1DEP
}

##############################################################################
## CLEAN UP
##############################################################################
function clean_up_map_aln {
    if $CLEAN; then
      CLEANDONE=`sbatch -J ${BASE}_bionano_clean --dependency=$CLEAN1DEP -o ${OUT}/clean.slurm.%A.out --mem=2g --time=1:00:00 -c 1 \
           --qos=${QOS} ${SCRIPTS}/clean.sh | awk '{print $4}'`
    fi
    export CLEANDONE
}





## LOOP
##############################################################################
## SUBMIT BATCH JOBS
##############################################################################
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
    if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
    OUT=${MAIN}/${SLURMOUTDIR}
    ASM=`readlink -f ${ASM}`
    ### PIPELINE
    CLEAN1DEP=afterok
    convert_asm
    map_align
    merge_maps
    score_map_alns
    map_alns_to_bdg
    clean_up_map_aln
  cd ../
done < $ASMFOFN




