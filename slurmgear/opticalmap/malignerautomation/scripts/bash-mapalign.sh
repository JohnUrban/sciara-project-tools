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
## FUNCTION:   FASTA QUERY FOFN TO SMOOTHED MAPS
##############################################################################
function convert_queries {
    D=query_maps
    if $HASFASTAFOFN; then
      if [ -d $D ]; then rm -r $D; fi
      mkdir $D
      cd $D
      ## Add names of future smooth mapped files to maps.fofn now so jobs are launched for them later
      i=0
      while read fastaloc; do
        let i++
        if [[ "$fastaloc" == *.fasta ]]; then BASE=`basename ${fastaloc} .fasta`; 
        elif [[ "$fastaloc" == *.fa ]]; then BASE=`basename ${fastaloc} .fa`;
        else BASE=query; fi
        ##OUT_PFX=${PWD}/fastaloc_${i}.${BASE}.${REC_ENZ}
        OUT_PFX=fastaloc_${i}.${BASE}.${REC_ENZ}
        make_insilico_map -o $OUT_PFX $fastaloc $REC_SEQ 2>> log
        smooth_maps_file -m $MIN_FRAG_SIZE ${OUT_PFX}.maps > ${OUT_PFX}.smoothed.maps 2>>log ;
      done < $FASTAFOFN >> ${MAPSFOFN}
      cd ../
    fi
}

##############################################################################
## FUNCTION:   FASTA ASM TO SMOOTHED MAPS
##############################################################################
function convert_asm {
    D=asm_map
    if $CONVERTASM; then
      if [ -d $D ]; then rm -r $D; fi
      mkdir $D
      cd $D
      # outputs
      ASM_OUT_PFX=${BASE}.${REC_ENZ}
      # convert the asm fasta file to the Maligner maps format and smooth the maps file by merging consecutive fragments that are less than 1kb
      make_insilico_map -o $ASM_OUT_PFX $ASM $REC_SEQ 2>> log
      smooth_maps_file -m $MIN_FRAG_SIZE ${ASM_OUT_PFX}.maps > ${ASM_OUT_PFX}.smoothed.maps 2>>log
      cd ../
    fi
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
     while read SMOOTH_MAPS; do
       mapfilebase=`basename $mapfilename`
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
     done < $MAPSFOFN 
     cd ../
    fi
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
       F1=`ls ${ALN}/*aln | head -n 1`
       head -n 1 $F1 > all.bionano.smoothed.maps.aln

       for f in ${ALN}/*.smoothed.maps*aln; do
           tail -n +2 $f >> all.bionano.smoothed.maps.aln
       done
       tail -n +2 all.bionano.smoothed.maps.aln | wc -l > num_alignments.txt
     cd ../
    fi
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
       tail -n +2 ${ALL} | cut -f 19 | awkSum > score.txt
     cd ../
    fi
}

##############################################################################
## BEDGRAPH
##############################################################################
function map_alns_to_bdg {
    if $MERGE; then
    D=merge
    if [ ! -d $D ]; then mkdir $D; fi
    cd $D
      G=${BASE}.genome
      faSize -detailed $ASM > $G    ###{BASE}.genome

      tail -n +2 $ALL | awk 'OFS="\t" { if ($10<0) print $2,0,$11; else print $2,$10,$11}' | sortBed -i - | genomeCoverageBed -i - -g ${BASE}.genome -bg > ${ALL}.bedGraph
      awk '{s+=($3-$2)}END{print s}' ${ALL}.bedGraph > span.txt
      awk '{s+=($3-$2)*$4}END{print s}' ${ALL}.bedGraph > total_base_cov.txt
      ## calculate metrics
      score=`head -n 1 score.txt`
      span=`cat span.txt`
      cov=`cat total_base_cov.txt`
      num=`cat num_alignments.txt`
      scorecov=`python -c "print 1e4*$score/$cov.0"`
      scorenum=`python -c "print $score/$num.0"`
      asmsize=`awk '{s+=$2}END{print s}' $G`
      spanasm=`python -c "print $span/$asmsize.0"`
      covasm=`python -c "print $cov/$asmsize.0"`
      covspan=`python -c "print $cov/$span.0"`
      covnum=`python -c "print $cov/$num.0"`
      ## populate allstats.txt
      echo -e score"\t"$score > allstats.txt
      echo -e span"\t"$span >> allstats.txt
      echo -e total_cov"\t"$cov >> allstats.txt
      echo -e num_aln"\t"$num >> allstats.txt
      echo -e 1e4xScore/Cov"\t"$scorecov >> allstats.txt
      echo -e Score/Num"\t"$scorenum >> allstats.txt
      echo -e Span/Asm"\t"$spanasm >> allstats.txt
      echo -e Cov/Asm"\t"$covasm >> allstats.txt
      echo -e Cov/Span"\t"$covspan >> allstats.txt
      echo -e Cov/Num"\t"$covnum >> allstats.txt
    cd ../
    fi
}

##############################################################################
## CLEAN UP
##############################################################################
function clean_up_map_aln {
    if $CLEAN; then
      cd merge/
        RM=false
        F=allstats.txt
        if [ -f $F ]; then X=`awk '$2!=""' $F | wc -l`; if [ $X -ge 9 ]; then RM=true; fi; fi
        if $RM; then 
          rm *aln *bedGraph *genome; 
          rm -r ../aln/ ../asm_map/
        fi
      cd ../

    fi
}





##############################################################################
## RUNN PIPELINE
##############################################################################
## FIRST CONVERT QUERIES IF NEED BE
QOS=${QOS1}
convert_queries
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
    if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
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




