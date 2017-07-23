#!/bin/bash


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
        echo ${PWD}/fastaloc_${i}.${BASE}.${REC_ENZ}.smoothed.maps
      done < $FASTAFOFN >> ${MAPSFOFN}
      ## SUBMIT JOB
      QCONVDEP=`sbatch -J convertqueries -o convertqueriestomaps.slurm.%A.out --mem=8g --time=4:00:00 -c 2 --qos=$QOS \
       --export=REC_ENZ=${REC_ENZ},REC_SEQ=${REC_SEQ},MIN_FRAG_SIZE=${MIN_FRAG_SIZE},FASTAFOFN=${FASTAFOFN},MAPSFOFN=${MAPSFOFN} \
       ${SCRIPTS}/fa2map-while.sh | awk '{print $4}'`
      cd ../
    fi
    export QCONVDEP
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
      if $HASFASTAFOFN; then DEPENDS=--dependency=afterok:${QCONVDEP} ; else DEPENDS=""; fi
      CONVDEP=`sbatch ${DEPENDS} -J ${BASE}_convertasm -o ${OUT}/convertasm.slurm.%A.out --mem=8g --time=2:00:00 -c 2 --qos=$QOS \
         --export=ASM_FASTA=${ASM},BASE=${BASE},REC_ENZ=${REC_ENZ},REC_SEQ=${REC_SEQ},MIN_FRAG_SIZE=${MIN_FRAG_SIZE} ${SCRIPTS}/fa2map-for-slurm-mapalign.sh | awk '{print $4}'`
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
          --export=ASM_MAP=${ASM_MAP},SMOOTH_MAPS=${mapfilename},REC_SEQ=${REC_SEQ},REC_ENZ=${REC_ENZ},MIN_FRAG_SIZE=${MIN_FRAG_SIZE},QUERY_MISS_PENALTY=${QUERY_MISS_PENALTY},REF_MISS_PENALTY=${REF_MISS_PENALTY},QUERY_MAX_MISSES=${QUERY_MAX_MISSES},REF_MAX_MISSES=${REF_MAX_MISSES},SD_RATE=${SD_RATE},MIN_SD=${MIN_SD},MAX_SCORE_PER_INNER_CHUNK=${MAX_SCORE_PER_INNER_CHUNK},MAX_ALIGNMENTS_PER_QUERY=${MAX_ALIGNMENTS_PER_QUERY},PATH=${PATH},PYTHONPATH=${PYTHONPATH} \
          ${SCRIPTS}/maligner_dp.for_slurm_mapalign.sh | awk '{print $4}'`
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


