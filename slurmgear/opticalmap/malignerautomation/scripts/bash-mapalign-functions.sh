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
        echo ${PWD}/fastaloc_${i}.${BASE}.${REC_ENZ} >> ${MAPSFOFN}
        OUT_PFX=fastaloc_${i}.${BASE}.${REC_ENZ}
        make_insilico_map -o $OUT_PFX $fastaloc $REC_SEQ 2>> log
        smooth_maps_file -m $MIN_FRAG_SIZE ${OUT_PFX}.maps > ${OUT_PFX}.smoothed.maps 2>>log ;
      done < $FASTAFOFN 
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
     while read SMOOTH_MAPS; do
       mapfilebase=`basename $SMOOTH_MAPS`
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



