#!/bin/bash

if [ $# -eq 0 ]; then echo "Provide previx to read fastx file"; exit; fi

PREID=${1} ## ONLY ARGUMENT TO SCRIPT
RELREADLOC=../${PREID}.fastq.gz

W=100
S=100
Q=0
M=32g
MUSCALE=0.4

PUFFERFISH=~/software/sciaratools/sciara-project-tools/pufferfish/pufferfish/pufferfish_main.py 
GCMED=~/software/sciaratools/sciara-project-tools/pufferfish/pufferfish/gc_median.py 
CNVCOR=~/software/sciaratools/sciara-project-tools/pufferfish/pufferfish/hmm_state_mean_correction.py 
PICARD=~/software/picardtools/picard-tools-2.1.1/picard.jar

###DATA=/gpfs/scratch/msokka/signalcorrection/data/genomes/hg19/
GENOME=/gpfs/scratch/jurban/nsseq/hg19_rdna_included.genome
BT2=/gpfs/scratch/jurban/nsseq/bt2index_rdna-included/hg19_index
GSEQ=/gpfs/scratch/jurban/nsseq/hg19_rdna-included.fa




### FILES MADE/USED LATER
READS=${PREID}.fastq.gz
FULLID=${PREID}.q${Q}.w${W}.s${S}
RD=${FULLID}.bedGraph
GCTABLE=${FULLID}.GCTABLE.txt

GCSUBTRACT_PRECN=${FULLID}.gc_subtract_precnnorm.bedGraph
GCSUBTRACT_PRECN=${FULLID}.gc_subtract_precnnorm.bedGraph
GCMEDS_PRECN=${FULLID}.gc_medians_control_precnnorm.bedGraph
GCMADS_PRECN=${FULLID}.gc_mads_control_precnnorm.bedGraph
Z_PRECN=${FULLID}.gc_zscores_precnnorm.bedGraph
FE_PRECN=${FULLID}.gc_FE_over_gcmedians_precnnorm.bedGraph
SIG_DIST_PRECN=${FULLID}.GC_DISTRIBUTION_precnnorm.txt
GCSTATS_PRECN=${FULLID}.GC_STATS_precnnorm.txt

STATES1=${FULLID}.cnvstates.bedGraph
STATEMEANS=${FULLID}.cnvstatemeans.bedGraph
NORMBDG=${FULLID}.cnv-normalized-gcmedfe.bedGraph
NORMRD=${FULLID}.cnv-normalized-read-depth.bedGraph
ZSCORES_POSTCN_ORIG=${FULLID}.gc_zscores_postcnnorm_origGC.bedGraph
ZSCORES_POSTCN_RECALC=${FULLID}.gc_zscores_postcnnorm_recalcGC.bedGraph

GCTABLE_POSTCN=${FULLID}.GCTABLE_postcnnorm.txt

GCSUBTRACT_POSTCN=${FULLID}.gc_subtract_postcnnorm.bedGraph
GCMEDS_POSTCN=${FULLID}.gc_medians_control_postcnnorm.bedGraph
ZSCORES_POSTCN=${FULLID}.gc_zscores_postcnnorm.bedGraph
FE_POSTCN=${FULLID}.gc_FE_over_gcmedians_postcnnorm.bedGraph
SIG_DIST_POSTCN=${FULLID}.GC_DISTRIBUTION_postcnnorm.txt
GCSTATS_POSTCN=${FULLID}.GC_STATS_postcnnorm.txt

STATES2=${FULLID}.cnvstates.altprotocol.bedGraph
STATEMEANS2=${FULLID}.cnvstatemeans.altprotocol.bedGraph
NORMBDG2=${FULLID}.cnv-normalized-RD1.bedGraph
NORMRD2=${FULLID}.cnv-normalized-RD2.bedGraph  ######## MIGHT NOT NEED THIS ONE
GCTABLE_POSTCN2=${FULLID}.GCTABLE_postcnnorm.altprotocol.txt
GCSUBTRACT_POSTCN2=${FULLID}.gc_subtract_postcnnorm.altprotocol.bedGraph
GCMEDS_POSTCN2=${FULLID}.gc_medians_control_postcnnorm.altprotocol.bedGraph
ZSCORES_POSTCN_RECALC2=${FULLID}.gc_zscores_postcnnorm.altprotocol.bedGraph
FE_POSTCN2=${FULLID}.gc_FE_over_gcmedians_postcnnorm.altprotocol.bedGraph
SIG_DIST_POSTCN2=${FULLID}.GC_DISTRIBUTION_postcnnorm.altprotocol.txt
GCSTATS_POSTCN2=${FULLID}.GC_STATS_postcnnorm.altprotocol.txt



echo Mkdir and go there.
mkdir ${PREID}
cd ${PREID}
ln -s ${RELREADLOC} .


echo MAP READS
${PUFFERFISH} mapreads --dry --threads 8 -b ${BT2} ${READS}
##${PUFFERFISH} mapreads --threads 8 -b ${BT2} ${READS}


echo GET COV
## Note: with the number of reads, I can easily get an avg of >= 1 read in 10 bp windows -- nonetheless I will try 100 first.....
${PUFFERFISH} getcov --dry -f $PICARD -g $GENOME -w $W -s $S -Q 0 -m $M --rmdup --clean ${PREID}.bam 
${PUFFERFISH} getcov -f $PICARD -g $GENOME -w $W -s $S -Q 0 -m $M --rmdup --clean ${PREID}.bam 


echo GET GC CONTENT OF BINS - PUT IN TABLE WITH NSSEQ VALUES
bedtools makewindows -w $W -s $S -g $GENOME | bedtools nuc -fi $GSEQ -bed - | sortBed -i - | cut -f 5 | paste <( sortBed -i $RD ) - > ${GCTABLE}


echo GET GC MEDIANS
${GCMED} -f ${GCTABLE} --gccol 5 --sigcol 4 --bdg ${GCSUBTRACT_PRECN} --control ${GCMEDS_PRECN} --mad ${GCMADS_PRECN} --zscore ${Z_PRECN} --medfe ${FE_PRECN} \
	--dist ${SIG_DIST_PRECN} > ${GCSTATS_PRECN}



echo FIND COPY NUMBER STATES
${PUFFERFISH} puffcn -7 --late ${FE_PRECN} --mu 0.125,0.1428571,0.1666667,0.2,0.25,0.3333333,0.5,1,2,3,4,5,6,7,8 --mu_scale ${MUSCALE} \
	--special_idx 7 --init_special 0.9 --leave_special_state 2000000 --leave_other 0.00001,0.000000000001 > ${STATES1}

echo GET BDG OF STATEMEANS
${CNVCOR} --signal ${FE_PRECN} --states ${STATES1} --levels ${STATEMEANS} --normbdg ${NORMBDG}

echo NORMALIZE READ DEPTH TO COPY NUMBER
paste $RD ${STATEMEANS} | awk 'OFS="\t" {print $1,$2,$3,$4/$8}' > ${NORMRD}

echo NORMALIZE CNV-NORM-RD TO GC WITH PRE-EXISTING GC MEDS AND MADS
paste ${NORMRD} ${GCMEDS_PRECN} ${GCMADS_PRECN} | awk '$12!=0' | awk 'OFS="\t" {print $1,$2,$3,($4-$8)/$12}' > ${ZSCORES_POSTCN_ORIG}

echo ALSO TRY RE-CALCULATING GC STATS FOR ZSCORES

echo 	GET GC CONTENT OF BINS - PUT IN TABLE WITH NSSEQ VALUES
cut -f 5 ${GCTABLE} | paste <( sortBed -i ${NORMRD} ) - > ${GCTABLE_POSTCN}

echo	GET GC STATS
${GCMED} -f ${GCTABLE_POSTCN} --gccol 5 --sigcol 4 --bdg ${GCSUBTRACT_POSTCN} --control ${GCMEDS_POSTCN} --zscore ${ZSCORES_POSTCN_RECALC} --medfe ${FE_POSTCN} \
	--dist ${SIG_DIST_POSTCN} > ${GCSTATS_POSTCN}



echo ALT PROCEDURE OF ESTIMATING COPY NUMBER W/O GCMEDFE -- MEDIAN NORM INSTEAD
echo FIND COPY NUMBER STATES
${PUFFERFISH} puffcn -4 -bw 500 --late ${RD} --mu 0.125,0.1428571,0.1666667,0.2,0.25,0.3333333,0.5,1,2,3,4,5,6,7,8 --mu_scale ${MUSCALE} \
	--special_idx 7 --init_special 0.9 --leave_special_state 2000000 --leave_other 0.00001,0.000000000001 > ${STATES2}


echo GET BDG OF STATEMEANS
${CNVCOR} --signal ${RD} --states ${STATES2} --levels ${STATEMEANS2} --normbdg ${NORMBDG2}

echo NORMALIZE READ DEPTH TO COPY NUMBER
paste $RD ${STATEMEANS2} | awk 'OFS="\t" {print $1,$2,$3,$4/$8}' > ${NORMRD2}

echo GET GCSTATS NOW -- AFTER CN NORM
cut -f 5 ${GCTABLE} | paste <( sortBed -i ${NORMRD2} ) - > ${GCTABLE_POSTCN2}

echo GET DONE
${GCMED} -f ${GCTABLE_POSTCN2} --gccol 5 --sigcol 4 --bdg ${GCSUBTRACT_POSTCN2} --control ${GCMEDS_POSTCN2} --zscore ${ZSCORES_POSTCN_RECALC2} --medfe ${FE_POSTCN2} \
	--dist ${SIG_DIST_POSTCN2} > ${GCSTATS_POSTCN2}




