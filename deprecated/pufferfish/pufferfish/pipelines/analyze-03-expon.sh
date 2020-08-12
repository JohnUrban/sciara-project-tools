#!/bin/bash

## NEED TO REMOVE CHRY, CHRM, AND GAPS FROM NS-SIGNAL AND GC STUFF ETC -- BEFORE CALCULATING MEDIANS, HMM, ETC

if [ $# -eq 0 ]; then echo "Provide previx to read fastx file"; exit; fi

PREID=${1} ## ONLY ARGUMENT TO SCRIPT
RELREADLOC=../${PREID}.fastq.gz

W=100
S=100
Q=0
M=32g
MAPPABILITY_BDG=/gpfs/scratch/jurban/nsseq/w100.s100.mappability.bedGraph

## AUTODETECT MUSCALE
##MUSCALE=0.4
#RANDOM_INTERVAL_LEN=3000
#N_RANDOM_INTERVALS=1000000
RANDOM_INTERVAL_LEN=30000
N_RANDOM_INTERVALS=100000
RANDBED=random_N${N_RANDOM_INTERVALS}_L${RANDOM_INTERVAL_LEN}.bed
RANDSTATSBED=random_N${N_RANDOM_INTERVALS}_L${RANDOM_INTERVAL_LEN}.stats.bed

PUFFERFISH=~/software/sciaratools/sciara-project-tools/pufferfish/pufferfish/pufferfish_main.py 
GCMED=~/software/sciaratools/sciara-project-tools/pufferfish/pufferfish/gc_median.py 
CNVCOR=~/software/sciaratools/sciara-project-tools/pufferfish/pufferfish/hmm_state_mean_correction.py 
PICARD=~/software/picardtools/picard-tools-2.1.1/picard.jar

DATA=/gpfs/scratch/msokka/signalcorrection/data/genomes/hg19
GENOME=${DATA}/hg19_rdna_excluded.genome
BT2=${DATA}/bt2index_rdna-excluded/hg19_index
GSEQ=${DATA}/hg19_rdna-excluded.fa
EXCLUDED_REGIONS=/gpfs/scratch/jurban/nsseq/excluded_regions.bed

### FILES MADE/USED LATER
READS=${PREID}.fastq.gz
BAM=${PREID}.bam
BED=${PREID}.bed
PILEUP=${PREID}.bedGraph
FULLID=${PREID}.q${Q}.w${W}.s${S}
RD=${FULLID}.bedGraph
GCTABLE=${FULLID}.GCTABLE.txt

GCSUBTRACT_PRECN=${FULLID}.gc_subtract_precnnorm.bedGraph
GCSUBTRACT_PRECN=${FULLID}.gc_subtract_precnnorm.bedGraph
GCMEDS_PRECN=${FULLID}.gc_medians_control_precnnorm.bedGraph
GCMADS_PRECN=${FULLID}.gc_mads_control_precnnorm.bedGraph
Z_PRECN=${FULLID}.gc_zscores_precnnorm.bedGraph

FE_PRECN=${FULLID}.gc_FE_over_gcmeans_precnnorm.bedGraph

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

WINDOWS=w${W}.s${S}.bedGraph
GCWINDOWS=gc.w${W}.s${S}.bedGraph
MAPWINDOWS=map.w${W}.s${S}.bedGraph


echo Mkdir and go there.
#mkdir ${PREID}
cd ${PREID}
#ln -s ${RELREADLOC} .


echo MAP READS
${PUFFERFISH} mapreads --dry --threads 8 -b ${BT2} ${READS}
##${PUFFERFISH} mapreads --threads 8 -b ${BT2} ${READS}

echo MACS2 RMDUP
## NOTE: this operation ruins the rDNA analysis -- that needs to be done separately -- possibly by extracting before this step and treating in parallel (copy number can be corrected separately -- possibly by rDNA * genome_mean/rDNA_mean or rDNA * genome_median/rDNA_median), but then GCnorm w/ GCmeds learned from whole genome)....
#macs2 filterdup -i ${BAM} -g hs --keep-dup 1 -o ${BED}

echo MACS2 PILEUP
#macs2 pileup -i ${BED} -o ${PILEUP} --extsize 350

echo PREPARING EXCLUDED REGIONS
sortBed -i ${EXCLUDED_REGIONS} | mergeBed -i - > excluded_regions.bed
EXCLUDED_REGIONS=excluded_regions.bed

echo SORTING PILEUP AND SUBTRACTING EXCLUDED REGIONS
#sortBed -i ${PILEUP} | bedtools subtract -sorted -a - -b ${EXCLUDED_REGIONS} > ${PILEUP}.1
##################bedtools subtract -sorted -a ${PILEUP} -b <(sortBed -i ${EXCLUDED_REGIONS}) > ${PILEUP}.1
#mv ${PILEUP}.1 ${PILEUP}

echo MAKING w$W s$S WINDOWS....
#bedtools makewindows -w $W -s $S -g $GENOME | sortBed -i - | intersectBed -sorted -v -a - -b <(sortBed -i ${EXCLUDED_REGIONS}) > $WINDOWS

echo OBTAINING GC INFO FOR WINDOWS
#bedtools nuc -fi $GSEQ -bed $WINDOWS | sortBed -i - | awk 'OFS="\t" {print $1,$2,$3,$5}' > $GCWINDOWS

echo MAPPING PILEUP VALUES ONTO w$W s$S WINDOWS FOR INITIAL READ DEPTH VALUES 
#bedtools map -c 4 -o sum -a ${WINDOWS} -b ${PILEUP} | sortBed -i - > ${RD}

## May want to consider correcting for mappability -- particularly if you do Qfiltering....
## http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig
## http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
## http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
## http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377
## If you have filtered for q>=2 (or something more stringent) -- do bin/median(mappability) -- if you did not filter for mapq -- try nothing or bin*median(mappability)
## I can see that NS-seq bins are negatively correlated with mappability bins -- i.e. when the region becomes less unique, more reads pileup... not great
## Probably the best fix is filtering q>=2 (or more stringent) b/c the effect of copy number gains can really send multi-sites out of whack....

echo OBTAINING MEDIAN MAPPABILITY FOR START SITES IN EACH w$W s$S WINDOW
#bedtools map -c 4 -o median -a ${WINDOWS} -b ${MAPPABILITY_BDG} | sortBed -i - > ${MAPWINDOWS}

echo SLIGHT CORRECTION OF READ DEPTH FOR MAPPABILITY AND FORMATION OF GCTABLE
## Should produce: chr, stat, end GC, SIG
#paste ${MAPWINDOWS} ${RD} ${GCWINDOWS} | awk '$4 != "." && $8!="." && $12 != "." {OFS="\t"; print $1,$2,$3,$12,$4*$8}' > ${GCTABLE}


echo ENSURING EXCLUDED REGIONS ARE REMOVED FROM ANALYSIS
#bedtools intersect -sorted -v -a ${GCTABLE} -b <(sortBed -i ${EXCLUDED_REGIONS}) > ${GCTABLE}.1
#mv ${GCTABLE}.1 ${GCTABLE}

echo CONSTRUCTING FINAL READ DEPTH FILE FROM GCTABLE FOR LATER USE
#cut -f 1,2,3,5 ${GCTABLE} > ${RD}


echo GET GC MEDIANS
#${GCMED} -f ${GCTABLE} --gccol 4 --sigcol 5 --meanfe ${FE_PRECN} --dist ${SIG_DIST_PRECN} > ${GCSTATS_PRECN}

echo LEARN MEAN RATIO OF STDEV:MEAN IN $N_RANDOM_INTERVALS WINDOWS OF LENGTH $RANDOM_INTERVAL_LEN BP FOR MUSCALE PARAMETER
#bedtools random -l ${RANDOM_INTERVAL_LEN} -n ${N_RANDOM_INTERVALS} -g $GENOME  | sortBed -i - | intersectBed -sorted -v -a - -b ${EXCLUDED_REGIONS} | awk 'OFS="\t" {print $1,$2,$3}' > ${RANDBED}
#bedtools map -c 4 -o mean,sstdev -a ${RANDBED} -b ${FE_PRECN} | awk '$4!=0 && $5!=0 && $4!="." && $5!="."' > ${RANDSTATSBED}
#MUSCALE=`awk '{SUM+=$5/$4}END{print SUM/NR}' ${RANDSTATSBED}`
#echo "MUSCALE = $MUSCALE"

## If wanted to round up (ceiling)
##MUSCALE=`awk '{SUM+=$5/$4}END{print int((SUM/NR)+1)}' ${RANDSTATSBED}`
MUSCALE=1

echo FIND COPY NUMBER STATES
## expect to leave copy num 1 once in 200Mb+ (MCF7 had 24 large amplicons spanning from 0.2-12Mb: https://www.ncbi.nlm.nih.gov/pubmed/12414653 )
## 30% of MCF7 have CN!=1 - w/ mean amplicons ~582kb and mean deletion ~942kb -- w/8x more amps than dels -- giving weighted avg of 622 kb ::http://journals.plos.org/plosone/art$
## --> there were at least 99-100 breakpoints in that study indicating a change every 5-10 Mb on avg
${PUFFERFISH} puffcn -7 --late ${FE_PRECN} --mu 0.125,0.1428571,0.1666667,0.2,0.25,0.3333333,0.5,1,2,3,4,5,6,7,8,9,10,11,12 --mu_scale ${MUSCALE} --emodel exponential \
	--special_idx 7 --init_special 0.9 --leave_special_state 1000000000 --leave_other 0.000000000000001,0.000000000000000000001 > ${STATES1}

echo GET BDG OF STATEMEANS
${CNVCOR} --signal ${FE_PRECN} --states ${STATES1} --levels ${STATEMEANS} --normbdg ${NORMBDG}

echo NORMALIZE READ DEPTH TO COPY NUMBER
paste $RD ${STATEMEANS} | awk 'OFS="\t" {print $1,$2,$3,$4/$8}' > ${NORMRD}


echo 	GET GC CONTENT OF BINS - PUT IN TABLE WITH NSSEQ VALUES
## Should produce: chr, start, end, GC, SIG
cut -f 4 ${GCTABLE} | paste <( sortBed -i ${NORMRD} ) - | awk 'OFS="\t" {print $1,$2,$3,$5,$4}' > ${GCTABLE_POSTCN}

echo	GET GC STATS
${GCMED} -f ${GCTABLE_POSTCN} --gccol 4 --sigcol 5 --medfe ${FE_POSTCN} --dist ${SIG_DIST_POSTCN} > ${GCSTATS_POSTCN}







exit









echo NORMALIZE CNV-NORM-RD TO GC WITH PRE-EXISTING GC MEDS AND MADS
paste ${NORMRD} ${GCMEDS_PRECN} ${GCMADS_PRECN} | awk '$12!=0' | awk 'OFS="\t" {print $1,$2,$3,($4-$8)/$12}' > ${ZSCORES_POSTCN_ORIG}




echo ALT PROCEDURE OF ESTIMATING COPY NUMBER W/O GCMEDFE -- MEDIAN NORM INSTEAD
echo FIND COPY NUMBER STATES
${PUFFERFISH} puffcn -4 -bw 500 --late ${RD} --mu 0.125,0.1428571,0.1666667,0.2,0.25,0.3333333,0.5,1,2,3,4,5,6,7,8 --mu_scale ${MUSCALE} \
	--special_idx 7 --init_special 0.9 --leave_special_state 2000000 --leave_other 0.00001,0.000000000001 > ${STATES2}


echo GET BDG OF STATEMEANS
${CNVCOR} --signal ${RD} --states ${STATES2} --levels ${STATEMEANS2} --normbdg ${NORMBDG2}

echo NORMALIZE READ DEPTH TO COPY NUMBER
paste $RD ${STATEMEANS2} | awk 'OFS="\t" {print $1,$2,$3,$4/$8}' > ${NORMRD2}

echo GET GCSTATS NOW -- AFTER CN NORM
## Should produce: chr,	start, end, GC,	SIG
cut -f 4 ${GCTABLE} | paste <( sortBed -i ${NORMRD2} ) - | awk 'OFS="\t" {print $1,$2,$3,$5,$4}' > ${GCTABLE_POSTCN2}

echo GET DONE
${GCMED} -f ${GCTABLE_POSTCN2} --gccol 4 --sigcol 5 --bdg ${GCSUBTRACT_POSTCN2} --control ${GCMEDS_POSTCN2} --zscore ${ZSCORES_POSTCN_RECALC2} --medfe ${FE_POSTCN2} \
	--dist ${SIG_DIST_POSTCN2} > ${GCSTATS_POSTCN2}




