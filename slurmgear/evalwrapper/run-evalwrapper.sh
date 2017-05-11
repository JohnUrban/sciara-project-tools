#!/bin/bash

## Input Assemblies FOFN
ASMFOFN=input.fofn
CLEANALL=true

## GETTING ABS PATH OF ASMFOFN
ASMFOFN=`readlink -f $ASMFOFN`

###################  SHORT READ  ###########################
# specify paths to lap read sample (LR1,LR2) and all reads (R1,R2)-- give dummy answers if will not be using (that will serve as place-holders)
LR1=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.1.fastq
LR2=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.2.fastq
R1=~/data/scratch/male-ilmn/data/ilmnraw/R1.fastq
R2=~/data/scratch/male-ilmn/data/ilmnraw/R2.fastq

## OPTIONS for what programs to use. 
## FILL IN   ==>"EvalThese"<==   BELOW WITH ONE OF THESE
ALL=eval.cfg
OnlyAle=eval.aleonly.cfg
OnlyBusco=eval.buscoOnly.cfg
OnlyLap=eval.laponly.cfg
OnlyReapr=eval.reapronly.cfg
OnlyReaprNoClean=eval.reapronly.noclean.cfg
OnlyReaprNoCleanAggressive=eval.reapronly.noclean.aggressive.cfg

## FILL IN WITH CORRECT VARIABLE
EvalThese=$ALL

## May need to adjust the following
SHORTSCRIPTS=/gpfs_home/jurban/software/sciaratools/sciara-project-tools/slurmgear/shortreadeval/scripts/
SHORTAUTO=${SHORTSCRIPTS}/auto-shortreadeval.sh
SHORTEVAL=${SHORTSCRIPTS}/eval.ARGS.sh
SHORTCONFIG=${SHORTSCRIPTS}/configs/${EvalThese}


############### BIONANO MALIGNER SECTION #################
BIONANOCLEAN=false
if $CLEANALL; then BIONANOCLEAN=true; fi

REC_ENZ=BssSI
REC_SEQ=CACGAG

BIONANOBASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/opticalmap/malignerautomation
BIONANOSCRIPTS=${BIONANOBASE}/scripts/
BIONANOCONFIGS=${BIONANOBASE}/configs/
BIONANOFOFNS=${BIONANOBASE}/fofns/

BIONANOCONFIG=${BIONANOCONFIGS}/maligner-config-sciara.cfg
MAPSFOFN=${BIONANOFOFNS}/bionanomaps.examp2.fofn
BIONANORUN=${BIONANOSCRIPTS}/auto-malign.sh




############### LONG READ SECTION #################
## LONG READ LOCATIONS
ONT=~/data/scratch/minion2016/fast5fastqs/allReadsFromAllONTlibsCombined.fastq
PACBIO=~/data/scratch/pac_bio_data/filt/all_subreads.fastq

## LONG2PE READ LOCATIONS
ONT1=/gpfs/data/sgerbi/jurban/scratch/minion2016/fast5fastqs/molreads/long2pe/ontmol-1.fastq
ONT2=/gpfs/data/sgerbi/jurban/scratch/minion2016/fast5fastqs/molreads/long2pe/ontmol-2.fastq
PACBIO1=/gpfs/data/sgerbi/jurban/scratch/pac_bio_data/filt/long2pe/pacbio-1.fastq
PACBIO2=/gpfs/data/sgerbi/jurban/scratch/pac_bio_data/filt/long2pe/pacbio-2.fastq

## RUN INFO LOCATIONS
LRBASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/longreadeval
LRSCRIPTS=${LRBASE}/scripts/
AUTOLR=${LRSCRIPTS}/auto-lrpipe.sh
LRCONFIGS=${LRBASE}/configs/
LRCONFIG=${LRCONFIGS}/longread-config-sciara.cfg

## OTHER OPTIONS
LRCLEAN=false
if $CLEANALL; then LRCLEAN=true; fi

############### TRANSCRIPT SECTION #################
TRANJOBPRE=transcript
TBLASTX=true
TRANSNJOBS=100
TRANSCLEAN=false
if $CLEANALL; then TRANSCLEAN=true; fi

TRANSBASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/transcripteval
TRANSSCRIPTS=${TRANSBASE}/scripts/
TRANSCONFIGS=${TRANSBASE}/configs/
TRANSFOFNS=${TRANSBASE}/fofns/

TRANSCONFIG=${TRANSCONFIGS}/trans-config-sciara.cfg ## does both blastn and tblastx
TRANSFOFN=${TRANSFOFNS}/
TRANSRUN=${TRANSSCRIPTS}/auto-trans.sh

TRANS1=~/data/illumina/generalTranscriptome/trinity/trinity_out_dir/Trinity.fasta
TRANS2=/gpfs/data/sgerbi/jurban/flies/dmel/dmel-all-transcript-r6.14.fasta
TRANS3=/gpfs/data/sgerbi/jurban/flies/anopheles_gambiae/anopheles-gambiae-pesttranscriptsagamp46.fa

TRANJOBPRE1=${TRANJOBPRE}_sciara_
TRANJOBPRE2=${TRANJOBPRE}_dmel_
TRANJOBPRE3=${TRANJOBPRE}_mosquito_

############### PEPTIDE SECTION #################
## Evaluate with peptides
PEPJOBPRE=peptide
PEPNJOBS=100
PEPCLEAN=false
if $CLEANALL; then PEPCLEAN=true; fi

PEPBASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/peptideval
PEPSCRIPTS=${PEPBASE}/scripts/
PEPCONFIGS=${PEPBASE}/configs/
PEPFOFNS=${PEPBASE}/fofns/

PEPCONFIG=${PEPCONFIGS}/peptide-config-sciara.cfg ## does both blastn and tblastx
PEPFOFN=${PEPFOFNS}/
PEPRUN=${PEPSCRIPTS}/auto-pep.sh

PEP2=/gpfs/data/sgerbi/jurban/flies/dmel/dmel-all-translation-r6.14.fasta 
PEP3=/gpfs/data/sgerbi/jurban/flies/anopheles_gambiae/anopheles-gambiae-pestpeptidesagamp46.fa 

PEPJOBPRE2=${PEPJOBPRE}_dmel_
PEPJOBPRE3=${PEPJOBPRE}_mosquito_

############### KNOWN SEQUENCES SECTION #################
## Also evaluate Known Seqs
## USE TRANS variables (e.g. TRANSSCRIPTS etc) for everything other than these 4 things
KNOWNJOBPRE=knownseqs_
KNOWNTBLASTX=false
KNOWNNJOBS=1
KNOWNCLEAN=false
if $CLEANALL; then KNOWNCLEAN=true; fi
KNOWNSEQS=/gpfs/data/sgerbi/jurban/sciaraknownseqs/allCoprophilaNTSeqOnNCBI.fa


############### RNASEQ SECTION #################
#$RNARUN $RNASCRIPTS $RNACONFIG $RNACLEAN $ASMFOFN $RNAFOFN

RNACLEAN=false
if $CLEANALL; then RNACLEAN=true; fi

RNABASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/rnaseqeval
RNASCRIPTS=${RNABASE}/scripts/
RNACONFIGS=${RNABASE}/configs/
RNAFOFNS=${RNABASE}/fofns/

RNACONFIG=${RNACONFIGS}/rnaseq-config-sciara.cfg
RNAFOFN=${RNAFOFNS}/reads.fofn
RNARUN=${RNASCRIPTS}/auto-rnaseqeval.sh

RNAFOFN=`readlink -f $RNAFOFN`



##############################################
##############################################
##############################################
################ EXECUTE #####################

mkdir shortread
cd shortread
bash $SHORTAUTO $ASMFOFN $LR1 $LR2 $R1 $R2 $EvalThese $SHORTSCRIPTS
cd ../

mkdir bionano
cd bionano
$BIONANORUN $BIONANOCLEAN $BIONANOCONFIG $ASMFOFN $MAPSFOFN $REC_ENZ $REC_SEQ $BIONANOSCRIPTS
cd ../

mkdir longread
cd longread 
bash $AUTOLR $LRCLEAN $LRCONFIG $ASMFOFN $LRSCRIPTS $ONT $PACBIO $ONT1 $ONT2 $PACBIO1 $PACBIO2
cd ../


mkdir blast_analyses
cd blast_analyses

mkdir transcriptome
cd transcriptome
$TRANSRUN $TRANSSCRIPTS $TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS1 $TRANSNJOBS $TBLASTX $TRANJOBPRE1
cd ../

mkdir dmel
cd dmel
$TRANSRUN $TRANSSCRIPTS $TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS2 $TRANSNJOBS $TBLASTX $TRANJOBPRE2
cd ../

mkdir anopheles
cd anopheles
$TRANSRUN $TRANSSCRIPTS $TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS2 $TRANSNJOBS $TBLASTX $TRANJOBPRE3
cd ../

mkdir dmel_peptides
cd dmel_peptides
$PEPRUN $PEPSCRIPTS $PEPCONFIG $PEPCLEAN $ASMFOFN $PEP2 $PEPNJOBS $PEPJOBPRE2
cd ../

mkdir anopheles_peptides
cd anopheles_peptides
$PEPRUN $PEPSCRIPTS $PEPCONFIG $PEPCLEAN $ASMFOFN $PEP2 $PEPNJOBS $PEPJOBPRE3
cd ../

mkdir knownseqs
cd knownseqs
$TRANSRUN $TRANSSCRIPTS $TRANSCONFIG $KNOWNCLEAN $ASMFOFN $KNOWNSEQS $KNOWNNJOBS $KNOWNTBLASTX $KNOWNJOBPRE
cd ../

#leave blast_analyses
cd ../

mkdir rnaseq
cd rnaseq
$RNARUN $RNASCRIPTS $RNACONFIG $RNACLEAN $ASMFOFN $RNAFOFN
cd ../

################ EXECUTE #####################
##############################################
##############################################
##############################################





