#!/bin/bash

## Virtuall the same as run-evalwrapper.sh
## However, Brown CCV changed from CentOS6 and SLURM15 to RedHat7 and Slurm17 -- and, most importantly, changed how to use condos (--qos vs --account), and what condos I can use.
## All of the "auto" and "pipeline" scripts have small modifications and were given the same name with "-new" appended.

## Input Assemblies FOFN
ASMFOFN=input.fofn
CLEANALL=true

## GETTING ABS PATH OF ASMFOFN
ASMFOFN=`readlink -f $ASMFOFN`

RENAME=true

## MIN CONTIG SIZE -- defaults to analyzing all (minlen=0).
MINLEN=0


##TODO:
## ADD OPTION: make all sequences uppercase


###################  SHORT READ  ###########################
# specify paths to lap read sample (LR1,LR2) and all reads (R1,R2)-- give dummy answers if will not be using (that will serve as place-holders)
LR1=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.1.fastq
LR2=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.2.fastq
R1=~/data/scratch/male-ilmn/data/ilmnraw/R1.fastq
R2=~/data/scratch/male-ilmn/data/ilmnraw/R2.fastq

## OPTIONS for what programs to use. 
## FILL IN   ==>"EvalThese"<==   BELOW WITH ONE OF THESE
ALL=eval.cfg
ALLDIRTY=eval.noclean.cfg
OnlyAle=eval.aleonly.cfg
OnlyBusco=eval.buscoOnly.cfg
OnlyLap=eval.laponly.cfg
OnlyReapr=eval.reapronly.cfg
OnlyReaprNoClean=eval.reapronly.noclean.cfg
OnlyReaprNoCleanAggressive=eval.reapronly.noclean.aggressive.cfg
OnlyPilon=pilon.eval.cfg
AleModule=alemodule.eval.cfg ## builds bt2, maps reads, runs ale


## FILL IN WITH CORRECT VARIABLE
EvalThese=$ALL

## May need to adjust the following
SHORTSCRIPTS=/gpfs_home/jurban/software/sciaratools/sciara-project-tools/slurmgear/shortreadeval/scripts/
SHORTAUTO=${SHORTSCRIPTS}/auto-shortreadeval-new.sh
SHORTEVAL=${SHORTSCRIPTS}/eval.ARGS-new.sh
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
BIONANORUN=${BIONANOSCRIPTS}/auto-malign-new.sh




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
AUTOLR=${LRSCRIPTS}/auto-lrpipe-new.sh
LRCONFIGS=${LRBASE}/configs/
LR_DEFAULT_CFG=longread-config-sciara.cfg
LR_ALE_MODULE_CFG=ale-module-longread-config-sciara.cfg

LRCONFIG=${LRCONFIGS}/${LR_DEFAULT_CFG}


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
OTHERSPP_TRANSCONFIG=${TRANSCONFIGS}/other-spp-trans-config-sciara.cfg
TRANSFOFN=${TRANSFOFNS}/
TRANSRUN=${TRANSSCRIPTS}/auto-trans-new.sh

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
PEPRUN=${PEPSCRIPTS}/auto-pep-new.sh

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

KNOWNCONFIG=${TRANSCONFIGS}/known-config-sciara.cfg 

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
RNARUN=${RNASCRIPTS}/auto-rnaseqeval-new.sh

RNAFOFN=`readlink -f $RNAFOFN`



############### BUSCO V3 SECTION #################
BUSCOV3CLEAN=true
BUSCOV3BASE=/users/jurban/software/sciaratools/sciara-project-tools/slurmgear/buscov3
BUSCOV3SCRIPTS=${BUSCOV3BASE}/scripts/
BUSCOV3CONFIGS=${BUSCOV3BASE}/configs/
BUSCOV3FOFNS=${BUSCOV3BASE}/fofns/
BUSCOV3CONFIG=${BUSCOV3CONFIGS}/buscov3-config-sciara.cfg
BUSCOV3RUN=${BUSCOV3SCRIPTS}/auto-buscov3-new.sh





##############################################
##############################################
##############################################
################ EXECUTE #####################
NEWDIR=eval_asms
if $RENAME || [ $MINLEN -gt 0 ]; then mkdir -p ${NEWDIR} ; fi

## IF OPTED; rename and/or set min contig length
if $RENAME; then echo renaming....;
 if [ $MINLEN -gt 0 ]; then echo ...also setting min contig length to $MINLEN ; fi
 while read fasta; do
  b=`basename $fasta`
  fasta_name_changer.py -f $fasta -r contig -n --key key-${b}.txt | extractFastxEntries.py --fa --stdin --minlen $MINLEN > ${NEWDIR}/${b}
 done < $ASMFOFN
elif [ $MINLEN -gt 0 ] && [ $RENAME != "true" ]; then echo ...setting min contig length to $MINLEN ... ;
 while read fasta; do
  b=`basename $fasta`
  extractFastxEntries.py --fa -f $fasta --minlen $MINLEN > ${NEWDIR}/${b}
 done < $ASMFOFN
fi

## CREATE AND USE UPDATED FOFN
if $RENAME || [ $MINLEN -gt 0 ]; then
  for f in ${NEWDIR}/*; do
    readlink -f $f;
  done > renamed.fofn
  ASMFOFN=`readlink -f renamed.fofn`
fi

#for f in ${NEWDIR}/*; do
#  readlink -f $f; 
#done > renamed.fofn
#ASMFOFN=`readlink -f renamed.fofn`


## BEGIN LAUNCHING JOBS
echo shortread
mkdir -p shortread
cd shortread
##bash $SHORTAUTO $ASMFOFN $LR1 $LR2 $R1 $R2 $EvalThese $SHORTSCRIPTS
bash $SHORTAUTO $ASMFOFN $LR1 $LR2 $R1 $R2 $SHORTCONFIG $SHORTSCRIPTS
cd ../

echo bionano
mkdir -p bionano
cd bionano
$BIONANORUN $BIONANOCLEAN $BIONANOCONFIG $ASMFOFN $MAPSFOFN $REC_ENZ $REC_SEQ $BIONANOSCRIPTS
cd ../

echo longread
mkdir -p longread
cd longread 
bash $AUTOLR $LRCLEAN $LRCONFIG $ASMFOFN $LRSCRIPTS $ONT $PACBIO $ONT1 $ONT2 $PACBIO1 $PACBIO2
cd ../

echo blast_analyses
mkdir -p blast_analyses
cd blast_analyses

echo transcriptome
mkdir -p transcriptome
cd transcriptome
$TRANSRUN $TRANSSCRIPTS $TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS1 $TRANSNJOBS $TBLASTX $TRANJOBPRE1
cd ../

echo dmel
mkdir -p dmel
cd dmel
$TRANSRUN $TRANSSCRIPTS $OTHERSPP_TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS2 $TRANSNJOBS $TBLASTX $TRANJOBPRE2
cd ../

echo anopheles
mkdir -p anopheles
cd anopheles
$TRANSRUN $TRANSSCRIPTS $OTHERSPP_TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS3 $TRANSNJOBS $TBLASTX $TRANJOBPRE3
cd ../

echo dmel_peptides
mkdir -p dmel_peptides
cd dmel_peptides
$PEPRUN $PEPSCRIPTS $PEPCONFIG $PEPCLEAN $ASMFOFN $PEP2 $PEPNJOBS $PEPJOBPRE2
cd ../

echo anopheles_peptides
mkdir -p anopheles_peptides
cd anopheles_peptides
$PEPRUN $PEPSCRIPTS $PEPCONFIG $PEPCLEAN $ASMFOFN $PEP3 $PEPNJOBS $PEPJOBPRE3
cd ../

echo knownseqs
mkdir -p knownseqs
cd knownseqs
$TRANSRUN $TRANSSCRIPTS $KNOWNCONFIG $KNOWNCLEAN $ASMFOFN $KNOWNSEQS $KNOWNNJOBS $KNOWNTBLASTX $KNOWNJOBPRE
cd ../

#leave blast_analyses
cd ../

echo rnaseq
mkdir -p rnaseq
cd rnaseq
$RNARUN $RNASCRIPTS $RNACONFIG $RNACLEAN $ASMFOFN $RNAFOFN
cd ../

echo buscov3
mkdir -p buscov3
cd buscov3
$BUSCOV3RUN $BUSCOV3SCRIPTS $BUSCOV3CONFIG $BUSCOV3CLEAN $ASMFOFN 
cd ../

################ EXECUTE #####################
##############################################
##############################################
##############################################





