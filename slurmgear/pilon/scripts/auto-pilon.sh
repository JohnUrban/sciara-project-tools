#!/bin/bash

if [ $# -eq 0 ]; then echo "
Arg1=SUFFIX to remove from basename of assembly (e.g. 1x so it can be replaced with 2x) -- use \"\" if nothing
Arg2=NEWSUFFIX to add to modified basename
Arg3=Logical true/false -- clean as it works?
Arg4=Config file
Arg5=FOFN - file of filenames -- paths to each assembly FASTA to polish - extension for all fasta files should be .fasta
Arg6=BUILDBT2 logical
Arg7=MAPREADS logical
Arg8=FLAG1 logical
Arg9=FLAG2 logical
Arg10=MARK mark duplicates logical
Arg11=PILON logical
Arg12=FIX bases,local,all,etc
Arg13=NOSTRAYS logical
Arg14=Pilon pipeline location
Arg15=R1 location
Arg16=R2 location
Arg17 scripts loc
Example:
bash $AUTO $SFX2RM $SFX2ADD $CLEAN $CONFIG $ASMFOFN $BUILDBT2 $MAPREADS $FLAG1 $FLAG2 $MARK $PILON $FIX $NOSTRAYS $PILONPIPE $R1 $R2 $SCRIPTS
"; exit; fi

## DEBUG CODE
#for i in {1..16}; do
# echo -e $i"\t" ${!i}
#done
#exit

SUFFIX=$1
NEWSUFFIX=$2
CLEAN=$3
CONFIG=$4
FOFN=$5
BUILDBT2=$6
MAPREADS=$7
FLAG1=$8
FLAG2=$9
MARK=${10}
PILON=${11}
FIX=${12}
NOSTRAYS=${13}
PILONPIPE=${14}
R1=${15}
R2=${16}
SCRIPTS=${17}

## DEBUG CODE
#for v in SUFFIX NEWSUFFIX CLEAN FOFN BUILDBT2 MAPREADS FLAG1 FLAG2 MARK PILON FIX NOSTRAYS PILONPIPE R1 R2 SCRIPTS; do
# echo -e $v "\t" ${!v}
#done
#exit

i=0
while read f; do
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  ref=`readlink -f $f` 
  b=`basename $f ${SUFFIX}.fasta`; 
  echo $b; 
  if [ ! -d $b ]; then mkdir $b; fi
  cd $b;
  asmout=${b}${NEWSUFFIX}
  $PILONPIPE $ref $asmout $QOS $CLEAN $CONFIG $BUILDBT2 $MAPREADS $FLAG1 $FLAG2 $MARK $PILON $FIX $NOSTRAYS $R1 $R2 $SCRIPTS
#  echo $f $ref $b 
  cd ../
done < $FOFN


