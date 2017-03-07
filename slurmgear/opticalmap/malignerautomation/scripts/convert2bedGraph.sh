#!/bin/bash

#need: ASM,BASE,ALL

G=${BASE}.genome
faSize -detailed $ASM > $G    ###{BASE}.genome

tail -n +2 $ALL | awk 'OFS="\t" { if ($10<0) print $2,0,$11; else print $2,$10,$11}' | sortBed -i - | genomeCoverageBed -i - -g ${BASE}.genome -bg > ${ALL}.bedGraph
awk '{s+=($3-$2)}END{print s}' ${ALL}.bedGraph > span.txt
awk '{s+=($3-$2)*$4}END{print s}' ${ALL}.bedGraph > total_base_cov.txt
## NOTE: total_base_cov was in error --	($3-$2)*4 instead of ($3-$2)*$4


## COMBINE ALL SCORES

## calculate metrics
score=`cat score.txt`
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
echo -e Cov/Span"\t"$covspan >>$1/allstats.txt
echo -e Cov/Num"\t"$covnum >> allstats.txt
