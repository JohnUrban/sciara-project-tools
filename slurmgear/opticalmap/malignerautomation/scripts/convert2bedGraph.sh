#!/bin/bash

#need: ASM,BASE,ALL

G=${BASE}.genome
faSize -detailed $ASM > $G    ###{BASE}.genome

N=`tail -n +2 ${ALL} | wc -l | awk '{print $1}'`

if [ $N -eq 0 ]; then
  echo 0 > span.txt
  echo 0 > total_base_cov.txt
else
  tail -n +2 $ALL | awk 'OFS="\t" { if ($10<0) print $2,0,$11; else print $2,$10,$11}' | sortBed -i - | genomeCoverageBed -i - -g ${BASE}.genome -bg > ${ALL}.bedGraph
  awk '{s+=($3-$2)}END{print s}' ${ALL}.bedGraph > span.txt
  awk '{s+=($3-$2)*$4}END{print s}' ${ALL}.bedGraph > total_base_cov.txt
fi
## NOTE: total_base_cov was in error --	($3-$2)*4 instead of ($3-$2)*$4


## COMBINE ALL SCORES

## calculate metrics
score=`head -n 1 score.txt`
span=`cat span.txt`
cov=`cat total_base_cov.txt`
num=`cat num_alignments.txt`

#Score/Cov ... was supposed to represent an avg score per covered base, but this metric does not seem to do well in practice -- mostly abandonded
if [ $cov -eq 0 ]; then scorecov=0; else scorecov=`python -c "print 1e4*$score/$cov.0"`; fi

#Score/Num = Avg Score
if [ $num -eq 0 ]; then scorenum=0; else scorenum=`python -c "print $score/$num.0"`; fi

## ASM size
asmsize=`awk '{s+=$2}END{print s}' $G`

## Span/AsmSize = proportion asm spanned by alignments
spanasm=`python -c "print $span/$asmsize.0"`

## Cov/AsmSize = Coverage normalized to AsmSize = expected cov over any given base in assembly
covasm=`python -c "print $cov/$asmsize.0"`

## Cov/Span = Coverage normalized to span length = expected coverage over bases in regions that are actually covered
if [ $span -eq 0 ]; then covspan=0; else covspan=`python -c "print $cov/$span.0"`; fi

## Cov/Num = Total amount of bases covered normalized to number of alignments = avg number of bases covered per alignment (proportional to avg mol length that was aligned)
if [ $num -eq 0 ]; then covnum=0; else covnum=`python -c "print $cov/$num.0"`; fi

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
