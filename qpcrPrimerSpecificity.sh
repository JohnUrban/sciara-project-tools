#!/bin/bash
BT2=$1
P1=$2
P2=$3
MAXINS=1000

if [ $# -eq 0 ]; then
  echo 
  echo Usage: 
  echo primerSpecificity.sh bt2index primer1.fa primer2.fa [optional: blastn_index genome_sequence.fa]
  echo
  echo primer1.fa can have all the fwd primers for multiple pairs.
  echo primer2.fa can have all the corresponding rev primers for the same primer pairs.
  echo The primer pairs should be on the same line in each file -- same exact as paired reads.
  echo If blastn index and genome.fa given as 4th and 5th args, it will also spit out info on BLAST alignment of the target.
  echo
  exit 0
fi

## could use "--no-head" instead of "grep -v ^@"
## could use "--omit-sec-seq" to omit sequence and quals (instead of cutting them out)

echo "MAPPING PRIMER 1 INDEPENDENTLY"
A=`bowtie2 -x $BT2 -U $P1 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
echo
echo "Primer 1 maps to" $A "site(s) independently."
echo

echo "MAPPING PRIMER 2 INDEPENDENTLY"
B=`bowtie2 -x $BT2 -U $P2 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
echo
echo "Primer 2 maps to" $B "site(s) independently."
echo

echo "MAPPING PRIMER 1 INDEPENDENTLY - LESS STRINGENT"
AL=`bowtie2 --rdg 2,2 --rfg 2,2 --mp 2 --very-sensitive -N 1 -x $BT2 -U $P1 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
echo
echo "Primer 1 maps to" $AL "site(s) independently - less stringently."
echo

echo "MAPPING PRIMER 2 INDEPENDENTLY - LESS STRINGENT"
BL=`bowtie2 --rdg 2,2 --rfg 2,2 --mp 2 --very-sensitive -N 1 -x $BT2 -U $P2 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
echo
echo "Primer 2 maps to" $BL "site(s) independently - less stringently."
echo



echo "STRINGENT PAIRING TEST"
AA=`bowtie2 --no-discordant --no-mixed -x $BT2 -1 $P1 -2 $P2 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
AA=`echo $AA/2 | bc`
echo
echo "There is/are" $AA "stringent site(s) with <10kb between them for this primer pair."
echo

echo "LESS STRINGENT PAIRING TEST"
BB=`bowtie2 --no-discordant --no-mixed --rdg 2,2 --rfg 2,2 --mp 2 --very-sensitive -N 1 -x $BT2 -1 $P1 -2 $P2 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
BB=`echo $BB/2 | bc`
echo
echo "There is/are" $BB "total site(s) with <10kb between them when requiring less stringent primer:template matching for this primer pair."
echo

echo "PAIRING TEST - PRIMER 1 WITH ITSELF"
CC=`bowtie2 --no-discordant --no-mixed -x $BT2 -1 $P1 -2 $P1 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
CC=`echo $CC/2 | bc`
echo
echo "There is/are" $CC "total site(s) with <10kb between them for this primer pair."
echo

echo "PAIRING TEST - PRIMER 2 WITH ITSELF"
DD=`bowtie2 --no-discordant --no-mixed -x $BT2 -1 $P2 -2 $P2 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
DD=`echo $DD/2 | bc`
echo
echo "There is/are" $DD "total site(s) with <10kb between them for this primer pair."
echo

echo "LOW STRINGENCY PAIRING TEST - PRIMER 1 WITH ITSELF"
CCL=`bowtie2 --no-discordant --no-mixed --rdg 2,2 --rfg 2,2 --mp 2 --very-sensitive -N 1 -x $BT2 -1 $P1 -2 $P1 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
CCL=`echo $CCL/2 | bc`
echo
echo "There is/are" $CCL "total site(s) with <10kb between them when requiring less stringent primer:template matching for this primer pair."
echo

echo "LOW STRINGENCY PAIRING TEST - PRIMER 2 WITH ITSELF"
DDL=`bowtie2 --no-discordant --no-mixed --rdg 2,2 --rfg 2,2 --mp 2 --very-sensitive -N 1 -x $BT2 -1 $P2 -2 $P2 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | cut -f 1,2,3,4,5,6,7,9,12,13,14,15,16,17,18,19,20 | tee /dev/stderr | wc -l`
DDL=`echo $DDL/2 | bc`
echo
echo "There is/are" $DDL "total site(s) with <10kb between them when requiring less stringent primer:template matching for this primer pair."
echo

if [ $# -eq 5 ]; then
 echo BLASTING PCR TARGET
 DIR=`dirname $2`
 BASE=`basename $2 | stringEdit - .fa`
 TARGET=${DIR}/${BASE}-target.fa
 bowtie2 --no-discordant --no-mixed -x $1 -1 $2 -2 $3 -f -a --maxins $MAXINS | samtools view -S -F 4 - | grep -v ^@ | grep ^p1 | awk 'OFS="\t" {print $3,$4-1,$4-1+$9}' | fastaFromBed -fi $5 -bed - -fo $TARGET
 echo
 BLAST=`blastn -db $4 -query $TARGET -outfmt 6 -culling_limit 1 | tee /dev/stderr | wc -l`
 echo
fi

echo FINAL REPORT:
echo INDEPENDENT STATS:
echo "Primer 1 maps to" $A "site(s) independently."
echo "Primer 1 maps to" $AL "site(s) independently - when requiring less stringency."
echo "Primer 2 maps to" $B "site(s) independently."
echo "Primer 2 maps to" $BL "site(s) independently - when requiring less stringency."

echo
echo PROPER PAIRING:
echo "There is/are" $AA "stringent site(s) with <10kb between them for this primer pair."
echo "There is/are" $BB "total site(s) with <10kb between them when requiring less stringent primer:template matching for this primer pair."
echo 
echo SELF PAIRING - PRIMER 1:
echo "There is/are" $CC "total low stringency site(s) with <10kb between them for this primer pair."
echo "There is/are" $CCL "total site(s) with <10kb between them when requiring less stringent primer:template matching for this primer pair."
echo
echo SELF PAIRING - PRIMER 2:
echo "There is/are" $DD "total low stringency site(s) with <10kb between them for this primer pair."
echo "There is/are" $DDL "total site(s) with <10kb between them when requiring less stringent primer:template matching for this primer pair."
echo
if [ $# -eq 5 ]; then
 echo "PCR TARGET BLAST RESULTS"
 echo The target PCR product matches $BLAST sites in the genome sequence given.
fi
