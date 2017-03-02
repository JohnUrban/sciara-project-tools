#!/bin/bash

for f in */merge/; do
 a=`echo $f | sed 's/\/merge\///'`
 score=`cat $f/score.txt`
 cov=`cat $f/total_base_cov.txt` ##WARNING
 span=`mergeBed -i $f/*bedGraph | awk '{s+=$3-$2}END{print s}'`
 numaln=`cat $f/num_alignments.txt`
 echo -e $a"\t"$score"\t"$cov"\t"$span"\t"$numaln
done
