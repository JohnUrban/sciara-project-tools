#!/bin/bash

source /users/jurban/data/software/hisat/source.sh

echo P, ${P}
echo READSFOFN, $READSFOFN
echo HIDX, $HIDX
echo PRE, $PRE
echo CLEAN, $CLEAN


## PROCESS READSFOFN TO GET READLISTS
NLINES=`cat $READSFOFN | wc -l`
i=0
while read line; do
  let i++
  if [ $i -eq $NLINES ]; then
    R1LIST+=`echo $line | awk '{print $1}'`;
    R2LIST+=`echo $line | awk '{print $2}'`;
  else
    R1LIST+=`echo $line | awk '{print $1","}'`;
    R2LIST+=`echo $line | awk '{print $2","}'`;
  fi
done < $READSFOFN

echo R1LIST, $R1LIST
echo R2LIST, $R2LIST
echo



## MAP READs, COLLECT STATS, NO INTERMEDIATE FILES
#hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS --reorder 2> mapreads.err | samtools view -bSh -q 0 -F 256 - | bedtools bamtobed -bedpe -i - 2> bamtobed.err | awk 'OFS="\n" {diff+=$1!=$4; same+=$1==$4; diff2+=$1!=$4&&$8>=2; atleastonemapped+=$1!="."||$4!="."; bothmapped+=$1!="."&&$4!="."; bothmappeddiff+=$1!="."&&$4!="."&&$1!=$4; total+=1; MAPQ+=$8}END{print "same_contig\t"same, "different_contig\t"diff, "different_contig_both_mapped\t"bothmappeddiff, "different_contig_qge2\t"diff2, "only_one_mapped\t"atleastonemapped-bothmapped, "both_mapped\t"bothmapped, "total_pairs\t"total, "pct_pairs_same_ctg\t"100.0*same/total, "pct_pairs_not_same_ctg\t"100.0*diff/total, "pct_pairs_with_1-2mappedmates_not_same_ctg\t"100.0*diff/(atleastonemapped), "pct_pairs_bothmapped_not_same_ctg\t"100.0*bothmappeddiff/bothmapped, "pct_pairs_not_same_ctg_qge2\t"100.0*diff2/total, "pct_pairs_with_2mappedmates_not_same_ctg_qge2\t"100.0*diff/(bothmapped), "sum_min_mapq\t"MAPQ, "avg_min_mapq_all_pairs\t"MAPQ/total, "avg_min_mapq_pairs_bothmatesmapped\t"MAPQ/bothmapped}' 

## OK - Making intermediate BAM so that I can have option of keeping it
hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS --reorder 2> mapreads.err > ${PRE}.bam
samtools view -bSh -q 0 -F 256 ${PRE}.bam | bedtools bamtobed -bedpe -i - 2> bamtobed.err | awk 'OFS="\n" {diff+=$1!=$4; same+=$1==$4; diff2+=$1!=$4&&$8>=2; atleastonemapped+=$1!="."||$4!="."; bothmapped+=$1!="."&&$4!="."; bothmappeddiff+=$1!="."&&$4!="."&&$1!=$4; total+=1; MAPQ+=$8}END{print "same_contig\t"same, "different_contig\t"diff, "different_contig_both_mapped\t"bothmappeddiff, "different_contig_qge2\t"diff2, "only_one_mapped\t"atleastonemapped-bothmapped, "both_mapped\t"bothmapped, "total_pairs\t"total, "pct_pairs_same_ctg\t"100.0*same/total, "pct_pairs_not_same_ctg\t"100.0*diff/total, "pct_pairs_with_1-2mappedmates_not_same_ctg\t"100.0*diff/(atleastonemapped), "pct_pairs_bothmapped_not_same_ctg\t"100.0*bothmappeddiff/bothmapped, "pct_pairs_not_same_ctg_qge2\t"100.0*diff2/total, "pct_bothmapped_not_same_ctg_qge2\t"100.0*diff2/bothmapped, "pct_bothmappeddiff_with_qge2\t"100.0*diff2/(bothmappeddiff), "sum_min_mapq\t"MAPQ, "avg_min_mapq_all_pairs\t"MAPQ/total, "avg_min_mapq_pairs_bothmatesmapped\t"MAPQ/bothmapped, "overall_alignment_rate\t"100.0*(bothmapped*2+(atleastonemapped-bothmapped))/(total*2)}' > info.txt
conc1=`grep "aligned concordantly exactly 1 time" mapreads.err | awk '{print $1}'`
conc2=`grep "aligned concordantly >1 times" mapreads.err | awk '{print $1}'`
conc=`echo $conc1 $conc2 | awk '{print $1+$2}'`
disc=`grep "aligned discordantly 1 time" mapreads.err | awk '{print $1}'`
nbothmapped=`grep -w both_mapped info.txt | awk '{print $2}'`
echo -e concordant_pairs"\t"${conc} >> info.txt
echo -e discordant_pairs"\t"${disc} >> info.txt
echo pct_of_mapped_pairs_that_are_concordant ${conc} ${nbothmapped} | awk 'OFS="\t" {print $1, 100.0*$2/$3}' >> info.txt
echo pct_of_mapped_pairs_that_are_discordant ${disc} ${nbothmapped} | awk 'OFS="\t" {print $1, 100.0*$2/$3}' >> info.txt

# CLEAN
if $CLEAN; then
  rm ${PRE}.bam;
  rm ../ht2/*
fi


############################################################################
## CAN IGNORE ALL BELOW
## ABOVE LINE CAPTURES ALL THAT I WANTED

## MAP READS
#hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS 2> mapreads.err | samtools view -bSh - | samtools sort -T ${PRE} > ${PRE}.bam
#samtools index ${PRE}.bam

## MAP READS
##hisat2 -p ${P} --dta -x ${HIDX} -1 $R1LIST -2 $R2LIST --rna-strandness $STRANDEDNESS 2> mapreads.err | samtools view -bSh - | samtools sort -n -T ${PRE} > ${PRE}.bam 


# MAP READS AND GET MAPQ STATS 
#echo num_reads sum_mapq avg_mapq | awk 'OFS="\t" {print $1,$2,$3}' > mapq.txt
#samtools view -F 4 ${PRE}.bam | awk 'OFS="\t" {s+=$5}END{print NR,s,s/NR}' >> mapq.txt
##TODO - Get AS too....
# GET INFO ON HOW MANY PAIRS MAP TO DIFFERENT CONTIGS
#echo same_ctg diff_ctg total_pairs pct_same pct_diff | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' > pairinfo.0.txt
#bedtools bamtobed -bedpe -i ${PRE}.bam 2> bamtobed.err | awk '{diff+=$1!=$4; same+=$1==$4; atleastonemapped+=$1!="."||$4!="."; bothmapped+=$1!="."&&$4!="."; total+=1}END{print same, diff, atleastonemapped-bothmapped, bothmapped, total, 100.0*same/total, 100.0*diff/total, 100.0*diff/(atleastonemapped)}'
##awk '{diff+=$1!=$4; same+=$1==$4; total+=1}END{print same, diff, total, 100.0*same/total, 100.0*diff/total}' >> pairinfo.0.txt
# GET INFO ON HOW MANY PAIRS MAP TO DIFFERENT CONTIGS - filtered
# This allows one to get rid of noise from multireads
# However, it also no longer counts pairs where 1 read mapped and 1 did not -- or where 1 read mapq > cutoff and other < cutoff
###for i in 2 10 20 30 40; do
#for i in 2 40; do
#  echo same_ctg diff_ctg total pct_same pct_diff | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' > pairinfo.${i}.txt
#  samtools view -bSh -q $i ${PRE}.bam | bedtools bamtobed -bedpe -i - 2> bamtobed.q${i}.err | awk '{diff+=$1!=$4; same+=$1==$4; atleastonemapped+=$1!="."||$4!="."; bothmapped+=$1!="."&&$4!="."; total+=1}END{print same, diff, atleastonemapped-bothmapped, bothmapped, total, 100.0*same/total, 100.0*diff/total, 100.0*diff/(atleastonemapped)}' >> pairinfo.${i}.txt &
#    #####awk '{diff+=$1!=$4; same+=$1==$4; total+=1}END{print same, diff, total, 100.0*same/total, 100.0*diff/total}' >> pairinfo.${i}.txt &
#done
#wait
## NEW WAY: one pass. Keep only primary alignments (-F 256). Keep all MAPQ - BEDPE uses lower MAPQ from the 2 mates in a pair. If one does not map, it will have MAPQ=0. 
#samtools view -bSh -q 0 -F 256 ${PRE}.bam
#awk 'OFS="\n" {diff+=$1!=$4; same+=$1==$4; diff2+=$1!=$4&&$8>=2; atleastonemapped+=$1!="."||$4!="."; bothmapped+=$1!="."&&$4!="."; bothmappeddiff+=$1!="."&&$4!="."&&$1!=$4; total+=1; MAPQ+=$8}END{print "same_contig\t"same, "different_contig\t"diff, "different_contig_both_mapped\t"bothmappeddiff, "different_contig_qge2\t"diff2, "only_one_mapped\t"atleastonemapped-bothmapped, "both_mapped\t"bothmapped, "total_pairs\t"total, "pct_pairs_same_ctg\t"100.0*same/total, "pct_pairs_not_same_ctg\t"100.0*diff/total, "pct_pairs_with_1-2mappedmates_not_same_ctg\t"100.0*diff/(atleastonemapped), "pct_pairs_bothmapped_not_same_ctg\t"100.0*bothmappeddiff/bothmapped, "pct_pairs_not_same_ctg_qge2\t"100.0*diff2/total, "pct_pairs_with_2mappedmates_not_same_ctg_qge2\t"100.0*diff/(bothmapped), "sum_min_mapq\t"MAPQ, "avg_min_mapq_all_pairs\t"MAPQ/total, "avg_min_mapq_pairs_bothmatesmapped\t"MAPQ/bothmapped}'


# CLEAN
#if $CLEAN; then 
#  rm ${PRE}.bam; 
#  rm ../ht2/*
#fi
