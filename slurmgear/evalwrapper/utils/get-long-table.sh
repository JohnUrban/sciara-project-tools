#!/bin/bash


if [ $# -eq 0 ] || [ $1 == "-h" ] || [ $1 == "-help" ] || [ $1 == "--help" ]; then
    echo "
    Usage: bash $0 FOFN

    ...where FOFN has list of all assemblies used in the assembly evaluations in subdirs you are trying to summarize.
    (( typically called input.fofn ))

    Alt Usage: bash $0 FOFN debug
    ...writing the word 'debug' as the second argument will create files that count the number of metrics from each function in 2-col tab-delim file.

    Alt Usage: bash $0 FOFN vizmat
    ...writing the word 'vizmat' as the second argument will create longtable files that have been pruned compared to default and are strictly numeric.
    ...this one may ultimately be used for LaTeX table conversion instead of default if desired as well....

    Alt Usage: bash $0 FOFN simple
    ...writing the word 'simple' as the second argument will create longtable files that have been pruned compared to both default and vizmat.
    ...like vizmat these are strictly numeric.
    ...only metrics that can easily be interpreted in isolation as better or worse when compared to another assembly were retained
    ...this one can be used for LaTeX table conversion BUT it is not recommended as a part of MANUAL INSPECTION
    ...the default or vizmat provide information that you can process as a human more readily by combining with other info to judge as better or worse or curious
    ...the point of simple is to be able to compare 2 asms and plot -1 or 1 (later on) for better or worse.
    ...use 'simplenames' to get names; use 'simplefactors' to get the multiplication vector to use against each pair of metric scores before ranking (in R)

    Alt Usage: bash $0 FOFN all
    ...writing the word 'all' as the second argument will create the default longtables, the vizmat longtables, and the debug files in one pass.

    Alt Usage: bash $0 FOFN names|vizmatnames|simplenames|simplefactors
    ...writing the word 'names' or 'vizmatnames' or 'simplenames' as the second argument will give a comma-sep list of metric names in the order that they appear in the corresponding long table form.
    ...writing the word 'simplefactors' returns a vector of 1 and -1 to be used to multiply against a vector of scores for the given metric (same order as simple) before ranking in R
    ...for now you still need a dummy word in position 1 -- e.g. 'input.fofn' -- but it will not use it.
    "
    exit
fi

##FOFN=input.fofn
FOFN=$1
NARG=$#
DEBUG=''
if [ $NARG -eq 2 ]; then DEBUG=$2; fi

## ASSUMES FOLLOWING DIRS
SHORT=shortread
BIONANO=bionano
LONG=longread


##NAME VARIABLES
VIZMATNAMES="Num_Contigs,Max_Contig_Size,NG50,LG50,E_G,Complete_BUSCOs,Complete_and_single-copy_BUSCOs,Complete_and_duplicated_BUSCOs,Fragmented_BUSCOs,Missing_BUSCOs,Bowtie2_Ilmn_Pct,Bowtie2_Pct_Aligned_Concordantly_Once,Bowtie2_Pct_Aligned_Concordantly_Multiple,Bowtie2_Pct_Did_NOT_Align_Concordantly,Bowtie2_Pct_DNAC_that_align_Discordantly,Bowtie2_Pct_align_discordantly,Bowtie2_Pct_didnt_align_Concord_nor_Discord,Bowtie2_Pct_SE_reads_from_DNACD_pairs_that_did_not_align,Bowtie2_Pct_SE_reads_from_DNACD_pairs_that_aligned_once,Bowtie2_Pct_SE_reads_from_DNACD_pairs_that_aligned_multiple,ALE_Score,ALE_Place_Avg,ALE_Insert_Avg,ALE_kmer_Avg,ALE_Depth_Score_Avg,ALE_Total_Mapped_Reads,ALE_Total_UnMapped_Reads,ALE_Total_Placed_Reads,LAP,REAPR_MBS,REAPR_Pct_EF,REAPR_Num_Errors,REAPR_FCD_Errors_within_a_contig,REAPR_FCD_Errors_over_a_gap,REAPR_Low_frag_cov_within_a_contig,REAPR_Low_frag_cov_over_a_gap,REAPR_Num_Warnings,REAPR_Low_Score_Regions,REAPR_Links,REAPR_Soft_Clip,REAPR_Collapsed_repeats,REAPR_Low_Read_Coverage,REAPR_Low_Perfect_Coverage,REAPR_Wrong_Orientation,REAPR_Broken_Num_Contigs,REAPR_Broken_Max_Contig_Length,REAPR_Broken_NG50,REAPR_Broken_LG50,REAPR_Broken_EG,REAPR_Broken_Num_Gaps,REAPR_Broken_Total_Gap_Len,FRC_Num_Features,FRC_Norm_Num_Features_(per_Mb),FRC_Mapped,FRC_Unmapped,FRC_Proper,FRC_Wrong_Dist,FRC_Wrong_Orientation,FRC_Wrong_Contig,FRC_Singleton,FRC_Mean_Cov,FRC_Spanning_Cov,FRC_Proper_Pairs_Cov,FRC_Wrong_Mate_Cov,FRC_Singleton_Mate_Cov,FRC_Different_Contig_Cov,BNG_Score,BNG_Span,BNG_Cov,BNG_Num,BNG_1e4*Score/Cov,BNG_Score/Num,BNG_Span/AsmSize,BNG_Cov/AsmSize,BNG_Cov/Span,BNG_Cov/Num,ONT_num_align,ONT_num_entries,ONT_num_uniq_align,ONT_num_uniq_entries,ONT_sum_mapq,ONT_pct_aln,ONT_ratio,ONT_avg_mapq,PacBio_num_align_,PacBio_num_entries,PacBio_num_uniq_align,PacBio_num_uniq_entries,PacBio_sum_mapq,PacBio_Pct_aln,PacBio_ratio,PacBio_avg_mapq,Sniffles_ONT_Num_SV,Sniffles_PacBio_Num_SV,Sniffles_Combined_Num_SV,Sniffles_ONT_Sum_SV,Sniffles_PacBio_Sum_SV,Sniffles_Combined_Sum_SV"
NAMES="Num_Contigs,Asm_Size,Max_Contig_Size,NG50,LG50,E_G,Complete_BUSCOs,Complete_and_single-copy_BUSCOs,Complete_and_duplicated_BUSCOs,Fragmented_BUSCOs,Missing_BUSCOs,Bowtie2_Ilmn_Pct,Bowtie2_Num_Paired,Bowtie2_Aligned_Concordantly_Once,Bowtie2_Aligned_Concordantly_Multiple,Bowtie2_Did_NOT_Align_Concordantly,Bowtie2_Percent_DNAC_that_align_Discordantly,Bowtie2_Number_and_Percent_Total_that_align_discordantly,Bowtie2_Pairs_didnt_align_Concord_nor_Discord,Bowtie2_Num_SE_reads_in_those_pairs,Bowtie2_Num_SE_reads_in_those_that_did_not_align,Bowtie2_Num_SE_reads_in_those_that_aligned_once,Bowtie2_Num_SE_reads_in_those_that_aligned_multiple,ALE_Score,ALE_Place_Avg,ALE_Insert_Avg,ALE_kmer_Avg,ALE_Depth_Score_Avg,ALE_Depth_Avg,ALE_Total_Reads,ALE_Total_Mapped_Reads,ALE_Total_UnMapped_Reads,ALE_Total_Placed_Reads,LAP,REAPR_MBS,REAPR_Pct_EF,REAPR_Num_Errors,REAPR_FCD_Errors_within_a_contig,REAPR_FCD_Errors_over_a_gap,REAPR_Low_frag_cov_within_a_contig,REAPR_Low_frag_cov_over_a_gap,REAPR_Num_Warnings,REAPR_Low_Score_Regions,REAPR_Links,REAPR_Soft_Clip,REAPR_Collapsed_repeats,REAPR_Low_Read_Coverage,REAPR_Low_Perfect_Coverage,REAPR_Wrong_Orientation,REAPR_Original_Asm_Length,REAPR_Broken_Asm_Length,REAPR_Original_Num_Contigs,REAPR_Broken_Num_Contigs,REAPR_Original_Mean_Contig_Length,REAPR_Broken_Mean_Contig_Length,REAPR_Original_Max_Contig_Length,REAPR_Broken_Max_Contig_Length,REAPR_Original_N50,REAPR_Broken_N50,REAPR_Original_NG50,REAPR_Broken_NG50,REAPR_Original_LG50,REAPR_Broken_LG50,REAPR_Original_EG,REAPR_Broken_EG,REAPR_Broken_Num_Gaps,REAPR_Broken_Total_Gap_Len,FRC_Num_Features,FRC_Norm_Num_Features_(per_Mb),FRC_Insert_Size_Mean,FRC_Insert_Size_Std_Dev,FRC_Reads,FRC_Mapped,FRC_Unmapped,FRC_Proper,FRC_Wrong_Dist,FRC_Zero_Qual,FRC_Wrong_Orientation,FRC_Wrong_Contig,FRC_Singleton,FRC_Mean_Cov,FRC_Spanning_Cov,FRC_Proper_Pairs_Cov,FRC_Wrong_Mate_Cov,FRC_Singleton_Mate_Cov,FRC_Different_Contig_Cov,BNG_Score,BNG_Span,BNG_Cov,BNG_Num,BNG_1e4*Score/Cov,BNG_Score/Num,BNG_Span/AsmSize,BNG_Cov/AsmSize,BNG_Cov/Span,BNG_Cov/Num,ONT_num_align,ONT_num_entries,ONT_num_uniq_align,ONT_num_uniq_entries,ONT_sum_mapq,ONT_pct_aln,ONT_ratio,ONT_avg_mapq,PacBio_num_align_,PacBio_num_entries,PacBio_num_uniq_align,PacBio_num_uniq_entries,PacBio_sum_mapq,PacBio_Pct_aln,PacBio_ratio,PacBio_avg_mapq,Sniffles_ONT_Num_SV,Sniffles_PacBio_Num_SV,Sniffles_Combined_Num_SV,Sniffles_ONT_Sum_SV,Sniffles_PacBio_Sum_SV,Sniffles_Combined_Sum_SV"
SIMPLENAMES="Num_Contigs,Max_Contig_Size,NG50,LG50,E_G,Complete_BUSCOs,Bowtie2_Ilmn_Pct,Bowtie2_Pct_Aligned_Concordantly,Bowtie2_Pct_Aligned_Disconcordantly,ALE_Score,ALE_Place_Avg,ALE_Insert_Avg,ALE_kmer_Avg,ALE_Depth_Score_Avg,LAP,REAPR_MBS,REAPR_Pct_EF,REAPR_FCD_Errors_within_a_contig,REAPR_FCD_Errors_over_a_gap,REAPR_Low_frag_cov_within_a_contig,REAPR_Low_frag_cov_over_a_gap,REAPR_Low_Score_Regions,REAPR_Links,REAPR_Soft_Clip,REAPR_Collapsed_repeats,REAPR_Low_Read_Coverage,REAPR_Low_Perfect_Coverage,REAPR_Wrong_Orientation,REAPR_Broken_Num_Contigs,REAPR_Broken_Max_Contig_Length,REAPR_Broken_NG50,REAPR_Broken_LG50,REAPR_Broken_EG,REAPR_Broken_Num_Gaps,REAPR_Broken_Total_Gap_Len,FRC_Num_Features,FRC_Norm_Num_Features_(per_Mb),FRC_Proper,FRC_Wrong_Dist,FRC_Wrong_Orientation,FRC_Wrong_Contig,FRC_Singleton,FRC_Mean_Cov,FRC_Spanning_Cov,FRC_Proper_Pairs_Cov,FRC_Wrong_Mate_Cov,FRC_Singleton_Mate_Cov,FRC_Different_Contig_Cov,BNG_Score,BNG_Span,BNG_Cov,BNG_Num,BNG_1e4*Score/Cov,BNG_Score/Num,BNG_Span/AsmSize,BNG_Cov/AsmSize,BNG_Cov/Span,BNG_Cov/Num,ONT_sum_mapq,ONT_pct_aln,ONT_ratio,ONT_avg_mapq,PacBio_sum_mapq,PacBio_Pct_aln,PacBio_ratio,PacBio_avg_mapq,Sniffles_ONT_Num_SV,Sniffles_PacBio_Num_SV,Sniffles_Combined_Num_SV,Sniffles_ONT_Sum_SV,Sniffles_PacBio_Sum_SV,Sniffles_Combined_Sum_SV"
SIMPLEFACTORS="1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,1,1,1"



function getsizestats { ##takes fasta
    b=`basename $1 .fasta`
    F=sizestats/${b}.tsv
    if [ ! -d sizestats ]; then mkdir sizestats; fi
    if [ ! -f $F ]; then faSize -detailed $1 | asm-stats.py -t | awk '{gsub(/,/,"\n"); print}' > $F ; fi
    awk 'NR==1 || NR==2 || NR==3 || NR==10 || NR==11 || NR==12 {print $1}' $F
}

function getsizestats_vizmat { ##takes fasta
    b=`basename $1 .fasta`
    F=sizestats/${b}.tsv
    if [ ! -d sizestats ]; then mkdir sizestats; fi
    if [ ! -f $F ]; then faSize -detailed $1 | asm-stats.py -t | awk '{gsub(/,/,"\n"); print}' > $F ; fi
    ## asmsize taken out b/c not meaningful in terms of an asm getting better or worse (Except in terms of going further away from expected size)
    awk 'NR==1 || NR==3 || NR==10 || NR==11 || NR==12 {print $1}' $F
}

function getbusco {
    F=$SHORT/${1}/busco/*/short*
    for l in "Complete BUSCOs" "Complete and single-copy BUSCOs" "Complete and duplicated BUSCOs" "Fragmented BUSCOs" "Missing BUSCOs"; do
        grep "${l}" $F | awk '{print $1}'
    done
}

function getbusco_simple {
    F=$SHORT/${1}/busco/*/short*
    ## Removed
    ## "Complete and single-copy BUSCOs" "Complete and duplicated BUSCOs" "Fragmented BUSCOs" "Missing BUSCOs"
    ## Since it is unclear w/o specific knowledge that having more or less of the complete buscos as single copy vs duplicated is "better" for an asm
    ## and since missing buscos gives no extra info over complete
    for l in "Complete BUSCOs"; do
        grep "${l}" $F | awk '{print $1}'
    done
}

function getbowtie2 {
    F=$SHORT/${1}/mreads/*.err
    grep "overall alignment rate" $F | awk '{sub(/%/,""); print $1}'
    pair=`grep "were paired; of these:" $F | awk '{print $1}'`
    echo $pair
    grep "aligned concordantly exactly 1 time" $F | awk '{print $1,$2}' | awk '{sub(/\ /,"_"); sub(/%/,"\\%"); print}'
    grep "aligned concordantly >1 times" $F | awk '{print $1,$2}' | awk '{sub(/\ /,"_"); sub(/%/,"\\%"); print}'
    grep "aligned concordantly 0 times" $F | grep -v "of these" | awk '{print $1,$2}' | awk '{sub(/\ /,"_"); sub(/%/,"\\%"); print}'
    grep "aligned discordantly 1 time" $F | awk '{print $2}' | awk '{sub(/\ /,"_"); sub(/%/,"\\%"); print}' #pct of pairs that aln conc 0 times that aln disc
    disc1=`grep "aligned discordantly 1 time" $F | awk '{print $1}'`
    pctdisc1=`echo $disc1 $pair | awk '{print 100.0*$1/$2}'`
    echo "${disc1}_(${pctdisc1}\%)" ## pct of all pairs that aln disc
    grep "pairs aligned 0 times concordantly or discordantly; of these:" $F | awk '{print $1}'
    grep "mates make up the pairs; of these:" $F | awk '{print $1}'
    grep "aligned 0 times" $F | awk '{print $1}'
    grep "aligned exactly 1 time" $F | awk '{print $1}'
    grep "aligned >1 times" $F | awk '{print $1}'
}

function getbowtie2_vizmat {
    F=$SHORT/${1}/mreads/*.err
    grep "overall alignment rate" $F | awk '{sub(/%/,""); print $1}'
    pair=`grep "were paired; of these:" $F | awk '{print $1}'`
    conc1=`grep "aligned concordantly exactly 1 time" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'`
    concmult=`grep "aligned concordantly >1 times" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'`
    echo $conc1
    echo $concmult
    #echo $conc1 $concmult | awk '{print $1+$2}   ## pct aligned concordantly
    grep "aligned concordantly 0 times" $F | grep -v "of these" | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'
    grep "aligned discordantly 1 time" $F | awk '{print $2}' | awk '{sub(/\ /,"_"); sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}' #pct of reads aln conc 0 times that aln disc
    disc1=`grep "aligned discordantly 1 time" $F | awk '{print $1}'`
    pctdisc1=`echo $disc1 $pair | awk '{print 100.0*$1/$2}'`
    echo ${pctdisc1} #pct of total pairs that aln disc
    unmap=`grep "pairs aligned 0 times concordantly or discordantly; of these:" $F | awk '{print $1}'`
    pctunmap=`echo $unmap $pair | awk '{print 100.0*$1/$2}'`
    echo $pctunmap
    grep "aligned 0 times" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'
    grep "aligned exactly 1 time" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'
    grep "aligned >1 times" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'
}


function getbowtie2_simple {
    ## 3 outputs: overall aln rate, %of all pairs conc, % of all pairs disc
    F=$SHORT/${1}/mreads/*.err
    grep "overall alignment rate" $F | awk '{sub(/%/,""); print $1}'
    pair=`grep "were paired; of these:" $F | awk '{print $1}'`
    conc1=`grep "aligned concordantly exactly 1 time" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'`
    concmult=`grep "aligned concordantly >1 times" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'`
    ## it is not immediately clear if more or less single vs multi is better -- easier to interpret %concordant marginalized over both
    echo $conc1 $concmult | awk '{print $1+$2}'   ## pct aligned concordantly
    ## removed pct DNAC since that does not give any new info wrt pct conc
    ## remove "pct of reads aln conc 0 times that aln disc" since this is not immediately interp as better/worse
    disc1=`grep "aligned discordantly 1 time" $F | awk '{print $1}'`
    pctdisc1=`echo $disc1 $pair | awk '{print 100.0*$1/$2}'`
    echo ${pctdisc1} #pct of total pairs that aln disc
    ## pct unmap does not nec give any more info 
    ## SE read are not immediately comparable/interpretable -- I will let ALE/LAP/REAPR/FRC deal with them
}



function getale {
    F=$SHORT/${1}/ale/*ALE.txt
    ## 1=ALE score, 2=ncont, 3=asmlen, 4=placeAvg, 5=insertAvg, 6=kmerAvg, 7=depthScoreAvg, 8=depthAvg, 9=totalreads, 10=totalMappedReads, 11=totalUnMappedReads, 12=total placed reads, 13=readlenAvg, 14=avgReadOvlap
    awk 'NR==1 || NR==4 || NR==5 || NR==6 || NR==7 || NR==8 || NR==9 || NR==10 || NR==11 || NR==12 {print $3}' $F
}

function getale_vizmat {
    F=$SHORT/${1}/ale/*ALE.txt
    ## 1=ALE score, 2=ncont, 3=asmlen, 4=placeAvg, 5=insertAvg, 6=kmerAvg, 7=depthScoreAvg, 8=depthAvg, 9=totalreads, 10=totalMappedReads, 11=totalUnMappedReads, 12=total placed reads, 13=readlenAvg, 14=avgReadOvlap
    awk 'NR==1 || NR==4 || NR==5 || NR==6 || NR==7 || NR==10 || NR==11 || NR==12 {print $3}' $F
}

function getale_simple {
    F=$SHORT/${1}/ale/*ALE.txt
    ## 1=ALE score, 2=ncont, 3=asmlen, 4=placeAvg, 5=insertAvg, 6=kmerAvg, 7=depthScoreAvg, 8=depthAvg, 9=totalreads, 10=totalMappedReads, 11=totalUnMappedReads, 12=total placed reads, 13=readlenAvg, 14=avgReadOvlap
    ## trimmed off mapped, unmapped, and placed -- will let those be reflected in BT2 fxn and as part of the other scores here.
    awk 'NR==1 || NR==4 || NR==5 || NR==6 || NR==7 {print $3}' $F
}

function getlap {
    F=$SHORT/${1}/lap/*.lapscore
    awk '{print $1}' $F | head -n 1
}


function getasmstats {
    F=sizestats/${b}.tsv
    paste -sd"\t" $F
}

function getreapr {
    D=$SHORT/${1}/reapr/output_directory/
    head -n 1 $D/per-base-mean-score.txt
    ## pctEF, numErrors, FCD, FCD_gap, frag_cov, frag_cov_gap, num warnings, lowscore, link, softclip, collapsed repeat, readcov, low perfect cov, readorientation
    awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39+$40+$41+$42, $39, $40, $41, $42, $43+$44+$45+$46+$47+$48+$49, $43, $44, $45, $46, $47, $48, $49}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'
    ## asm bases, brasm bases, nseqs, br nseqs, mean len, br meanlen, longest, br longest, n50, br n50
    awk 'NR==2 {OFS="\t"; print $2, $20, $3, $21, $4, $22, $5, $23, $6, $24}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
    stats=`getasmstats $1`
    bstats=`awk '{gsub(/,/,"\t"); print }' ${D}/broken_assembly.sizestats.csv`
    echo $stats | awk '{print $10}' ##ng50
    echo $bstats | awk '{print $10}' ##broken ng50
    echo $stats | awk '{print $11}' ## lg50
    echo $bstats | awk '{print $11}' ##broken lg50
    echo $stats | awk '{print $12}' ##eg
    echo $bstats | awk '{print $12}' ## broken eg
    ## ngaps, gaplen
    awk 'NR==2 {OFS="\t"; print $36, $37}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
}

function getreapr_vizmat {
    D=$SHORT/${1}/reapr/output_directory/
    head -n 1 $D/per-base-mean-score.txt
    ## pctEF, numErrors, FCD, FCD_gap, frag_cov, frag_cov_gap, num warnings, lowscore, link, softclip, collapsed repeat, readcov, low perfect cov, readorientation
    awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39+$40+$41+$42, $39, $40, $41, $42, $43+$44+$45+$46+$47+$48+$49, $43, $44, $45, $46, $47, $48, $49}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'
    ## br nseqs,  br longest, 
    awk 'NR==2 {OFS="\t"; print $21, $23}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
    stats=`getasmstats $1`
    bstats=`awk '{gsub(/,/,"\t"); print }' ${D}/broken_assembly.sizestats.csv`
    echo $bstats | awk '{print $10}' ##broken ng50
    echo $bstats | awk '{print $11}' ##broken lg50
    echo $bstats | awk '{print $12}' ## broken eg
    ## in broken asm: ngaps, gaplen
    awk 'NR==2 {OFS="\t"; print $36, $37}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
}

function getreapr_simple {
    ## FOR REAPR SIMPLE -- rm sumError and sumWarn -- just use indiv err and warns
    D=$SHORT/${1}/reapr/output_directory/
    head -n 1 $D/per-base-mean-score.txt
    ## pctEF, FCD, FCD_gap, frag_cov, frag_cov_gap, lowscore, link, softclip, collapsed repeat, readcov, low perfect cov, readorientation
    awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'
    ## br nseqs,  br longest, 
    awk 'NR==2 {OFS="\t"; print $21, $23}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
    stats=`getasmstats $1`
    bstats=`awk '{gsub(/,/,"\t"); print }' ${D}/broken_assembly.sizestats.csv`
    echo $bstats | awk '{print $10}' ##broken ng50
    echo $bstats | awk '{print $11}' ##broken lg50
    echo $bstats | awk '{print $12}' ## broken eg
    ## in broken asm: ngaps, gaplen
    awk 'NR==2 {OFS="\t"; print $36, $37}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
}



function getfrc {
    F1=$SHORT/${1}/frc/*gff
    F2=$SHORT/${1}/frc/*frc_assemblyTable.csv
    grep -c -v ^# $F1
    ##InsertSizeMean,InsertSizeStd,READS,MAPPED,UNMAPPED,PROPER,WRONG_DIST,ZERO_QUAL,WRONG_ORIENTATION,WRONG_CONTIG,SINGLETON,MEAN_COVERAGE,SPANNING_COVERAGE,PROPER_PAIRS_COVERAGE,WRONG_MATE_COVERAGE,SINGLETON_MATE_COV,DIFFERENT_CONTIG_COV
    tail -n 1 $F2 | awk '{gsub(/,/,"\n"); print}' | tail -n +3
}

function getfrc_vizmat {
    F1=$SHORT/${1}/frc/*gff
    F2=$SHORT/${1}/frc/*frc_assemblyTable.csv
    grep -c -v ^# $F1
    ## 1=BAM,2=LIB_TYPE,
    ##3=InsertSizeMean,4=InsertSizeStd,5=READS,
    ## 6=MAPPED,7=UNMAPPED,8=PROPER,9=WRONG_DIST,    --> 10=ZERO_QUAL, <--- excluding b/c usually 0
    ## 11=WRONG_ORIENTATION,12=WRONG_CONTIG,13=SINGLETON,    14=MEAN_COVERAGE, 
    ## 15=SPANNING_COVERAGE,16=PROPER_PAIRS_COVERAGE,17=WRONG_MATE_COVERAGE,18=SINGLETON_MATE_COV,19=DIFFERENT_CONTIG_COV
    tail -n 1 $F2 | awk 'OFS="\t" {gsub(/,/,"\t"); print $6,$7,$8,$9,$11,$12,$13,$14,$15,$16,$17,$18,$19}' | awk '{gsub(/\t/,"\n"); print}'
}

function getfrc_simple {
    ## FOR FRC SIMPLE -- rm 6,7, 10 
    ## was thnking about removing 13 and 14 -- but kept
    F1=$SHORT/${1}/frc/*gff
    F2=$SHORT/${1}/frc/*frc_assemblyTable.csv
    grep -c -v ^# $F1
    ## 1=BAM,2=LIB_TYPE,
    ##3=InsertSizeMean,4=InsertSizeStd,5=READS,
    ## 6=MAPPED,7=UNMAPPED,8=PROPER,9=WRONG_DIST,    --> 10=ZERO_QUAL, <--- excluding b/c usually 0
    ## 11=WRONG_ORIENTATION,12=WRONG_CONTIG,13=SINGLETON,    14=MEAN_COVERAGE, 
    ## 15=SPANNING_COVERAGE,16=PROPER_PAIRS_COVERAGE,17=WRONG_MATE_COVERAGE,18=SINGLETON_MATE_COV,19=DIFFERENT_CONTIG_COV
    tail -n 1 $F2 | awk 'OFS="\t" {gsub(/,/,"\t"); print $8,$9,$11,$12,$13,$14,$15,$16,$17,$18,$19}' | awk '{gsub(/\t/,"\n"); print}'
}


function calcmalignerstats {
    #Takes $D for merge directory as $1
    if [ ! -f $1/score.txt ]; then tail -n +2 $1/all.bionano.smoothed.maps.aln | awk '{s+=$19}END{print s}' > $1/score.txt; fi
    if [ ! -f $1/span.txt ]; then awk '{s+=($3-$2)}END{print s}' $1/all.bionano.smoothed.maps.aln.bedGraph > $1/span.txt; fi
    if [ ! -f $1/total_base_cov.txt ]; then awk '{s+=($3-$2)*$4}END{print s}' $1/all.bionano.smoothed.maps.aln.bedGraph > $1/total_base_cov.txt; fi
    if [ ! -f $1/num_alignments.txt ]; then tail -n +2 $1/all.bionano.smoothed.maps.aln | wc -l > $1/num_alignments.txt; fi
    G=$1/*.genome
    asmsize=`awk '{s+=$2}END{print s}' $G`

    ## calculate metrics
    score=`cat $1/score.txt`
    span=`cat $1/span.txt`
    cov=`cat $1/total_base_cov.txt`
    num=`cat $1/num_alignments.txt`
    scorecov=`python -c "print 1e4*$score/$cov.0"`
    scorenum=`python -c "print $score/$num.0"`
    asmsize=`awk '{s+=$2}END{print s}' $G`
    spanasm=`python -c "print $span/$asmsize.0"`
    covasm=`python -c "print $cov/$asmsize.0"`
    covspan=`python -c "print $cov/$span.0"`
    covnum=`python -c "print $cov/$num.0"`

    ## populate allstats.txt
    echo -e score"\t"$score > $1/allstats.txt
    echo -e span"\t"$span >> $1/allstats.txt
    echo -e total_cov"\t"$cov >> $1/allstats.txt
    echo -e num_aln"\t"$num >> $1/allstats.txt
    echo -e 1e4xScore/Cov"\t"$scorecov >> $1/allstats.txt
    echo -e Score/Num"\t"$scorenum >> $1/allstats.txt
    echo -e Span/Asm"\t"$spanasm >> $1/allstats.txt
    echo -e Cov/Asm"\t"$covasm >> $1/allstats.txt
    echo -e Cov/Span"\t"$covspan >>$1/allstats.txt
    echo -e Cov/Num"\t"$covnum >> $1/allstats.txt
}

function getmaligner {
    D=${BIONANO}/${1}/merge
    if [ ! -f ${D}/allstats.txt ]; then 
        calcmalignerstats  ${D}; 
    fi
    awk '{print $2}' ${D}/allstats.txt
}

function getsniffles {
    D=${LONG}/${1}/snifflestats
    PRE=${D}/ont
    ORDEREDONT="${PRE}.numaln  ${PRE}.numentries ${PRE}.numuniqaln ${PRE}.numuniqentries ${PRE}.sum.mapq ${PRE}.pctaln ${PRE}.alnratio ${PRE}.avg.mapq"
    PRE=${D}/pacbio
    ORDEREDPB="${PRE}.numaln  ${PRE}.numentries ${PRE}.numuniqaln ${PRE}.numuniqentries ${PRE}.sum.mapq ${PRE}.pctaln ${PRE}.alnratio ${PRE}.avg.mapq"
    ORDEREDSV="${D}/ontnumsv ${D}/pbnumsv ${D}/combnumsv ${D}/ontsumsv ${D}/pbsumsv ${D}/combsumsv"
    ORDEREDFILES="$ORDEREDONT $ORDEREDPB $ORDEREDSV"

    for f in ${ORDEREDFILES}; do
      b=`basename $f`
      ##echo -e $b"\t"`cat $f`
      echo `cat $f`
    done
}

function getsniffles_simple {
    D=${LONG}/${1}/snifflestats
    PRE=${D}/ont
    ## numalign and numuniqalign do not give extra info over pctaln and ratio
    ## numentries and numuniqentries are not useful here necessarily
    ## summapq and avgmapq are obv correlated -- but sum rewards more alignments given same avg... so I guess keep both... having said that, it is then also correlated with pctaln and maybe I should just use avg....
    ORDEREDONT="${PRE}.sum.mapq ${PRE}.pctaln ${PRE}.alnratio ${PRE}.avg.mapq"
    PRE=${D}/pacbio
    ORDEREDPB="${PRE}.sum.mapq ${PRE}.pctaln ${PRE}.alnratio ${PRE}.avg.mapq"
    ORDEREDSV="${D}/ontnumsv ${D}/pbnumsv ${D}/combnumsv ${D}/ontsumsv ${D}/pbsumsv ${D}/combsumsv"
    ORDEREDFILES="$ORDEREDONT $ORDEREDPB $ORDEREDSV"

    for f in ${ORDEREDFILES}; do
      b=`basename $f`
      ##echo -e $b"\t"`cat $f`
      echo `cat $f`
    done
}


function getall {
    line=$1
    b=$2
    getsizestats $line
    getbusco $b
    getbowtie2 $b  
    getale $b
    getlap $b
    getreapr $b
    getfrc $b
    getmaligner $b
    getsniffles $b
}

function getall_vizmat {
    line=$1
    b=$2
    getsizestats_vizmat $line
    getbusco $b
    getbowtie2_vizmat $b  
    getale_vizmat $b
    getlap $b
    getreapr_vizmat $b
    getfrc_vizmat $b
    getmaligner $b
    getsniffles $b
}

function getsimple {
    line=$1
    b=$2
    getsizestats_vizmat $line
    getbusco_simple $b
    getbowtie2_simple $b  
    getale_simple $b
    getlap $b
    getreapr_simple $b
    getfrc_simple $b
    getmaligner $b
    getsniffles_simple $b
}


function debugmode {
    line=$1
    b=$2
    SIZE=`getsizestats $line | wc -l`
    BUSCO=`getbusco $b | wc -l`
    BT2=`getbowtie2 $b  | wc -l`
    ALE=`getale $b | wc -l`
    LAP=`getlap $b | wc -l`
    REAP=`getreapr $b | wc -l`
    FRC=`getfrc $b | wc -l`
    BNG=`getmaligner $b | wc -l`
    SNIF=`getsniffles $b | wc -l`
    for var in SIZE BUSCO BT2 ALE LAP REAP FRC BNG SNIF; do
        echo -e $var"\t"${!var}
    done
}

function getnames {
    echo $NAMES
}

function getnames_vizmat {
    echo $VIZMATNAMES 
}

function getnames_simple {
    echo $SIMPLENAMES
}

function getnames_simple {
    echo $SIMPLENAMES
}

function getfactors_simple {
    echo $SIMPLEFACTORS
}

function main {
    LONGTABLE=longtables
    if [ ! -d $LONGTABLE ]; then mkdir $LONGTABLE; fi
    if [ $NARG -eq 2 ] && [ $DEBUG == "debug" ]; then
        ##echo DEBUG
        while read line; do
            b=`basename $line .fasta`
            debugmode $line $b > $LONGTABLE/$b.longtable.debugmode
        done < $FOFN
    elif [ $NARG -eq 2 ] && [ $DEBUG == "vizmat" ]; then
        ##echo VIZMAT
        while read line; do
            b=`basename $line .fasta`
            getall_vizmat $line $b > $LONGTABLE/$b.longtable.vizmat
        done < $FOFN
    elif [ $NARG -eq 2 ] && [ $DEBUG == "simple" ]; then
        ##echo SIMPLE
        while read line; do
            b=`basename $line .fasta`
            getsimple $line $b > $LONGTABLE/$b.longtable.simple
        done < $FOFN
    elif [ $NARG -eq 2 ] && [ $DEBUG == "all" ]; then
        ##echo ALL
        while read LINE; do
            B=`basename $LINE .fasta`
            getall $LINE $B > $LONGTABLE/$b.longtable
            getall_vizmat $LINE $B > $LONGTABLE/$b.longtable.vizmat
            getsimple $LINE $B > $LONGTABLE/$b.longtable.simple
            debugmode $LINE $B > $LONGTABLE/$b.longtable.debugmode
        done < $FOFN
    elif [ $NARG -eq 2 ] && [ $DEBUG == "names" ]; then
        ##echo NAMES
        getnames
    elif [ $NARG -eq 2 ] && [ $DEBUG == "vizmatnames" ]; then
        ##echo VIZMAT NAMES
        getnames_vizmat
    elif [ $NARG -eq 2 ] && [ $DEBUG == "simplenames" ]; then
        ##echo SIMPLE NAMES
        getnames_simple
    elif [ $NARG -eq 2 ] && [ $DEBUG == "simplefactors" ]; then
        ##echo SIMPLE NAMES
        getfactors_simple
    else
        ##echo LONGTABLE (LaTeX) 
        while read line; do
            b=`basename $line .fasta`
            getall $line $b > $LONGTABLE/$b.longtable
        done < $FOFN
    fi
}


### EXECUTE
main
