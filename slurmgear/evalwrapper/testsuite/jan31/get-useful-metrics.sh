#!/bin/bash

## John Urban, 2016- 2017; v1
## For collecting metrics out of Sciara assembly directories...
## This code is very much hacked together and will need to be revised for more flexibility...

function help {
    echo "Usage: $0 -i:g:o:dfncavht
        -i with argument = path to FOFN (defaults to input.fofn)
        ..where FOFN has list of all assemblies used in the assembly evaluations in subdirs you are trying to summarize.
          (( typically called input.fofn ))
        -g with argument = expected genome size (default = 292000000)
        -o with argument = outdirectory name (default = tables)
        -d debug mode: create files that count the number of metrics from each function in 2-col tab-delim file
        -f get factors -- get the multiplication vector to use against each pair of metric scores before ranking (in R)
        -n get names -- returns names of all metrics...
        -c get categories -- returns categories of all metrics...
        -a get all -- populates out-directory with the score.txt, all.scores.tab, debug.tsv, names.text, and factors.text
        -v verbose -- does nothing at the moment....
        -h help - returns this message; also returns this when no arguments given
        -t test - just a useful option for testing things during development

"
}


if [ $# -eq 0 ]; then help; exit; fi

##Defaults:
FOFN=input.fofn
NARG=$#
DEBUG=false
EXPGENSIZE=292000000
FACTORS=false
NAMES=false
CATEGORIES=false
ALL=false
VERBOSE=false
HELP=false
OUTDIR="tables"
TEST=false

## ASSUMES FOLLOWING DIRS -- TODO -- add options to change these
SHORT=shortread
BIONANO=bionano
LONG=longread
BLAST=blast_analyses
TRANSCRIPTOME=${BLAST}/transcriptome
KNOWN=${BLAST}/knownseqs
DMEL=${BLAST}/dmel
DMELPEP=${BLAST}/dmel_peptides
MOSQ=${BLAST}/anopheles
MOSQPEP=${BLAST}/anopheles_peptides
RNA=rnaseq
BUSCOV3=buscov3

while getopts "i:g:o:dfncavht" arg; do
    case $arg in
        i) FOFN=$OPTARG;;
        g) EXPGENSIZE=$OPTARG;;
        o) OUTDIR=$OPTARG;;
        d) DEBUG=true;;
        f) FACTORS=true;;
        n) NAMES=true;;
        c) CATEGORIES=true;;
        a) ALL=true;;
        v) VERBOSE=true;;
        h) HELP=true;;
        t) TEST=true;;
        *) help;;
    esac
done

## HELP
if ${HELP}; then help; exit; fi




##NAME VARIABLES
SIMPLENAMES="Num_Contigs,Max_Contig_Size,NG50,LG50,E_G,Complete_BUSCOs,Bowtie2_Ilmn_Pct,Bowtie2_Pct_Aligned_Concordantly,Bowtie2_Pct_Aligned_Disconcordantly,ALE_Score,ALE_Place_Avg,ALE_Insert_Avg,ALE_kmer_Avg,ALE_Depth_Score_Avg,LAP,REAPR_MBS,REAPR_Pct_EF,REAPR_FCD_Errors_within_a_contig,REAPR_FCD_Errors_over_a_gap,REAPR_Low_frag_cov_within_a_contig,REAPR_Low_frag_cov_over_a_gap,REAPR_Low_Score_Regions,REAPR_Links,REAPR_Soft_Clip,REAPR_Collapsed_repeats,REAPR_Low_Read_Coverage,REAPR_Low_Perfect_Coverage,REAPR_Wrong_Orientation,REAPR_Broken_Num_Contigs,REAPR_Broken_Max_Contig_Length,REAPR_Broken_NG50,REAPR_Broken_LG50,REAPR_Broken_EG,REAPR_Broken_Num_Gaps,REAPR_Broken_Total_Gap_Len,FRC_Num_Features,FRC_Norm_Num_Features_(per_Mb),FRC_Proper,FRC_Wrong_Dist,FRC_Wrong_Orientation,FRC_Wrong_Contig,FRC_Singleton,FRC_Mean_Cov,FRC_Spanning_Cov,FRC_Proper_Pairs_Cov,FRC_Wrong_Mate_Cov,FRC_Singleton_Mate_Cov,FRC_Different_Contig_Cov,BNG_Score,BNG_Span,BNG_Cov,BNG_Num,BNG_1e4*Score/Cov,BNG_Score/Num,BNG_Span/AsmSize,BNG_Cov/AsmSize,BNG_Cov/Span,BNG_Cov/Num,ONT_sum_mapq,ONT_pct_aln,ONT_ratio,ONT_avg_mapq,PacBio_sum_mapq,PacBio_Pct_aln,PacBio_ratio,PacBio_avg_mapq,Sniffles_ONT_Num_SV,Sniffles_PacBio_Num_SV,Sniffles_ONT_Sum_SV,Sniffles_PacBio_Sum_SV"

## R ranks smallest to biggest as 1 to N. When treating 1 as best rank, this means that things that are better as they have higher/bigger values need to be multiplied by -1, which reverse the direction of ranking... e.g. 1,5,10 becomes -10,-5,-1.
SIMPLEFACTORS="1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,1,1,1"

## ECHO TO STDERR : "ERRCHO"
function errcho {
 echo $@ 1>&2
}

##VERBOSITY
function verbose { 
 msg=${1}
 if $VERBOSE; then errcho ${msg}; fi
}

function getsizestats { ##takes fasta
    ## 6 outputs: ncontigs, asmlen, maxlen, ng50, lg50, eg
    b=`basename $1 .fasta`
    F=sizestats/${b}.tsv
    if [ ! -d sizestats ]; then mkdir sizestats; fi
    if [ ! -f $F ]; then faSize -detailed $1 | asm-stats.py -t -G $G | awk '{gsub(/,/,"\n"); print}' > $F ; fi
    awk 'NR==1 || NR==2 || NR==3 || NR==10 || NR==11 || NR==12 {print $1}' $F
}



function getsizestats_vizmat { ##takes fasta
    ## 5 outputs: ncontigs, maxlen, ng50, lg50, eg
    ## factors: 1,-1,-1,1,-1
    b=`basename $1 .fasta`
    F=sizestats/${b}.tsv
    if [ ! -d sizestats ]; then mkdir sizestats; fi
    if [ ! -f $F ]; then faSize -detailed $1 | asm-stats.py -t | awk '{gsub(/,/,"\n"); print}' > $F ; fi
    ## asmsize taken out b/c not meaningful in terms of an asm getting better or worse (Except in terms of going further away from expected size)
    awk 'NR==1 || NR==3 || NR==10 || NR==11 || NR==12 {print $1}' $F
}

function getbusco_simple {
    ## factors: -1
    ## 1 output: complete
    ##F=$SHORT/${1}/busco/*/short*
    F=$1
    ## Removed
    ## "Complete and single-copy BUSCOs" "Complete and duplicated BUSCOs" "Fragmented BUSCOs" "Missing BUSCOs"
    ## Since it is unclear w/o specific knowledge that having more or less of the complete buscos as single copy vs duplicated is "better" for an asm
    ## and since missing buscos gives no extra info over complete
    for l in "Complete BUSCOs"; do
        grep "${l}" $F | awk '{print $1}'
    done
}

function getbuscov3_simple {
    F=$1
    complete=`grep "Complete BUSCOs" $F | awk '{print $1}'`
    fragment=`grep "Fragmented BUSCOs" $F | awk '{print $1}'`    
    echo $complete
    echo $complete+$fragment | bc
    ## factors: -1,-1
}

function getbowtie2_simple {
    ## 2 outputs: overall aln rate, %of all pairs conc
    ## no longer do this as it is hard to use for directo comparison: % of all pairs disc 
    ## --- e.g. if 50% of pairs align with 40% concordant and 10% discordant vs 100% pairs with 80% concordant and 20% discordant -- the first asm would be rewarded for having fewer discordant... though it only has fewer b/c fewer aligned...
    ##F=$SHORT/${1}/mreads/*.err
    F=$1
    grep "overall alignment rate" $F | awk '{sub(/%/,""); print $1}'
    pair=`grep "were paired; of these:" $F | awk '{print $1}'`
    conc1=`grep "aligned concordantly exactly 1 time" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'`
    concmult=`grep "aligned concordantly >1 times" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'`
    ## it is not immediately clear if more or less single vs multi is better -- easier to interpret %concordant marginalized over both
    echo $conc1 $concmult | awk '{print $1+$2}'   ## pct aligned concordantly
    ## removed pct DNAC since that does not give any new info wrt pct conc
    ## remove "pct of reads aln conc 0 times that aln disc" since this is not immediately interp as better/worse
    #disc1=`grep "aligned discordantly 1 time" $F | awk '{print $1}'`
    #pctdisc1=`echo $disc1 $pair | awk '{print 100.0*$1/$2}'`
    #echo ${pctdisc1} #pct of total pairs that aln disc
    ## pct unmap does not nec give any more info 
    ## SE read are not immediately comparable/interpretable -- I will let ALE/LAP/REAPR/FRC deal with them
    ##factors: -1,-1,
}



function getale_simple {
    F=$1
    #F=$SHORT/${1}/ale/*ALE.txt
    ## 1=ALE score, 2=ncont, 3=asmlen, 4=placeAvg, 5=insertAvg, 6=kmerAvg, 7=depthScoreAvg, 8=depthAvg, 9=totalreads, 10=totalMappedReads, 11=totalUnMappedReads, 12=total placed reads, 13=readlenAvg, 14=avgReadOvlap
    ## trimmed off mapped, unmapped, and placed -- will let those be reflected in BT2 fxn and as part of the other scores here.
    awk 'NR==1 || NR==4 || NR==5 || NR==6 || NR==7 {print $3}' $F
    ## factors: -1,-1,-1,-1,-1
}

function getale_simple_long {
    F=$1
    #F=$SHORT/${1}/ale/*ALE.txt
    ## 1=ALE score, 2=ncont, 3=asmlen, 4=placeAvg, 5=insertAvg, 6=kmerAvg, 7=depthScoreAvg, 8=depthAvg, 9=totalreads, 10=totalMappedReads, 11=totalUnMappedReads, 12=total placed reads, 13=readlenAvg, 14=avgReadOvlap
    ## trimmed off mapped, unmapped, and placed -- will let those be reflected in BT2 fxn and as part of the other scores here.
    awk 'NR==1 || NR==4 || NR==6 || NR==7 {print $3}' $F
}

function getlap {
    #F=$SHORT/${1}/lap/*.lapscore
    F=$1
    awk '{print $1}' $F | head -n 1
}


function getasmstats {
    ## WARN: assumes getsizestats was already run
    F=sizestats/${b}.tsv
    paste -sd"\t" $F
}

function getreapr_simple {
    ## FOR REAPR SIMPLE -- rm sumError and sumWarn -- just use indiv err and warns
    ##D=$SHORT/${1}/reapr/output_directory/
    D=$1
    ###GGG=$2
    head -n 1 $D/per-base-mean-score.txt
    ## pctEF, FCD, FCD_gap, frag_cov, frag_cov_gap, lowscore, link, softclip, collapsed repeat, readcov, low perfect cov, readorientation
    ##awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'
    ## get rid of FCD/cov over gap metrics
    ## pctEF, FCD, frag_cov, lowscore, link, softclip, collapsed repeat, readcov, low perfect cov, readorientation
    awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39, $41, $43, $44, $45, $46, $47, $48, $49}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'
    
    ## Broken Asm -size stats -- removed for now b/c it really only seems to matter when aggressive breaking is used since most long-read contigs do not have Ns...
    ## br nseqs,  br longest, 
    #awk 'NR==2 {OFS="\t"; print $21, $23}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
    #stats=`getasmstats $1`
    #bstats=`awk '{gsub(/,/,"\t"); print }' ${D}/broken_assembly.sizestats.csv`
    #echo $bstats | awk '{print $10}' ##broken ng50
    #echo $bstats | awk '{print $11}' ##broken lg50
    #echo $bstats | awk '{print $12}' ## broken eg

    ## Below is only useful when not aggressively breaking the assembly as it basically tells you how much it would break it up.
    ## After it is aggressively broken though, these are usually zeroed out -- and therefore not useful to rank. All in all, either use (i) broken asm size stats OR (ii) reapr's gap stats
    ## ngaps, gaplen
    awk 'NR==2 {OFS="\t"; print $36, $37}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    

    # PAST OUTPUT = #1=mbs, 2=pctEF, 3=FCD, 4=frag_cov, 5=lowscore, 6=link, 7=softclip, 8=collapsed repeat, 9=readcov, 10=low perfect cov, 11=readorientation, 12=br nseqs,  13=broken longest, 14=broken ng50, 15=broken lg50, 16=broken eg, 17=ngaps, 18=gaplen
    # NEW OUTPUT = #1=mbs, 2=pctEF, 3=FCD, 4=frag_cov, 5=lowscore, 6=link, 7=softclip, 8=collapsed repeat, 9=readcov, 10=low perfect cov, 11=readorientation, 12=ngaps, 13=gaplen
}

function getreapr_simple_longread {
    D=$1
    head -n 1 $D/per-base-mean-score.txt
    ## From next awk command: removed $48 -> formerly output #10 --> low perfect cov --> b/c perfect mapping not used for long2pe
    awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39, $41, $43, $44, $45, $46, $47, $49}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'
    awk 'NR==2 {OFS="\t"; print $36, $37}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
    # OUTPUT = #1=mbs, 2=pctEF, 3=FCD, 4=frag_cov, 5=lowscore, 6=link, 7=softclip, 8=collapsed repeat, 9=readcov, 10=readorientation, 11=ngaps, 12=gaplen
}

function getreapr_simple_longread_lowcov {
    D=$1
    head -n 1 $D/per-base-mean-score.txt
    ## As with PB version -- From next awk command: removed $48 -> formerly output #10 --> low perfect cov --> b/c perfect mapping not used for long2pe
    ## In addition, From next awk command: removed $41 --> formerly output #4 --> frag_cov (frag cov too low) --> b/c ONT coverage too low for default minimum of 1 where there were too many false positives in ecoli analysis -- ONT min frag cov was set to 0 here to ignore this...
    awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39, $43, $44, $45, $46, $47, $49}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'
    awk 'NR==2 {OFS="\t"; print $36, $37}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
    # OUTPUT = #1=mbs, 2=pctEF, 3=FCD, 4=lowscore, 5=link, 6=softclip, 7=collapsed repeat, 8=readcov, 9=readorientation, 10=ngaps, 11=gaplen
}



function getfrc_simple {
    ## WARN: assumes getasmstats was already run
    ## 6=MAPPED
    ## numfrc, normnumfrc, 8=PROPER,9=WRONG_DIST, 11=WRONG_ORIENTATION,12=WRONG_CONTIG,13=SINGLETON, 14=MEAN_COVERAGE,15=SPANNING_COVERAGE,16=PROPER_PAIRS_COVERAGE,17=WRONG_MATE_COVERAGE,18=SINGLETON_MATE_COV,19=DIFFERENT_CONTIG_COV
    ## For now, it doesnt give the raw counts of reads that mapped in a proper pair, w/ wrong dit, w/ wrong orientation, ...etc
    ##   instead it gives all those normalized to number of mapped reads.
    ## The problem with using raw counts is best illustrated with an example:
    ##  If a bad asm has very few reads mapping and they all map with wrong dist or wrong orientn, those counts will still be lower than a good asm with tons more reads mapping and a lower pct mapping wrong...
    ##  When normalized it becomes obvious how to interp 100% wrong mapped vs 10%
    ##  Another example is an assembly with
    F1=$1
    F2=$2
    numfrc=`grep -c -v ^# $F1`
    asmlen=`getasmstats | awk '{print $2}'`
    normfrc=`echo $numfrc $asmlen | awk '{print 1e6*$1/$2}'`
    echo $numfrc
    echo $normfrc
    ## pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wrong contig, pct of mapped as singleton....
    tail -n 1 $F2 | awk 'OFS="\t" {gsub(/,/,"\t"); print 100.0*$8/$6, 100.0*$9/$6, 100.0*$11/$6, 100.0*$12/$6, 100.0*$13/$6}' | awk '{gsub(/\t/,"\n"); print}'
}


function calcmalignerstats {
    #Takes $D for merge directory as $1
    if [ ! -f $1/score.txt ]; then tail -n +2 $1/all.bionano.smoothed.maps.aln | awk '{s+=$19}END{print s}' > $1/score.txt; fi
    if [ ! -f $1/span.txt ]; then awk '{s+=($3-$2)}END{print s}' $1/all.bionano.smoothed.maps.aln.bedGraph > $1/span.txt; fi
    if [ ! -f $1/total_base_cov.txt ]; then awk '{s+=($3-$2)*$4}END{print s}' $1/all.bionano.smoothed.maps.aln.bedGraph > $1/total_base_cov.txt; fi
    if [ ! -f $1/num_alignments.txt ]; then tail -n +2 $1/all.bionano.smoothed.maps.aln | wc -l > $1/num_alignments.txt; fi
    GENOMEFILE=$1/*.genome
    asmsize=`awk '{s+=$2}END{print s}' $GENOMEFILE`

    ## calculate metrics
    score=`cat $1/score.txt`
    span=`cat $1/span.txt`
    cov=`cat $1/total_base_cov.txt`
    num=`cat $1/num_alignments.txt`
    scorecov=`python -c "print 1e4*$score/$cov.0"`
    scorenum=`python -c "print $score/$num.0"`
    asmsize=`awk '{s+=$2}END{print s}' $GENOMEFILE`
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
    grep score ${D}/allstats.txt | awk '{print $2}'
    grep span ${D}/allstats.txt | awk '{print $2}'
    grep total_cov ${D}/allstats.txt | awk '{print $2}'
    grep num_aln ${D}/allstats.txt | awk '{print $2}'
    grep "Score/Num" ${D}/allstats.txt | awk '{print $2}'
    grep "Span/Asm" ${D}/allstats.txt | awk '{print $2}'
    grep "Cov/Asm" ${D}/allstats.txt | awk '{print $2}'
    grep "Cov/Span" ${D}/allstats.txt | awk '{print $2}'
    grep "Cov/Num" ${D}/allstats.txt | awk '{print $2}'

    # OUTPUT:
    # Formerly: score, span, total_cov, num, 1e4*Score/Cov, avgscore=Score/Num, pctAsmSpanned=Span/Asm, AvgCovPerBaseInAsm=Cov/Asm, AvgCovPerSpannedBase=Cov/Span, AvgMolSize=Cov/NumAln
    # Now: 1e4*Score/Cov was removed b/c it did not seem to be able to rank assemblies well in E.coli testing -- the others did better
    ##Now:  score, span, total_cov, num, avgscore=Score/Num, pctAsmSpanned=Span/Asm, AvgCovPerBaseInAsm=Cov/Asm, AvgCovPerSpannedBase=Cov/Span, AvgMolSize=Cov/NumAln
}


function getlongreadstats_simple {
    for VAR in pctalnlen pctmatch_allreads pctmatch_alnreads pctmatch_alns pct_reads_that_aln alnratio pct_multi_ctg avgalnlen_allaln avgalnlen_allalnread avgalnlen_allread avgmapq_per_aln avgmapq_per_read avgmapq_with_unmap avgalnscore_per_aln avgalnscore_aln_reads avgalnscore_allreads; do
      grep -w $VAR ${1} | awk '{print $2}'
    done
}

function getsniffles_simple {
    #D=${LONG}/${1}/stats
    # give path to sniffles_x dir
    D=${1}
    DEL=${D}/numdel
    DUP=${D}/numdup
    INS=${D}/numins
    INV=${D}/numinv
    TRA=${D}/numtra
    NUM=${D}/*numsv.tmp
    cat $DEL $DUP $INS $INV $TRA | awk -v "s=0" '{s+=$1}END{print s}' > $NUM
    SUM=${D}/*sumsv
    ORDEREDFILES="${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}"

    for f in ${ORDEREDFILES}; do
      ##bn=`basename $f`
      ##echo -e $bn"\t"`cat $f`
      echo `cat $f`
    done
}




function getpilon {
 PILONDIR=${1}
 #G=${2}

 #get asm size
 ASMSIZE=`getasmstats | awk '{print $2}'`

 # pct confirmed - good to see proportion of assembly in good standing, but can penalize assemblies that indeed have more confirmed bases at the expense of a larger asm size that makes pct go down
 awk 'NR==2 {print $3}' ${PILONDIR}/confirmed.txt

 # pct of expected genome size confirmed -- normalizes all to expected genome size to see pct of G confirmed
 awk -v "G=${EXPGENSIZE}" 'NR==2 {print 100.0*$1/G}' ${PILONDIR}/confirmed.txt
 
 # number SNPs found per Mb
 awk 'NR==2 {print $1}' ${PILONDIR}/vars.txt

 # number small INSs found per Mb
 awk 'NR==2 {print $2}' ${PILONDIR}/vars.txt

 # number small DELs found per Mb
 awk 'NR==2 {print $4}' ${PILONDIR}/vars.txt

 # number SNPs+INSs+DELs found per Mb
 awk 'NR==2 {print $1+$2+$4}' ${PILONDIR}/vars.txt

 # total num alt alleles in VCF per Mb
 totalnumaltalleles=`cat ${PILONDIR}/num_alt_alleles.txt`
 echo 1000000.0*${totalnumaltalleles}/${ASMSIZE} | bc -l

 # Pct of asm flagged as Large collapsed region -- Note this can be inflated for assemblies with bubbles kept if Pilon is identifying 10kb+ stretches with >= 2*Global_Mean_cov -- perhaps not if it is local cov on contig....
 LCRLenSum=`grep "^Large collapsed region" ${PILONDIR}/pilon.err | awk -v "s=0" '{s+=$6}END{print s}'`
 echo 100.0*${LCRLenSum}/${ASMSIZE} | bc -l
 
 # Num LCRs per Mb
 NumLCR=`grep -c "^Large collapsed region" ${PILONDIR}/pilon.err`
 echo 1000000.0*${NumLCR}/${ASMSIZE} | bc -l


}


function getblastanalysis {
  #identity_len,    sum_bit,   len,     num_queries_ge_1_hit
  awk 'NR==2 {OFS="\n"; print $1,$2,$3,$4}' ${1}
}


function getrnaseqinfo {
  for var in pct_pairs_not_same_ctg pct_pairs_with_1-2mappedmates_not_same_ctg pct_pairs_bothmapped_not_same_ctg pct_pairs_not_same_ctg_qge2 pct_bothmapped_not_same_ctg_qge2 pct_bothmappeddiff_with_qge2 avg_min_mapq_all_pairs avg_min_mapq_pairs_bothmatesmapped overall_alignment_rate pct_of_mapped_pairs_that_are_concordant ; do
    grep -w $var ${1} | awk '{print $2}'
  done
}


function getsimple {
     ## 201 outputs total
     verbose ${b}

     ##SizeL (5 outputs) ncontigs, maxlen, ng50, lg50, eg
     verbose "size..."
     getsizestats_vizmat $line

     ##short read: (37 outputs total)
     ## bt2: (2 outputs) overall aln rate, %of all pairs conc, % of all pairs disc
     verbose "short bt2..."
     getbowtie2_simple $SHORT/${b}/mreads/*.err 

     ## ALE: (5 outputs) ALE score, ALE placeAvg, ALE insertAvg, ALE kmerAvg, ALE depthScoreAvg
     verbose "short ale..."
     getale_simple $SHORT/${b}/ale/*ALE.txt 

     ## LAP: (1 output) lap score
     verbose "short lap..."
     getlap $SHORT/${b}/lap/*.lapscore

     ## REAPR: (13 outputs) MBS, pctEF, FCD, frag_cov, lowscore, link, softclip, collapsed repeat, readcov, low perfect cov, readorientation, ngaps, gaplen
     verbose "short reapr..."
     getreapr_simple $SHORT/${b}/reapr/output_directory/ 

     ## FRC: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wrong contig, pct of mapped as singleton
     verbose "short frc..."
     getfrc_simple $SHORT/${b}/frc/*gff $SHORT/${b}/frc/*frc_assemblyTable.csv

     ## PILON: (9 outputs) pct confirmed, pct Gsize confirmed, nSNPs/Mb, nSmallINS/Mb, nSmallDEL/Mb, n SNPs+INSs+DELs found per Mb, total num alt alleles in VCF per Mb, Pct of asm flagged as Large collapsed region, Num LCRs per Mb
     verbose "short pilon..."
     getpilon ${SHORT}/${b}/pilonvars/ ${EXPGENSIZE}

     #PACBIO (46 outputs total)
     #PacBio - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns, pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     verbose "pacbio aln..."
     getlongreadstats_simple ${LONG}/${b}/mreads/pacbio*simple.txt

     ## SNIFFLES: (7 outputs) ${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}
     verbose "pacbio sniffles..."
     getsniffles_simple ${LONG}/${b}/sniffles_pb

     ## ALE LONG: (4 outputs) 1=ALE score, 4=placeAvg, 6=kmerAvg, 7=depthScoreAvg
     verbose "pacbio ale..."
     getale_simple_long $LONG/${b}/ale_pb/*ALE.txt

     ## FRC LONG: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wron...
     verbose "pacbio frc..."
     getfrc_simple $LONG/${b}/frc_pb/*gff $LONG/${b}/frc_pb/*_assemblyTable.csv

     ## REAPR LONG HIGH COV: (12 outputs) 1=mbs, 2=pctEF, 3=FCD, 4=frag_cov, 5=lowscore, 6=link, 7=softclip, 8=collapsed repeat, 9=readcov, 10=readorientation, 11=ngaps, 12=gaplen
     verbose "pacbio reapr..."
     getreapr_simple_longread $LONG/${b}/reapr_pb/output_directory/ 

     # NANOPORE (45 outputs total)
     #ONT - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) - same as above
     verbose "ont aln..."
     getlongreadstats_simple ${LONG}/${b}/mreads/ont*simple.txt

     ##	SNIFFLES: (7 outputs) -- same as above
     verbose "ont sniffles..."
     getsniffles_simple ${LONG}/${b}/sniffles_ont

     ## ALE LONG: (4 outputs) -- same as above
     verbose "ont ale..."
     getale_simple_long $LONG/${b}/ale_ont/*ALE.txt

     ## FRC LONG: (7 outputs) -- same as above
     verbose "ont frc..."
     getfrc_simple $LONG/${b}/frc_ont/*gff $LONG/${b}/frc_ont/*_assemblyTable.csv

     ##	REAPR LONG LOW COV: (11 outputs - fragcov removed) 1=mbs, 2=pctEF, 3=FCD, 4=lowscore, 5=link, 6=softclip, 7=collapsed repeat, 8=readcov, 9=readorientation, 10=ngaps, 11=gaplen
     verbose "ont reapr..."
     getreapr_simple_longread_lowcov $LONG/${b}/reapr_ont/output_directory/ 

     #BioNano (9 outputs) score, span, total_cov, num, avgscore=Score/Num, pctAsmSpanned=Span/Asm, AvgCovPerBaseInAsm=Cov/Asm, AvgCovPerSpannedBase=Cov/Span, AvgMolSize=Cov/NumAln
     verbose "bionano..."
     getmaligner $b

     # Old Busco: (1 outout) number complete
     verbose "old busco...."
     getbusco_simple $SHORT/${b}/busco/*/short*

     # New Busco - (6 clades, 2 outputs each, 12 outputs total) look at both Num Complete  and   NumComplete+NumFrags
     verbose "newbusco..."
     for lineage in eukaryota metazoa arthropoda insecta endopterygota diptera; do
       getbuscov3_simple ${BUSCOV3}/${b}/${lineage}/*/short*     
     done

     # Known Seq: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     verbose "known...."
     getblastanalysis ${KNOWN}/${b}/analysis/total.txt

     # Transcriptome (8 outputs total)
     verbose "transcriptome..."
     ## BLASTN: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${TRANSCRIPTOME}/${b}/analysis/total.txt
     ##TBLASTX: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${TRANSCRIPTOME}/${b}/analysis-tbx/total.txt

     # Dmel (12 outputs total)
     verbose "dmel...."
     ##	BLASTN:	(4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${DMEL}/${b}/analysis/total.txt
     ##TBLASTX:	(4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${DMEL}/${b}/analysis-tbx/total.txt
     ##TBLASTN:	(4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${DMELPEP}/${b}/analysis/total.txt

     # Anopheles (12 outputs total)
     verbose "anopheles..."
     ## BLASTN: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${MOSQ}/${b}/analysis/total.txt
     ##TBLASTX:        (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${MOSQ}/${b}/analysis-tbx/total.txt
     ##TBLASTN:        (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     getblastanalysis ${MOSQPEP}/${b}/analysis/total.txt

     # RNA-seq (10 outputs) pct_pairs_not_same_ctg, pct_pairs_with_1-2mappedmates_not_same_ctg, pct_pairs_bothmapped_not_same_ctg, pct_pairs_not_same_ctg_qge2, pct_bothmapped_not_same_ctg_qge2, pct_bothmappeddiff_with_qge2, avg_min_mapq_all_pairs, avg_min_mapq_pairs_bothmatesmapped, overall_alignment_rate, pct_of_mapped_pairs_that_are_concordant
     verbose "rna-seq..."
     getrnaseqinfo ${RNA}/${b}/mreads/info.txt

     #END
     verbose " "
}


########################################################################## FACTORS ##################################################################################

function getfactors {
     ## 201 outputs total
     FAC=""
     ##SizeL (5 outputs) ncontigs, maxlen, ng50, lg50, eg
     FAC+="1,-1,-1,1,-1,"

     ##short read: (37 outputs total)
     ## bt2: (2 outputs) overall aln rate, %of all pairs conc
     FAC+="-1,-1,"

     ## ALE: (5 outputs) ALE score, ALE placeAvg, ALE insertAvg, ALE kmerAvg, ALE depthScoreAvg
     FAC+="-1,-1,-1,-1,-1,"

     ## LAP: (1 output) lap score
     FAC+="-1,"

     ## REAPR: (13 outputs) MBS, pctEF, | FCD, frag_cov, lowscore, link, | softclip, collapsed_repeat, readcov,  low perfect cov, | readorientation, ngaps, gaplen
     FAC+="-1,-1,"
     FAC+="1,1,1,1,"
     FAC+="1,1,1,1,"
     FAC+="1,1,1,"

     ## FRC: (7 outputs) num feat, norm num feat | pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wrong contig, pct of mapped as singleton
     FAC+="1,1,"
     FAC+="-1,1,1,1,1,"
     
     ## PILON: (9 outputs) pct confirmed, pct Gsize confirmed, nSNPs/Mb, nSmallINS/Mb, nSmallDEL/Mb, n SNPs+INSs+DELs found per Mb, total num alt alleles in VCF per Mb, Pct of asm flagged as Large collapsed region, Num LCRs per Mb
     FAC+="-1,-1,1,1,1,1,1,1,1,"

     #PACBIO (46 outputs total)
     #PacBio - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns,
     FAC+="-1,-1,-1,-1,"
     ## pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, 
     FAC+="-1,1,1,-1,"
     ## avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, 
     FAC+="-1,-1,-1,-1,"
     ## avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     FAC+="-1,-1,-1,-1,"

     ## SNIFFLES: (7 outputs) ${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}
     FAC+="1,1,1,1,1,1,1,"

     ## ALE LONG: (4 outputs) 1=ALE score, 4=placeAvg, 6=kmerAvg, 7=depthScoreAvg
     FAC+="-1,-1,-1,-1,"

     ## FRC LONG: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wron...
     FAC+="1,1,"
     FAC+="-1,1,1,1,1,"

     ## REAPR LONG HIGH COV: (12 outputs) 1=mbs, 2=pctEF, | 3=FCD, 4=frag_cov, 5=lowscore, 6=link, | 7=softclip, 8=collapsed repeat, 9=readcov, | 10=readorientation, 11=ngaps, 12=gaplen
     FAC+="-1,-1,"
     FAC+="1,1,1,1,"
     FAC+="1,1,1,"
     FAC+="1,1,1,"

     # NANOPORE (45 outputs total)
     #ONT - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) - same as above
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns,
     FAC+="-1,-1,-1,-1,"
     ## pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, 
     FAC+="-1,1,1,-1,"
     ## avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, 
     FAC+="-1,-1,-1,-1,"
     ## avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     FAC+="-1,-1,-1,-1,"

     ##	SNIFFLES: (7 outputs) -- same as above
     FAC+="1,1,1,1,1,1,1,"

     ## ALE LONG: (4 outputs) -- same as above
     FAC+="-1,-1,-1,-1,"

     ## FRC LONG: (7 outputs) -- same as above
     FAC+="1,1,"
     FAC+="-1,1,1,1,1,"

     ##	REAPR LONG LOW COV: (11 outputs - fragcov removed) 1=mbs, 2=pctEF, | 3=FCD, 4=lowscore, 5=link, | 6=softclip, 7=collapsed repeat, 8=readcov, | 9=readorientation, 10=ngaps, 11=gaplen
     FAC+="-1,-1,"
     FAC+="1,1,1,"
     FAC+="1,1,1,"
     FAC+="1,1,1,"

     #BioNano (9 outputs) score, span, total_cov, num, avgscore=Score/Num, pctAsmSpanned=Span/Asm, AvgCovPerBaseInAsm=Cov/Asm, AvgCovPerSpannedBase=Cov/Span, AvgMolSize=Cov/NumAln
     FAC+="1,-1,-1,-1,1,-1,-1,-1,-1,"

     # Old Busco: (1 outout) number complete
     FAC+="-1,"

     # New Busco - (6 clades, 2 outputs each, 12 outputs total) look at both Num Complete  and   NumComplete+NumFrags
     FAC+="-1,-1,"
     FAC+="-1,-1,"
     FAC+="-1,-1,"
     FAC+="-1,-1,"
     FAC+="-1,-1,"
     FAC+="-1,-1,"

     # Known Seq: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     FAC+="-1,-1,-1,-1,"

     # Transcriptome (8 outputs total)
     FAC+="-1,-1,-1,-1,"
     FAC+="-1,-1,-1,-1,"

     # Dmel (12 outputs total)
     FAC+="-1,-1,-1,-1,"
     FAC+="-1,-1,-1,-1,"
     FAC+="-1,-1,-1,-1,"

     # Anopheles (12 outputs total)
     FAC+="-1,-1,-1,-1,"
     FAC+="-1,-1,-1,-1,"
     FAC+="-1,-1,-1,-1,"

     # RNA-seq (10 outputs) pct_pairs_not_same_ctg, pct_pairs_with_1-2mappedmates_not_same_ctg, pct_pairs_bothmapped_not_same_ctg, pct_pairs_not_same_ctg_qge2, 
     FAC+="1,1,1,1,"
     ## pct_bothmapped_not_same_ctg_qge2, pct_bothmappeddiff_with_qge2, avg_min_mapq_all_pairs, avg_min_mapq_pairs_bothmatesmapped, 
     FAC+="1,1,-1,-1,"
     ## overall_alignment_rate, pct_of_mapped_pairs_that_are_concordant
     FAC+="-1,-1"

     #END
     echo $FAC | awk '{gsub(/,/,"\n"); print}'
}

########################################################################## NAMES ##################################################################################

function getnames {
     ## 201 outputs total
     FAC=""
     ##SizeL (5 outputs) ncontigs, maxlen, ng50, lg50, eg
     FAC+="nContigs,maxCtgLen,NG50,LG50,EG,"

     ##short read: (37 outputs total)
     ## bt2: (2 outputs) overall aln rate, %of all pairs conc
     FAC+="Ilmn_PctAln,Ilmn_PctConc,"

     ## ALE: (5 outputs) ALE score, ALE placeAvg, ALE insertAvg, ALE kmerAvg, ALE depthScoreAvg
     FAC+="Ilmn_ALE,Ilmn_ALE_place,Ilmn_ALE_insert,Ilmn_ALE_kmer,Ilmn_ALE_depth,"

     ## LAP: (1 output) lap score
     FAC+="Ilmn_LAP,"

     ## REAPR: (13 outputs) MBS, pctEF, | FCD, frag_cov, lowscore, link, | softclip, collapsed_repeat, readcov,  low perfect cov, | readorientation, ngaps, gaplen
     FAC+="Ilmn_REAPR_MBS,Ilmn_REAPR_pctEF,Ilmn_REAPR_FCD,Ilmn_REAPR_fragcov,Ilmn_REAPR_lowscore,Ilmn_REAPR_link,Ilmn_REAPR_softclip,Ilmn_REAPR_collapsed,Ilmn_REAPR_readcov,Ilmn_REAPR_lowperfectcov,Ilmn_REAPR_readori,Ilmn_REAPR_ngaps,Ilmn_REAPR_gaplen,"

     ## FRC: (7 outputs) num feat, norm num feat | pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wrong contig, pct of mapped as singleton
     FAC+="Ilmn_FRC_numfeat,Ilmn_FRC_normnumfeat,"
     FAC+="Ilmn_FRC_proper,Ilmn_FRC_wrongdist,Ilmn_FRC_wrongori,Ilmn_FRC_wrongcontig,Ilmn_FRC_singleton,"
     
     ## PILON: (9 outputs) pct confirmed, pct Gsize confirmed, nSNPs/Mb, nSmallINS/Mb, nSmallDEL/Mb, n SNPs+INSs+DELs found per Mb, total num alt alleles in VCF per Mb, Pct of asm flagged as Large collapsed region, Num LCRs per Mb
     FAC+="Ilmn_Pilon_PctConfirmed,Ilmn_Pilon_pctGsizeConfirmed,Ilmn_Pilon_snps,Ilmn_Pilon_ins,Ilmn_Pilon_del,Ilmn_Pilon_all,Ilmn_Pilon_totalAlt,Ilmn_Pilon_pctAsmLenLCR,Ilmn_Pilon_nLCRperMb,"

     #PACBIO (46 outputs total)
     #PacBio - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns,
     FAC+="PacBio_pctalnlen,PacBio_pctmatch_allreads,PacBio_pctmatch_alnreads,PacBio_pctmatch_alns,"

     ## pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, 
     FAC+="PacBio_PctAln,PacBio_alnratio,PacBio_pct_multi_ctg,PacBio_avgalnlen_allaln,"

     ## avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, 
     FAC+="PacBio_avgalnlen_allalnread,PacBio_avgalnlen_allread,PacBio_avgmapq_per_aln,PacBio_avgmapq_per_read,"

     ## avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     FAC+="PacBio_avgmapq_with_unmap,PacBio_avgalnscore_per_aln,PacBio_avgalnscore_aln_reads,PacBio_avgalnscore_allreads,"

     ## SNIFFLES: (7 outputs) ${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}
     FAC+="PacBio_Sniffles_Del,PacBio_Sniffles_Dup,PacBio_Sniffles_Ins,PacBio_Sniffles_Inv,PacBio_Sniffles_Tra,PacBio_Sniffles_total,PacBio_Sniffles_len,"

     ## ALE LONG: (4 outputs) 1=ALE score, 4=placeAvg, 6=kmerAvg, 7=depthScoreAvg
     FAC+="PacBio_ALE,PacBio_ALE_place,PacBio_ALE_kmer,PacBio_ALE_depth,"

     ## FRC LONG: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wron...
     FAC+="PacBio_FRC_numfeat,PacBio_FRC_normnumfeat,"
     FAC+="PacBio_FRC_proper,PacBio_FRC_wrongdist,PacBio_FRC_wrongori,PacBio_FRC_wrongcontig,PacBio_FRC_singleton,"

     ## REAPR LONG HIGH COV: (12 outputs) 1=mbs, 2=pctEF, | 3=FCD, 4=frag_cov, 5=lowscore, 6=link, | 7=softclip, 8=collapsed repeat, 9=readcov, | 10=readorientation, 11=ngaps, 12=gaplen
     FAC+="PacBio_REAPR_MBS,PacBio_REAPR_pctEF,PacBio_REAPR_FCD,PacBio_REAPR_fragcov,PacBio_REAPR_lowscore,PacBio_REAPR_link,PacBio_REAPR_softclip,PacBio_REAPR_collapsed,PacBio_REAPR_readcov,PacBio_REAPR_readori,PacBio_REAPR_ngaps,PacBio_REAPR_gaplen,"


     # NANOPORE (45 outputs total)
     #ONT - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns,
     FAC+="ONT_pctalnlen,ONT_pctmatch_allreads,ONT_pctmatch_alnreads,ONT_pctmatch_alns,"

     ## pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, 
     FAC+="ONT_PctAln,ONT_alnratio,ONT_pct_multi_ctg,ONT_avgalnlen_allaln,"

     ## avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, 
     FAC+="ONT_avgalnlen_allalnread,ONT_avgalnlen_allread,ONT_avgmapq_per_aln,ONT_avgmapq_per_read,"

     ## avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     FAC+="ONT_avgmapq_with_unmap,ONT_avgalnscore_per_aln,ONT_avgalnscore_aln_reads,ONT_avgalnscore_allreads,"

     ## SNIFFLES: (7 outputs) ${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}
     FAC+="ONT_Sniffles_Del,ONT_Sniffles_Dup,ONT_Sniffles_Ins,ONT_Sniffles_Inv,ONT_Sniffles_Tra,ONT_Sniffles_total,ONT_Sniffles_len,"

     ## ALE LONG: (4 outputs) 1=ALE score, 4=placeAvg, 6=kmerAvg, 7=depthScoreAvg
     FAC+="ONT_ALE,ONT_ALE_place,ONT_ALE_kmer,ONT_ALE_depth,"

    ## FRC LONG: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wron...
     FAC+="ONT_FRC_numfeat,ONT_FRC_normnumfeat,"
     FAC+="ONT_FRC_proper,ONT_FRC_wrongdist,ONT_FRC_wrongori,ONT_FRC_wrongcontig,ONT_FRC_singleton,"

     ##	REAPR LONG LOW COV: (11 outputs - fragcov removed) 1=mbs, 2=pctEF, | 3=FCD, 4=lowscore, 5=link, | 6=softclip, 7=collapsed repeat, 8=readcov, | 9=readorientation, 10=ngaps, 11=gaplen
     FAC+="ONT_REAPR_MBS,ONT_REAPR_pctEF,ONT_REAPR_FCD,ONT_REAPR_lowscore,ONT_REAPR_link,ONT_REAPR_softclip,ONT_REAPR_collapsed,ONT_REAPR_readcov,ONT_REAPR_readori,ONT_REAPR_ngaps,ONT_REAPR_gaplen,"

     #BioNano (9 outputs) score, span, total_cov, num, avgscore=Score/Num, pctAsmSpanned=Span/Asm, AvgCovPerBaseInAsm=Cov/Asm, AvgCovPerSpannedBase=Cov/Span, AvgMolSize=Cov/NumAln
     FAC+="BioNano_score,BioNano_span,BioNano_cov,BioNano_num,BioNano_avgscore,BioNano_pctAsmSpan,BioNano_avgCovPerBaseInAsm,BioNano_avgCovPerSpannedBase,BioNano_avgAlignedMolSize,"

     # Old Busco: (1 outout) number complete
     FAC+="Buscov1_numComplete,"

     # New Busco - (6 clades, 2 outputs each, 12 outputs total) look at both Num Complete  and   NumComplete+NumFrags
     #eukaryota metazoa arthropoda insecta endopterygota diptera
     FAC+="Buscov3_euk_C,Buscov3_euk_CF,"
     FAC+="Buscov3_met_C,Buscov3_met_CF,"
     FAC+="Buscov3_art_C,Buscov3_art_CF,"
     FAC+="Buscov3_ins_C,Buscov3_ins_CF,"
     FAC+="Buscov3_end_C,Buscov3_end_CF,"
     FAC+="Buscov3_dip_C,Buscov3_dip_CF,"

     # Known Seq: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     FAC+="Known_blastn_identity,Known_blastn_bit,Known_blastn_alnlen,Known_blastn_hits,"

     # Transcriptome (8 outputs total)
     FAC+="Trans_blastn_identity,Trans_blastn_bit,Trans_blastn_alnlen,Trans_blastn_hits,"
     FAC+="Trans_tblastx_identity,Trans_tblastx_bit,Trans_tblastx_alnlen,Trans_tblastx_hits,"

     # Dmel (12 outputs total)
     FAC+="Dmel_blastn_identity,Dmel_blastn_bit,Dmel_blastn_alnlen,Dmel_blastn_hits,"
     FAC+="Dmel_tblastx_identity,Dmel_tblastx_bit,Dmel_tblastx_alnlen,Dmel_tblastx_hits,"
     FAC+="Dmel_tblastn_identity,Dmel_tblastn_bit,Dmel_tblastn_alnlen,Dmel_tblastn_hits,"

     # Anopheles (12 outputs total)
     FAC+="Agam_blastn_identity,Agam_blastn_bit,Agam_blastn_alnlen,Agam_blastn_hits,"
     FAC+="Agam_tblastx_identity,Agam_tblastx_bit,Agam_tblastx_alnlen,Agam_tblastx_hits,"
     FAC+="Agam_tblastn_identity,Agam_tblastn_bit,Agam_tblastn_alnlen,Agam_tblastn_hits,"

     # RNA-seq (10 outputs) pct_pairs_not_same_ctg, pct_pairs_with_1-2mappedmates_not_same_ctg, pct_pairs_bothmapped_not_same_ctg, pct_pairs_not_same_ctg_qge2, 
     FAC+="RNA_pct_pairs_not_same_ctg,RNA_pct_pairs_with_1-2mappedmates_not_same_ctg,RNA_pct_pairs_bothmapped_not_same_ctg,RNA_pct_pairs_not_same_ctg_qge2,"
     ## pct_bothmapped_not_same_ctg_qge2, pct_bothmappeddiff_with_qge2, avg_min_mapq_all_pairs, avg_min_mapq_pairs_bothmatesmapped, 
     FAC+="RNA_pct_bothmapped_not_same_ctg_qge2,RNA_pct_bothmappeddiff_with_qge2,RNA_avg_min_mapq_all_pairs,RNA_avg_min_mapq_pairs_bothmatesmapped,"
     ## overall_alignment_rate, pct_of_mapped_pairs_that_are_concordant
     FAC+="RNA_PctAln,RNA_PctMappedConc"

     #END
     echo $FAC | awk '{gsub(/,/,"\n"); print}'
}




########################################################################## Categories ##################################################################################
function awko {
 i=1
 while [ $i -le $2 ]; do echo $1 | awk 'OFS="\t" {print $1,$2,$3,$4}' ; let i++ ; done
}

function getcategories {
     ## 201 outputs total
     ##HEADER: outputs 4 columns: broad category, less broad category, was part of basic analysis, was part of thesis analysis
     awko "cat1 cat2 basic thesis" 1

     ##SizeL (5 outputs) ncontigs, maxlen, ng50, lg50, eg
     awko "contiguity contiguity 0 1" 2
     awko "contiguity contiguity 1 1" 1
     awko "contiguity contiguity 0 1" 2

     ##short read: (37 outputs total)
     ## bt2: (2 outputs) overall aln rate, %of all pairs conc
     awko "ilmn ilmn_bt2 0 1" 1
     awko "ilmn ilmn_bt2 0 0" 1

     ## ALE: (5 outputs) ALE score, ALE placeAvg, ALE insertAvg, ALE kmerAvg, ALE depthScoreAvg
     awko "ilmn ilmn_ale 1 1" 1
     awko "ilmn ilmn_ale 0 0" 4 

     ## LAP: (1 output) lap score
     awko "ilmn ilmn_lap 1 1" 1 

     ## REAPR: (13 outputs) MBS, pctEF, | FCD, frag_cov, lowscore, link, | softclip, collapsed_repeat, readcov,  low perfect cov, | readorientation, ngaps, gaplen
     awko "ilmn ilmn_reapr 0 1" 1
     awko "ilmn ilmn_reapr 1 1" 1 
     awko "ilmn ilmn_reapr 0 0" 11

     ## FRC: (7 outputs) num feat, norm num feat | pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wrong contig, pct of mapped as singleton
     awko "ilmn ilmn_frc 1 1" 1
     awko "ilmn ilmn_frc 0 1" 1
     awko "ilmn ilmn_frc 0 0" 5
     
     ## PILON: (9 outputs) pct confirmed, pct Gsize confirmed, nSNPs/Mb, nSmallINS/Mb, nSmallDEL/Mb, n SNPs+INSs+DELs found per Mb, total num alt alleles in VCF per Mb, Pct of asm flagged as Large collapsed region, Num LCRs per Mb
     awko "ilmn ilmn_pilon 0 0" 9

     #PACBIO (46 outputs total)
     #PacBio - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns,
     ## pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, 
     ## avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, 
     ## avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     awko "pacbio pacbio_alnstats 0 0" 4
     awko "pacbio pacbio_alnstats 0 1" 2
     awko "pacbio pacbio_alnstats 0 0" 10

     ## SNIFFLES: (7 outputs) ${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}
     awko "pacbio pacbio_sniffles 0 0" 5
     awko "pacbio pacbio_sniffles 0 1" 2

     ## ALE LONG: (4 outputs) 1=ALE score, 4=placeAvg, 6=kmerAvg, 7=depthScoreAvg
     awko "pacbio pacbio_ale 0 0" 4

     ## FRC LONG: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wron...
     awko "pacbio pacbio_frc 0 0" 7

     ## REAPR LONG HIGH COV: (12 outputs) 1=mbs, 2=pctEF, | 3=FCD, 4=frag_cov, 5=lowscore, 6=link, | 7=softclip, 8=collapsed repeat, 9=readcov, | 10=readorientation, 11=ngaps, 12=gaplen
     awko "pacbio pacbio_reapr 0 0" 12

     # NANOPORE (45 outputs total)
     #ONT - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns,
     ## pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, 
     ## avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, 
     ## avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     awko "ont ont_alnstats 0 0" 4
     awko "ont ont_alnstats 0 1" 2
     awko "ont ont_alnstats 0 0" 10

     ## SNIFFLES: (7 outputs) ${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}
     awko "ont ont_sniffles 0 0" 5
     awko "ont ont_sniffles 0 1" 2

     ## ALE LONG: (4 outputs) 1=ALE score, 4=placeAvg, 6=kmerAvg, 7=depthScoreAvg
     awko "ont ont_ale 0 0" 4

     ## FRC LONG: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wron...
     awko "ont ont_frc 0 0" 7

     ##	REAPR LONG LOW COV: (11 outputs - fragcov removed) 1=mbs, 2=pctEF, | 3=FCD, 4=lowscore, 5=link, | 6=softclip, 7=collapsed repeat, 8=readcov, | 9=readorientation, 10=ngaps, 11=gaplen
     awko "ont ont_reapr 0 0" 11

     #BioNano (9 outputs) score, span, total_cov, num, avgscore=Score/Num, pctAsmSpanned=Span/Asm, AvgCovPerBaseInAsm=Cov/Asm, AvgCovPerSpannedBase=Cov/Span, AvgMolSize=Cov/NumAln
     awko "bng bng_grp1 0 1" 4 
     awko "bng bng_grp2 0 0" 5 

     # Old Busco: (1 outout) number complete
     awko "busco busco_v1 1 1" 1

     # New Busco - (6 clades, 2 outputs each, 12 outputs total) look at both Num Complete  and   NumComplete+NumFrags
     #eukaryota metazoa arthropoda insecta endopterygota diptera
     awko "busco busco_v3 0 0" 12

     # Known Seq: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     awko "known known 0 0" 4

     # Transcriptome (8 outputs total)
     awko "trans trans_blastn 0 0" 4
     awko "trans trans_tblastx 0 0" 4

     # Dmel (12 outputs total)
     awko "dmel dmel_blastn 0 0" 4
     awko "dmel dmel_tblastx 0 0" 4
     awko "dmel dmel_tblastn 0 0" 4

     # Anopheles (12 outputs total)
     awko "agam agam_blastn 0 0" 4
     awko "agam agam_tblastx 0 0" 4
     awko "agam agam_tblastn 0 0" 4

     # RNA-seq (10 outputs) pct_pairs_not_same_ctg, pct_pairs_with_1-2mappedmates_not_same_ctg, pct_pairs_bothmapped_not_same_ctg, pct_pairs_not_same_ctg_qge2, 
     ## pct_bothmapped_not_same_ctg_qge2, pct_bothmappeddiff_with_qge2, avg_min_mapq_all_pairs, avg_min_mapq_pairs_bothmatesmapped, 
     ## overall_alignment_rate, pct_of_mapped_pairs_that_are_concordant
     awko "rna rna_diffctg 0 0" 6
     awko "rna rna_mapq 0 0" 2
     awko "rna rna_bt2 0 0" 2

     #END
}


##########################################################################DEBUG##################################################################################


function debugmsg {
  echo $@ | awk 'OFS="\t" {print $1,$2,$3, $2-$3}'
}

function debug_counts {
     ## Expect: 201 outputs total
     echo name nmetrics expected difference | awk 'OFS="\t" {print $1,$2,$3,$4}'

     ##SizeL (5 outputs) ncontigs, maxlen, ng50, lg50, eg
     SIZE=`getsizestats_vizmat $line | wc -l`
     debugmsg size $SIZE 5 

     ##short read: (37 outputs total)
     ## bt2: (2 outputs) overall aln rate, %of all pairs conc, % of all pairs disc
     BT2=`getbowtie2_simple $SHORT/${b}/mreads/*.err | wc -l`
     debugmsg shortbt2 $BT2 3

     ## ALE: (5 outputs) ALE score, ALE placeAvg, ALE insertAvg, ALE kmerAvg, ALE depthScoreAvg
     ALE=`getale_simple $SHORT/${b}/ale/*ALE.txt | wc -l`
     debugmsg shortale $ALE 5

     ## LAP: (1 output) lap score
     LAP=`getlap $SHORT/${b}/lap/*.lapscore | wc -l`
     debugmsg shortlap $LAP 1

     ## REAPR: (13 outputs) MBS, pctEF, FCD, frag_cov, lowscore, link, softclip, collapsed repeat, readcov, low perfect cov, readorientation, ngaps, gaplen
     REAPR=`getreapr_simple $SHORT/${b}/reapr/output_directory/ | wc -l`
     debugmsg shortreapr $REAPR 13

     ## FRC: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wrong contig, pct of mapped as singleton
     FRC=`getfrc_simple $SHORT/${b}/frc/*gff $SHORT/${b}/frc/*frc_assemblyTable.csv | wc -l`
     debugmsg shortfrc $FRC 7

     ## PILON: (9 outputs) pct confirmed, pct Gsize confirmed, nSNPs/Mb, nSmallINS/Mb, nSmallDEL/Mb, n SNPs+INSs+DELs found per Mb, total num alt alleles in VCF per Mb, Pct of asm flagged as Large collapsed region, Num LCRs per Mb
     PILON=`getpilon ${SHORT}/${b}/pilonvars/ 292000000 | wc -l`
     debugmsg shortpilon $PILON 9

     #PACBIO (45 outputs total)
     #PacBio - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) pctalnlen, pctmatch_allreads, pctmatch_alnreads, pctmatch_alns, pct_reads_that_aln, alnratio, pct_multi_ctg, avgalnlen_allaln, avgalnlen_allalnread, avgalnlen_allread, avgmapq_per_aln, avgmapq_per_read, avgmapq_with_unmap, avgalnscore_per_aln, avgalnscore_aln_reads, avgalnscore_allreads
     ALN=`getlongreadstats_simple ${LONG}/${b}/mreads/pacbio*simple.txt | wc -l`
     debugmsg pacbioaln $ALN 16

     ## SNIFFLES: (7 outputs) ${DEL} ${DUP} ${INS} ${INV} ${TRA} ${NUM} ${SUM}
     SNIF=`getsniffles_simple ${LONG}/${b}/sniffles_pb | wc -l`
     debugmsg pacbiosnif $SNIF 7

     ## ALE LONG: (4 outputs) 1=ALE score, 4=placeAvg, 6=kmerAvg, 7=depthScoreAvg
     ALE=`getale_simple_long $LONG/${b}/ale_pb/*ALE.txt | wc -l`
     debugmsg pacbioale $ALE 4

     ## FRC LONG: (7 outputs) num feat, norm num feat, pct of mapped that are proper, pct of mapped that have wrong dist, pct of mapped that have wrong ori, pct of mapped that have wron...
     FRC=`getfrc_simple $LONG/${b}/frc_pb/*gff $LONG/${b}/frc_pb/*_assemblyTable.csv | wc -l`
     debugmsg pacbiofrc $FRC 7

     ## REAPR LONG HIGH COV: (12 outputs) 1=mbs, 2=pctEF, 3=FCD, 4=frag_cov, 5=lowscore, 6=link, 7=softclip, 8=collapsed repeat, 9=readcov, 10=readorientation, 11=ngaps, 12=gaplen
     REAPR=`getreapr_simple_longread $LONG/${b}/reapr_pb/output_directory/ | wc -l`
     debugmsg pacbioreapr $REAPR 12

     # NANOPORE (44 outputs total)
     #ONT - aln stats, sniffles, ALE, FRC, REAPR
     ## ALN STATS: (16 outputs) - same as above
     ALN=`getlongreadstats_simple ${LONG}/${b}/mreads/ont*simple.txt | wc -l`
     debugmsg ontaln $ALN 16

     ##	SNIFFLES: (7 outputs) -- same as above
     SNIF=`getsniffles_simple ${LONG}/${b}/sniffles_ont | wc -l`
     debugmsg ontsnif $SNIF 7

     ## ALE LONG: (4 outputs) -- same as above
     ALE=`getale_simple_long $LONG/${b}/ale_ont/*ALE.txt | wc -l`
     debugmsg ontale $ALE 4

     ## FRC LONG: (7 outputs) -- same as above
     FRC=`getfrc_simple $LONG/${b}/frc_ont/*gff $LONG/${b}/frc_ont/*_assemblyTable.csv | wc -l`
     debugmsg ontfrc $FRC 7

     ##	REAPR LONG LOW COV: (11 outputs - fragcov removed) 1=mbs, 2=pctEF, 3=FCD, 4=lowscore, 5=link, 6=softclip, 7=collapsed repeat, 8=readcov, 9=readorientation, 10=ngaps, 11=gaplen
     REAPR=`getreapr_simple_longread_lowcov $LONG/${b}/reapr_ont/output_directory/ | wc -l`
     debugmsg ontreapr $REAPR 11

     #BioNano (9 outputs) score, span, total_cov, num, avgscore=Score/Num, pctAsmSpanned=Span/Asm, AvgCovPerBaseInAsm=Cov/Asm, AvgCovPerSpannedBase=Cov/Span, AvgMolSize=Cov/NumAln
     BNG=`getmaligner $b | wc -l`
     debugmsg bionano $BNG 9

     # Old Busco: (1 outout) number complete
     OLDBUS=`getbusco_simple $SHORT/${b}/busco/*/short* | wc -l`
     debugmsg oldbusco $OLDBUS 1
     
     # New Busco - (6 clades, 2 outputs each, 12 outputs total) look at both Num Complete  and   NumComplete+NumFrags
     NEWBUS=`for lineage in eukaryota metazoa arthropoda insecta endopterygota diptera; do getbuscov3_simple ${BUSCOV3}/${b}/${lineage}/*/short* ; done | wc -l`
     debugmsg newbusco $NEWBUS 12

     # Known Seq: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     KNOWNANS=`getblastanalysis ${KNOWN}/${b}/analysis/total.txt | wc -l`
     debugmsg known $KNOWNANS 4

     # Transcriptome (8 outputs total)
     ## BLASTN: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     TRANS1ANS=`getblastanalysis ${TRANSCRIPTOME}/${b}/analysis/total.txt | wc -l`
     ##TBLASTX: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     TRANS2ANS=`getblastanalysis ${TRANSCRIPTOME}/${b}/analysis-tbx/total.txt | wc -l`
     TRANSANS=`echo $TRANS1ANS $TRANS2ANS | awk '{print $1+$2}'`
     debugmsg trans $TRANSANS 8

     # Dmel (12 outputs total)
     ##	BLASTN:	(4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     DMEL1ANS=`getblastanalysis ${DMEL}/${b}/analysis/total.txt | wc -l`
     ##TBLASTX:	(4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     DMEL2ANS=`getblastanalysis ${DMEL}/${b}/analysis-tbx/total.txt | wc -l`
     ##TBLASTN:	(4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     DMEL3ANS=`getblastanalysis ${DMELPEP}/${b}/analysis/total.txt | wc -l`
     DMELANS=`echo $DMEL1ANS $DMEL2ANS $DMEL3ANS | awk '{print $1+$2+$3}'`
     debugmsg dmel $DMELANS 12

     # Anopheles (12 outputs total)
     ## BLASTN: (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     MOS1ANS=`getblastanalysis ${MOSQ}/${b}/analysis/total.txt | wc -l`
     ##TBLASTX:        (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     MOS2ANS=`getblastanalysis ${MOSQ}/${b}/analysis-tbx/total.txt | wc -l`
     ##TBLASTN:        (4 outputs) identity_len,    sum_bit,   len,     num_queries_ge_1_hit
     MOS3ANS=`getblastanalysis ${MOSQPEP}/${b}/analysis/total.txt | wc -l`
     MOSANS=`echo $MOS1ANS $MOS2ANS $MOS3ANS | awk '{print $1+$2+$3}'`
     debugmsg mos $MOSANS 12

     # RNA-seq (10 outputs) pct_pairs_not_same_ctg, pct_pairs_with_1-2mappedmates_not_same_ctg, pct_pairs_bothmapped_not_same_ctg, pct_pairs_not_same_ctg_qge2, pct_bothmapped_not_same_ctg_qge2, pct_bothmappeddiff_with_qge2, avg_min_mapq_all_pairs, avg_min_mapq_pairs_bothmatesmapped, overall_alignment_rate, pct_of_mapped_pairs_that_are_concordant
     RNAANS=`getrnaseqinfo ${RNA}/${b}/mreads/info.txt | wc -l`
     debugmsg rna $RNAANS 10

     #END
}




##TODO-continue to update as above is updated - jun19
function debugsimple {
    line=$1
    b=$2
    SIZE=`getsizestats_vizmat $line | wc `
    BUSCO=`getbusco_simple $SHORT/${b}/busco/*/short* | wc -l`
    BT2=`getbowtie2_simple $SHORT/${b}/mreads/*.err | wc -l`
    ALE=`getale_simple $SHORT/${b}/ale/*ALE.txt | wc -l`
    LAP=`getlap $SHORT/${b}/lap/*.lapscore | wc -l`
    REAPR=`getreapr_simple $SHORT/${b}/reapr/output_directory/ | wc -l`
    FRC=`getfrc_simple $SHORT/${b}/frc/*gff $SHORT/${b}/frc/*frc_assemblyTable.csv | wc -l`
    BNG=`getmaligner $b | wc -l`
    MAPONT=`getlongreadstats_simple ${LONG}/${b}/mreads/ont*simple.txt | wc -l`
    SNIFFONT=`getsniffles_simple ${LONG}/${b}/sniffles_ont | wc -l`
    MAPPB=`getlongreadstats_simple ${LONG}/${b}/mreads/pacbio*simple.txt | wc -l`
    SNIFFPB=`getsniffles_simple ${LONG}/${b}/sniffles_pb | wc -l`

    for var in SIZE BUSCO BT2 ALE LAP REAP FRC BNG SNIFFONT SNIFFPB; do
        echo -e $var"\t"${!var}
    done
    
}




function getnames_simple {
    echo $SIMPLENAMES
}

function getfactors_simple {
    echo $SIMPLEFACTORS
}

function main_debug {
    ##echo DEBUG
    while read line; do
        b=`basename $line .fasta`
        ## > del 2>&1 
        debug_counts $line $b > $OUTDIR/$b.debug.tsv
    done < $FOFN
}

function main_getmetrics {
    while read line; do
        if [ ${line: -3} == ".fa" ]; then b=`basename $line .fa`
        elif [ ${line: -6} == ".fasta" ]; then b=`basename $line .fasta`
        fi
        ## > del 2>&1 
        #getsniffles_simple ${LONG}/${b}/sniffles_ont
        getsimple $line $b > $OUTDIR/${b}.scores.txt
    done < $FOFN
    cd ${OUTDIR}
    ls *.scores.txt | awk '{sub(/\.scores\.txt/,""); print}' | paste -sd"\t" - | cat - <(paste *txt) > all.scores.tab
    cd ../
}


function main_test {
    ## this function is modified as needed to test things
    while read line; do
        b=`basename $line .fasta`
        getpilon ${SHORT}/${b}/pilonvars/ 292000000 ${EXPGENSIZE}
        getpilon ${SHORT}/${b}/pilonvars/ ${EXPGENSIZE}
    done < $FOFN
}

function main {
    OUTDIR=tables
    if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi
    if $TEST; then
        main_test ## this function is modified as needed to test things
    elif $ALL; then
        verbose METRICS....
        main_getmetrics
        verbose DEBUG....
        main_debug
        verbose NAMES...
        getnames > ${OUTDIR}/names.text
        verbose FACTORS....
        getfactors > ${OUTDIR}/factors.text
        verbose CATEGORIES....
        getcategories > ${OUTDIR}/categories.text
        verbose combining descriptors....
        paste <( echo names | cat - ${OUTDIR}/names.text) ${OUTDIR}/categories.text <( echo scale | cat - ${OUTDIR}/factors.text) > ${OUTDIR}/descriptors.text
        verbose "ALL DONE!"
    elif $DEBUG; then
        main_debug
    elif $NAMES ; then
        getnames
    elif $FACTORS ; then
        getfactors
    elif $CATEGORIES ; then
        getcategories
    else 
        main_getmetrics
    fi
}


### EXECUTE
main
