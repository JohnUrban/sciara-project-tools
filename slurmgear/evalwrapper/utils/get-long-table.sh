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

    Alt Usage: bash $0 FOFN all
    ...writing the word 'all' as the second argument will create the default longtables, the vizmat longtables, and the debug files in one pass.
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
    grep "aligned concordantly exactly 1 time" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'
    grep "aligned concordantly >1 times" $F | awk '{print $2}' | awk '{sub(/%/,""); sub(/\(/,""); sub(/\)/,""); print}'
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

function getale {
    F=$SHORT/${1}/ale/*ALE.txt
    awk 'NR==1 || NR==4 || NR==5 || NR==6 || NR==7 || NR==8 || NR==9 || NR==10 || NR==11 || NR==12 {print $3}' $F
}

function getale_vizmat {
    F=$SHORT/${1}/ale/*ALE.txt
    awk 'NR==1 || NR==4 || NR==5 || NR==6 || NR==7 || NR==10 || NR==11 || NR==12 {print $3}' $F
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
    awk 'NR==2 {$2, $20, $3, $21, $4, $22, $5, $23, $6, $24}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
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
    awk 'NR==2 {$21, $23}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
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
    D=$LONG/$1/snifflestats
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
    elif [ $NARG -eq 2 ] && [ $DEBUG == "all" ]; then
        ##echo ALL
        while read line; do
            b=`basename $line .fasta`
            getall $line $b > $LONGTABLE/$b.longtable
            getall_vizmat $line $b > $LONGTABLE/$b.longtable.vizmat
            debugmode $line $b > $LONGTABLE/$b.longtable.debugmode
        done < $FOFN
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
