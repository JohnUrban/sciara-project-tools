#!/bin/bash


if [ $# -eq 0 ] || [ $1 == "-h" ] || [ $1 == "-help" ] || [ $1 == "--help" ]; then
    echo "
    Usage: bash $0 FOFN

    ...where FOFN has list of all assemblies used in the assembly evaluations in subdirs you are trying to summarize.
    (( typically called input.fofn ))

    Alt Usage: bash $0 FOFN debug
    ...writing the word 'debug' as the second argument will create files that have more information.
    "
    exit
fi

##FOFN=input.fofn
FOFN=$1


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
    grep "aligned discordantly 1 time" $F | awk '{print $2}' | awk '{sub(/\ /,"_"); sub(/%/,"\\%"); print}'
    disc1=`grep "aligned discordantly 1 time" $F | awk '{print $1}'`
    pctdisc1=`echo $disc1 $pair | awk '{print 100.0*$1/$2}'`
    echo "${disc1}_(${pctdisc1}\%)"
    grep "pairs aligned 0 times concordantly or discordantly; of these:" $F | awk '{print $1}'
    grep "mates make up the pairs; of these:" $F | awk '{print $1}'
    grep "aligned 0 times" $F | awk '{print $1}'
    grep "aligned exactly 1 time" $F | awk '{print $1}'
    grep "aligned >1 times" $F | awk '{print $1}'
}

function getale {
    F=$SHORT/${1}/ale/*ALE.txt
    awk 'NR==1 || NR==4 || NR==5 || NR==6 || NR==7 || NR==8 || NR==9 || NR==10 || NR==11 || NR==12 {print $3}' $F
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
    awk 'NR==2 {OFS="\t"; print 100*$38/$2, $39+$40+$41+$42, $39, $40, $41, $42, $43+$44+$45+$46+$47+$48+$49, $43, $44, $45, $46, $47, $48, $49, $2, $20, $3, $21, $4, $22, $5, $23, $6, $24}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
    stats=`getasmstats $1`
    bstats=`awk '{gsub(/,/,"\t"); print }' ${D}/broken_assembly.sizestats.csv`
    echo $stats | awk '{print $10}'
    echo $bstats | awk '{print $10}'
    echo $stats | awk '{print $11}'
    echo $bstats | awk '{print $11}'
    echo $stats | awk '{print $12}'
    echo $bstats | awk '{print $12}'
    awk 'NR==2 {OFS="\t"; print $36, $37}' ${D}/05.summary.report.tsv | awk '{gsub(/\t/,"\n"); print}'    
  
}


function getfrc {
    F1=$SHORT/${1}/frc/*gff
    F2=$SHORT/${1}/frc/*frc_assemblyTable.csv
    grep -c -v ^# $F1
    tail -n 1 $F2 | awk '{gsub(/,/,"\n"); print}' | tail -n +3
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
    SNIF=`getsniffles $b | wc `
    for var in SIZE BUSCO BT2 ALE LAP REAP FRC BNG SNIF; do
        echo -e $var"\t"${!var}
    done
}

function main {
    LONGTABLE=longtables
    if [ ! -d $LONGTABLE ]; then mkdir $LONGTABLE; fi
    if [ $# -eq 2 ] && [ $2 == "debug" ]; then
        while read line; do
            b=`basename $line .fasta`
            debugmode $line $b > $LONGTABLE/$b
        done < $FOFN
    else
        while read line; do
            b=`basename $line .fasta`
            getall $line $b > $LONGTABLE/$b
        done < $FOFN
    fi
}


### EXECUTE
main
