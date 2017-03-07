#!/bin/bash

FOFN=input.fofn
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

while read line; do
  b=`basename $line .fasta`
  echo $b
  getsizestats $line
  getbusco $b
  getbowtie2 $b  
  getale $b
  getlap $b
  getreapr $b
  getfrc $b
  echo
done < $FOFN
