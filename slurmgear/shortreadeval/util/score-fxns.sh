#!/bin/bash

export PATH=/users/jurban/scratch/male-ilmn/long_read_evals/scripts:$PATH

#function checknames {
# F=alescores.${1}.txt
# awk '{print $1}' $F > $F.del;
# c=0
# for f in busco.${1}.txt numfrc.${1}.txt lapscores.${1}.txt reapr-errorfree.${1}.txt reapr-perbasescore.${1}.txt; do
#   awk '{print $1}' $f > $f.del
#   s=`diff $F.del $f.del | wc -l`
#   c=$(( $c+$s ))
#   rm $f.del
# done
# rm $F.del
# echo $c
#}

function checknames {
 all="pctmap.${1}.txt alescores.${1}.txt busco.${1}.txt numfrc.${1}.txt lapscores.${1}.txt reapr-errorfree.${1}.txt reapr-perbasescore.${1}.txt"
 paste $all | cut -f 1,3,5,7,9,11,13 > del
 t=`cat del | wc -l`
 c=0
 for i in {2..7}; do
   n=`awk -v "i=$i" '$1 == $i' del | wc -l`
   s=`echo $t-$n | bc`
   c=$(( $c+$s ))
 done
# rm del
 echo $c
}



function combinescores {
 c=`checknames $1`
 if [ $c -eq 0 ]; then
   paste pctmap.${1}.txt alescores.${1}.txt busco.${1}.txt lapscores.${1}.txt reapr-errorfree.${1}.txt reapr-perbasescore.${1}.txt numfrc.${1}.txt | awk 'OFS="\t" {print $1,$2,$4,$6,$8,$10,$12,$14}' > combined-scores.${1}.txt
 else echo Cannot Combine $1
 fi
}


function checkorder {
## check order over many
## 
#paste ${1}*/combined*txt | cut -f 1,8,15,22,29,36,43,50,57 | sed 's/.quiver/\t/g' | cut -f 1,3,5,7,9,11,13,15,17 > del
#paste $@ | cut -f 1,8,15,22,29,36,43,50,57 | sed 's/.quiver/\t/g' | cut -f 1,3,5,7,9,11,13,15,17 > del
paste $@ | cut -f 1,9,17,25,33,41,49,57,65 | sed 's/.quiver/\t/g' | cut -f 1,3,5,7,9,11,13,15,17,19 > del
 t=`cat del | wc -l`
 c=0
 for i in {2..12}; do
   if [ $i -le $# ]; then
     n=`awk -v "i=$i" '$1 == $i' del | wc -l`
     s=`echo $t-$n | bc`
     c=$(( $c+$s ))
   fi
 done
# rm del
 echo $c
}

#checkx {
#paste $@ | cut -f 1,9,17,25,33,41,49,57,65 | sed 's/.quiver/\t/g' | cut -f 2,4,6,8,10,12,14,16,18,20 > del
#}

#name, pct_ilmn_mapped, ale, busco_pct_complete, lap, reapr errfree, reapr base mean, num frc
#   all="pctmap.${1}.txt alescores.${1}.txt busco.${1}.txt lapscores.${1}.txt reapr-errorfree.${1}.txt reapr-perbasescore.${1}.txt numfrc.${1}.txt" 

function unitescores {
 c=`checkorder $@`
 if [ $c -eq 0 ]; then
   ## Get all pctmap over time
   paste $@ | cut -f 1,2,10,18,26,34,42,50,58,66 > pctmapped.united.txt
   paste $@ | cut -f 1,3,11,19,27,35,43,51,59,67 > ale.united.txt
   paste $@ | cut -f 1,4,12,20,28,36,44,52,60,68 > busco.united.txt
   paste $@ | cut -f 1,5,13,21,29,37,45,53,61,69 > lap.united.txt
   paste $@ | cut -f 1,6,14,22,30,38,46,54,62,70 > reapr-errorfree.united.txt
   paste $@ | cut -f 1,7,15,23,31,39,47,55,63,71 > reapr-mean.united.txt
   paste $@ | cut -f 1,8,16,24,32,40,48,56,64,72 > frc.united.txt
      
 fi
}


function simplyscore {
 getalescores.sh | sort -k1,1 > alescores.${1}.txt
 get-busco.sh | sort -k1,1 > busco.${1}.txt
 getnumfrc.sh | sort -k1,1 > numfrc.${1}.txt
 getlaps.sh | sort -k1,1 > lapscores.${1}.txt
 getreapr-errorfree.sh | sort -k1,1 > reapr-errorfree.${1}.txt
 getreapr-perbase-score.sh | sort -k1,1 > reapr-perbasescore.${1}.txt
 get-pctmap.sh | sort -k1,1 > pctmap.${1}.txt
 combinescores $1
}

function score {
 cd $1
 simplyscore $1
 cd ../ ;
}


function scoredir {
 for D in quiver0 initial_consensus miniasm; do
  if [ -d $D ]; then
   echo $D
   score $D
  fi
 done

 for i in {1..7}; do
  echo quiver${i}
  score quiver${i}
 done
}

