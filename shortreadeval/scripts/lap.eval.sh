#!/bin/bash


PATH=~/software/sciaratools/sciara-project-tools/lap/1.1/aligner:/users/jurban/software/sciaratools/sciara-project-tools/:$PATH

##PROBS
 if [ -f ${BASE}.prob ]; then
   c=`cat ${BASE}.prob | wc -l`
   if [ $c -lt 10 ]; then ##it is just a ghost file, so make it for reals
     calc_prob.py -p $P -a $REF  -q -1 $R1 -2 $R2 -X 800 -I 0 -o fr -m 432 -t 75 -b $BT2 > ${BASE}.prob
   fi
 else ##DoesNotExist so create
  date; echo prob file DoesNotExist so create
  calc_prob.py -p $P -a $REF  -q -1 $R1 -2 $R2 -X 800 -I 0 -o fr -m 432 -t 75 -b $BT2 > ${BASE}.prob
  date; echo prob file DoesNotExist so created it...
 fi
 
## LAPSCORE
 if [ -f ${BASE}.lapscore ]; then
   c=`cat ${BASE}.lapscore | wc -l`
   if [ $c -lt 1 ]; then ##it is just a ghost file, so make it for reals
      sum_prob.py -i ${BASE}.prob > ${BASE}.lapscore
   fi
 else ##DoesNotExist so create
   sum_prob.py -i ${BASE}.prob > ${BASE}.lapscore
 fi
