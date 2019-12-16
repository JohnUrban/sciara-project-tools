#!/bin/bash

## CHECK TO MAKE SURE template/create-config-and-launch.sh was minimally edited
ANS=`cat template/create-config-and-launch.sh | head | grep ^REF | awk 'NR==1 {print $1}'`
if [ $ANS == "REF=\"\"" ]; then echo "NEED TO MINIMALLY EDIT THE REF IN template/create-config-and-launch.sh" ; exit; fi

## MK DIRS
while read contig; do
  echo ${contig}
  cp -r template ${contig}
done < names.txt 
