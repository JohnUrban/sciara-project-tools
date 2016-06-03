#!/bin/bash

if [ $# -eq 0 ]; then echo "
Usage: sh copy.sh FROMdir TOdir
"
exit
fi


#show-to-dir -- mkdir if not present
if [ -d $2 ]; then echo ls $2; else mkdir $2; ls $2; fi

echo
cp $1/*fasta $2
cp $1/unitigging/0-mercounts/*.ms22.estMerThresh.err $2
cp $1/unitigging/4-unitigger/unitigger.sh $2
ls $2
echo
