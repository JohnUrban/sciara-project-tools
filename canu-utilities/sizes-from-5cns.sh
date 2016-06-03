#!/bin/bash

if [ $# -eq 0 ]; then echo "
Usage:
sh script.sh path-to-unitigging/5-consensus/[--optional--consensus.jobid_]
"
exit
fi

grep "Working on unitig" $1*.out | awk 'match($0, /length\ [0-9]*/) { print substr($0, RSTART, RLENGTH)}' | awk '{print $1"\t"$2}'
