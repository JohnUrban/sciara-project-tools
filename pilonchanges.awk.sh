#!/bin/bash

function msg {
  echo "Usage:
	pilonchanges.awk.sh pilon.changes

	Output:
		numMisMatch numInsertion numDeletion TotalCount TotalLinesSeen
"
}

if [ $# -eq 0 ]; then msg; exit; fi



awk 'OFS="\t" {if($3 == "." && $4 != ".") i+=1; else if($3 != "." && $4 == ".") d+=1; else if($3 != "." && $4 != ".") m +=1}END{print m,i,d,m+i+d,NR}' $1
