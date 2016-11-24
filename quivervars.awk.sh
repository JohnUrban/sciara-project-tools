#!/bin/bash

function msg {
  echo "Usage:
	quivervar.awk.sh pilon.changes

	Output:
		numMisMatch numInsertion numDeletion TotalCount Total(NonHeader)LinesSeen propMisMatch propIns propDel

	Uses TotalCount to get proportions. TotalLinesSeen is provided as a QC that TotalCount equals TotalLinesSeen.
"
}

if [ $# -eq 0 ]; then msg; exit; fi


awk 'OFS="\t" {if ( $1 ~ /^#/ ) h+=1; else if($3 == "substitution") m+=1; else if($3 == "insertion") i+=1; else if($3 == "deletion") d+=1} END{print m,i,d,m+i+d,NR-h,m/(m+i+d),i/(m+i+d),d/(m+i+d)}' $f
