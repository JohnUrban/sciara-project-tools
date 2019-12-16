#!/bin/bash

##DEFAULTS
MINCOV=0 ## Recommend filtering for >= 25X
MINQV=0  ## Recommend trying >= 20  (1 in 100 expected to be wrong) for higher quality motif analyses (this does not apply to modified_base)
MIN10LOGP=0 ## Recommend filtering for >= 20 corresponding to p-valueof 0.01 -- this should actually already be lowest score in the file if used 0.01 as p-val cutoff
		## However -- somehow there are -log10(p)s < 20...
		## Yet anything labeled only as modified_base is > 20
		## Therefore, all candidate mods are prob >20 before identification, and it becomes less so later....
			## I can see that the -log10(p) and identificationQv are NOT the same -- so it is not just simply replaced...
FASTA=""

## HELP
if [ $# -eq 0 ]; then echo -e " midbase.sh path-to=mobase.gff mod-type \n\t ...where modtpype in m5C, m4C, 6mA, and modified_base \n \t Optional args trailing: MINCOV, MINQV, MIN10LOGP BACKGROUND_FASTA\n\t Background calculated from the 41 bp centered on each mod by default. Provide FASTA for diff bg. It will take kmers from both strands. \n"; exit; fi

## GET VARS
GFF=$1
MOD=$2
if [ $# -ge 3 ]; then MINCOV=$3; fi
if [ $# -ge 4 ]; then MINQV=$4; fi
if [ $# -ge 5 ]; then MIN10LOGP=$5; fi
if [ $# -ge 6 ]; then FASTA=$6; fi

SCRIPTSDIR=`dirname ${0}`

#grep -v ^# $GFF | awk -v mod=$MOD '$3==mod {gsub(/[=|;]/,"\t"); print $12}' | python ${SCRIPTSDIR}/midbase.py 
#exit

if [ ${MOD} == "modified_base" ]; then 
  grep -v ^# $GFF | awk -v mod=$MOD '$3==mod {gsub(/[=|;]/,"\t"); print}' | awk -v "C=${MINCOV}" -v "P=${MIN10LOGP}" '$10>=C && $6>=P {print $12}' | python ${SCRIPTSDIR}/midbase.py ${FASTA}
else
  grep -v ^# $GFF | awk -v mod=$MOD '$3==mod {gsub(/[=|;]/,"\t"); print}' | awk -v "C=${MINCOV}" -v "Q=${MINQV}" -v "P=${MIN10LOGP}" '$10>=C && $6>=P && $NF>=Q {print $12}' | python ${SCRIPTSDIR}/midbase.py ${FASTA}
fi


exit

## NOTES ON GFF FILE BELOW
#Column	Description
#Seqid	Reference sequence tag for this observation. Same as refName in the .csv file.
#Source	Name of tool -- "kinModCall".
#Type	Modification type – a generic tag "modified_base" is used for unidentified bases. For identified bases, m6A, m4C, and m5C are used.
#Start	Location of modification.
#End	Location of modification.
#Score	-10 log (p-value) score for the detection of this event. Analogous to a Phred quality score. 
#		A value of 20 is the minimum default threshold for this file, and corresponds to a p-value of 0.01. A score of 30 corresponds to a p-value of 0.001.
#Strand	Native sample strand where kinetics were generated. “+” is the strand of the original FASTA and “-” is the reverse complement of the strand. 
#		Note that in the .csv file these are marked "0" and "1" respectively.
#Phase	Not applicable.
#Attributes	Contains extra fields. 
#	IPDRatio is traditional IPD Ratio, 
#	context is the reference sequence -20bp to +20bp around the modification plus the base at this location as the 21st character, and 
#		sequencing coverage of that position. 
#		Context is always written in 5’ -> 3’ orientation of the template strand. 
#	If the modification type is determined and coverage at that position is at least 10x, 
#		there are several additional metrics: 
#			an estimate of the fraction of modified reads and the upper and lower 95% confidence intervals are included as ‘frac’, ‘fracLow’, and ‘fracUp’, 
#			and a confidence score for the modification type is reported as identificationQv (calculated like the Score).
#
