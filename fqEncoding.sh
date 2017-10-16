#!/bin/bash
set -e
set -u
set -o pipefail

source /Users/jurban/software/sratoolkit/source.sh

function help {
    echo "
        Modified from MedHat awk example: https://www.biostars.org/p/63225/

        Usage: ${0} -f:n:hd
        -f [required] path to fastq/fq (can be .gz)
        -n number of records to look at [1000]
        Help/Debug
        -h help - returns this message; also returns this when no arguments given
        -d debug
"
}



##############################################################################
## TRIGGER HELP IF NO ARGS
##############################################################################
if [ $# -eq 0 ]; then help; exit; fi


##############################################################################
## DEFAULTS
##############################################################################
DIR=./
HELP=false
DEBUG=false
FASTQ=false
n=1000

##############################################################################
## GET OPTS
##############################################################################
while getopts "f:n:hd" arg; do
    case $arg in 
        f) FASTQ=$OPTARG;;
        n) n=$OPTARG;;
        h) HELP=true;;
        d) DEBUG=true;;
        *) help; exit;;
    esac
done


##############################################################################
## TRIGGER HELP IF OPTED FOR
##############################################################################
if ${HELP} || [ ${FASTQ} == false ] ; then help; exit; fi


##############################################################################
## FUNCTIONS
##############################################################################

function gzipped {
  if [ ${FASTQ: -3} == ".gz" ]; then echo true; else echo false; fi
}

function decode {
  awk '{if(NR%4==0) printf("%s",$0);}' |  \
  od -A n -t u1 | \
  awk 'BEGIN {min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding";}'
}  
  

##############################################################################
## EXECUTE
##############################################################################
N=`echo $n*4 | bc`
if gzipped ; then gunzip -c $FASTQ | head -n $N ; else head -n $N $FASTQ ; fi | decode


