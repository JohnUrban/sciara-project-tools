#!/bin/bash
module load blasr

## ECHO OUT VARIABLES
for var in BASE ASM PROTOCOL MAKEFAKEQUALS GAPSONLY SPANONLY THREADS RUNSETUP RUNMAPPING RUNSUPPORT RUNEXTRACTION RUNASSEMBLY RUNOUTPUT; do 
   echo -e ${var} ${!var}; 
done

## DEFINE FUNCTIONS
function start {
echo Start ${1} >> JellyApply.log
date >> JellyApply.log
}

function end {
echo End ${1} >> JellyApply.log
date >> JellyApply.log
echo >> JellyApply.log
}

function run {
 #arg1 = step name
 #arg2 = path to protocol.xml
 #arg3 = quotes around any args passed to step -- eg. "-p 2 -x 5"
 step=${1}
 protocol=${2}
 args=${3}
 echo args $args
 echo argsnew ${@:3}
 start ${1}
 echo "Jelly.py ${step} ${protocol} ${@:3} 2>>JellyApply.log"
 Jelly.py ${step} ${protocol} ${@:3} 2>>JellyApply.log 
 end ${1}
}


touch JellyApply.log
if ${MAKEFAKEQUALS}; then
  start "making quals"
  BASE=`basename ${ASM} .fasta` 
  fakeQuals.py ${ASM} ${BASE}.qual
  end "making quals"
fi

## if using fastq, may want the following as a solution
## fastq2faqual.py --fastq /path/to/polished_assembly.fastq --fa --qual --out ${outpre}



## Docment where gaps are
start "documenting gaps before running PBJelly"
scf-N-to-BED.py ${ASM} | awk 'OFS="\t" {print $1,$2,$3,$3-$2}' > gaps.bed
end "documenting gaps before running PBJelly"


## SOME ARGUMENTS
SETUPARGS=""
##MAPARGS=`echo x ${THREADS} | awk '{print "-"$1" \"--nproc="$2"\""}'`
MAPARGS=`echo x ${THREADS} | awk '{print "-"$1" \"-n "$2"\""}'`
SUPPORTARGS=""
EXTRACTARGS=""
ASSEMBLYARGS=`echo x ${THREADS} | awk '{print "-"$1" \"--nproc="$2"\""}'`
OUTPUTARGS=""

TMPSUPP=""
if ${GAPSONLY}; then TMPSUPP+="--capturedOnly "; fi
if ${SPANONLY}; then TMPSUPP+="--spanOnly "; fi
if ${GAPSONLY} || ${SPANONLY}; then SUPPORTARGS+=`echo x ${TMPSUPP} | awk '{print "-"$1" \""$2" "$3"\""}'` ; fi

## ECHO OUT VARIABLES
for var in SETUPARGS MAPARGS SUPPORTARGS EXTRACTARGS ASSEMBLYARGS OUTPUTARGS; do
   echo -e ${var} ${!var}; 
done


## RUNNING
if ${RUNSETUP}; then run setup ${PROTOCOL} ${SETUPARGS}; fi
##if ${RUNMAPPING}; then run mapping ${PROTOCOL} ${MAPARGS} ; fi
if ${RUNMAPPING}; then 
  start mapping; 
  Jelly.py mapping ${PROTOCOL} -x "-n ${THREADS}" 2>>JellyApply.log ;
  end mapping; 
fi

#######if ${RUNSUPPORT}; then run support ${PROTOCOL} ${SUPPORTARGS}; fi
if ${RUNSUPPORT}; then 
    start support
    if ${GAPSONLY} && ${SPANONLY}; then 
       Jelly.py support ${PROTOCOL} -x "--capturedOnly --spanOnly" 2>>JellyApply.log ;
    elif ${GAPSONLY}; then 
      Jelly.py support ${PROTOCOL} -x "--capturedOnly" 2>>JellyApply.log ;
    elif ${SPANONLY}; then 
      Jelly.py support ${PROTOCOL} -x "--spanOnly" 2>>JellyApply.log ;
    else
      Jelly.py support ${PROTOCOL} 2>>JellyApply.log
    fi
    end support
fi
if ${RUNEXTRACTION}; then run extraction ${PROTOCOL} ${EXTRACTARGS}; fi
if ${RUNASSEMBLY}; then run assembly ${PROTOCOL} ${ASSEMBLYARGS} ; fi
if ${RUNOUTPUT}; then run output ${PROTOCOL} ${OUTPUTARGS} ; fi

