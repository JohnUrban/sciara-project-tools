## Highest prob of needing changing
REF=""					  			## CHANGE
REFWINDOWS=""							## CHANGE - EXAMPLE = "contig_20:0-24000,contig_20:25000-17453741" -- leave blank if none to do --- turns "--identify" on for these
PROBLEM_REFWINDOWS=""						## CHANGE - EXAMPLE = "contig_20:24000-25000" -- leave blank if none to do --- turns "--identify" off for these

## CATCH ISSUES HERE
if [ -z ${REF} ]; then echo "NEED TO SPECIFY REF !!!"; exit; fi
if [ -z ${REFWINDOWS} ] && [ -z ${PROBLEM_REFWINDOWS} ]; then echo "NEED TO SPECIFY REFWINDOWS, PROBLEM_REFWINDOWS or both !!!"; exit; fi

## Optionally change these
QOS=biomed-condo
PVALUE=0.01
IDENTIFY=7
THREADS=1 ## NOTE - this script automatically doubles this number for launching the job if both REFWINDOWS and PROBLEM_REFWINDOWS are non-empty as those 2 jobs run in parallel will both get THREADS (therefore, need 2*THREADS)
TIME=48:00:00
MEM=30g


## ABSPATHING
REF=`readlink -f ${REF}`
REFBASE=`basename ${REF} .fasta`


## VARS THAT MAKE ASSUMPSOINS
TOPDIR=`readlink -f ../../..`					## ASSUMPTION ON FILE HIERARCHY AND WHERE YOU ARE IN IT
INPUT=${TOPDIR}/${REFBASE}/out_all.cmp.h5 			## ASSUMES INPUT SHOULD ALREADY BE IN THIS REF_BASENAME_DIR
KC=`readlink -f ./kinetics-contig-refwindows.sh`		## ASSUMES YOU ARE USING THE KC SCRIPT IN THIS DIR





for var in REF REFBASE REFWINDOWS PROBLEM_REFWINDOWS QOS PVALUE IDENTIFY THREADS TIME MEM TOPDIR INPUT KC; do
  echo ${var} ${!var}
done > vars.cfg

## LAUNCH
#sbatch --mem=${MEM} --time=${TIME} -c ${THREADS} -J ${REFBASE}_${contig}_ipd_redo --qos=${QOS} \
#      --export=PVALUE=${PVALUE},IDENTIFY=${IDENTIFY},THREADS=${THREADS},INPUT=${INPUT},REF=${REF},contigName=${contig} ${KC}


## DOUBLE THREADS COUNT IF RUNNING BOTH JOB TYPES
if [ ! -z ${REFWINDOWS} ] && [ ! -z ${PROBLEM_REFWINDOWS} ]; then THREADS=`echo 2*${THREADS} | bc`; fi

sbatch --mem=${MEM} --time=${TIME} -c ${THREADS} -J ${REFBASE}_ipd_redo --account=${QOS} ${KC}




