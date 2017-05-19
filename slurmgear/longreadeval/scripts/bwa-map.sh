#!/bin/bash

source ~/software/bwa/source.sh

BAM=${TYPE}.bam
bwa mem -t $MTHREADS -M -x $TYPE $BWAIDX $FASTQ | samtools sort -T $TYPE --threads $MTHREADS -o ${BAM}


## stats

## ont2d --> ont
if [ $TYPE == "ont2d" ] && [ -z $PRE ]; then PRE=ont; elif [ -z $PRE ]; then PRE=${TYPE}; fi


## Gets number of alignments -- however each long read can be split into >1 aln, at least with bwa
samtools view -c -F 4 $BAM > ${PRE}.numaln &

##Gets number of uniq read names that aligned
samtools view -F 4 $BAM | awk '{print $1}' | sort | uniq | wc -l > ${PRE}.numuniqaln &


## wait for last 2
wait


## Gets number of alignments + non-aln -- i.e. number of entries in BAM
samtools view -c $BAM > ${PRE}.numentries &

##Gets number of uniq read names that aligned or did not align --- num uniq entries -- should be equivalent to number of reads in input Fastq file -- should be same for all comparisons
samtools view $BAM | awk '{print $1}' | sort | uniq | wc -l > ${PRE}.numuniqentries &

## wait for last 2
wait



## Sums MAPQ regardless of everything
samtools view $BAM | awk '{s+=$5}END{print s}' > ${PRE}.sum.mapq &
wait


## PCTs and RATIOs and Avg MAPQ
## PCTs\
nualn=`cat ${PRE}.numuniqaln`
nuent=`cat ${PRE}.numuniqentries`
pctaln=`echo $nualn/$nuent | bc -l`
echo $pctaln > ${PRE}.pctaln

## RATIO
naln=`cat ${PRE}.numaln`
alnratio=`echo $naln/$nualn | bc -l`
echo $alnratio > ${PRE}.alnratio

## Avg MAPQ
mapqsum=`cat ${PRE}.sum.mapq`
mapqavg=`echo $mapqsum/$naln | bc -l`
echo $mapqavg > ${PRE}.avg.mapq



## Get things on a per read basis....
samtools view ${BAM} | python -c "
import re, sys
from collections import defaultdict
import numpy as np

readlengths = defaultdict(list)
alnlengths = defaultdict(list)
editdists = defaultdict(list)
alnscores = defaultdict(list)
mapqs = defaultdict(list)
contigs = defaultdict(set)
UN_alnscores = defaultdict(list)
UN_mapqs = defaultdict(list)
names = set([])
for line in sys.stdin:
    line = line.strip().split()
    name = line[0]
    names.add(name)
    readlengths[name].append( sum([int(e) for e in re.findall('(\d+)[MIS=X]', line[5])])  )
    if line[1] != "*":
        alnlengths[name].append( sum([int(e) for e in re.findall('(\d+)[MI=X]', line[5])]) )
        editdists[name].append( int(line[11].split(":")[2])  )
        alnscores[name].append( int(line[13].split(":")[2])  )
        mapqs[name].append( line[4]  )
        contigs[name].add( line[2]  )
    else:
        UN_alnscores[name].append( 0  )
        UN_mapqs[name].append( 0  )

    out = open('per-read.txt','w')
    num_0_ctg = 0
    num_1_ctg = 0
    num_multi_ctg = 0
    readlensum = 0
    alnlensum = 0
    mapqsum = 0
    editsum = 0
    mapqsum2 = 0
    editsum2 = 0
    numaln = 0
    numaln_un = 0
    numunaln = 0
    for name in list(names):
        ALN = name in alnlengths
        UN = name in UN_alnscores
        readlen = readlengths[name][0]
        if ALN:
            nctg = len( contigs[name] )
            alnlen = sum( alnlengths[name] )
            edit = sum( editdists[name] )
            alnscore = np.mean( alnscores[name] )
            mapq = np.mean( mapqs[name] )
            if UN:
                alnscore2 = np.mean( alnscores[name]+UN_alnscores[name] )
                mapq2 = np.mean( mapqs[name]+UN_mapqs[name] )
            else:
                alnscore2 = alnscore
                mapq2 = mapq
        elif UN:
            nctg = 0
            alnlen = 0
            edit = readlen
            alnscore = 0
            mapq = 0
            alnscore2 = 0
            mapq2 = 0
        out.write( ("\t").join([str(e) for e in [name, nctg, alnlen, edit, alnscore, mapq, alnscore2, mapq2]]) + '\n' )
        #additional analysis
        if nctg == 0: 
            num_0_ctg += 1
        elif nctg == 1: 
            num_1_ctg += 1
        elif nctg > 1: 
            num_multi_ctg += 1
        readlensum += readlen
        alnlensum += alnlen
        mapqsum += mapq
        editsum += edit
        mapqsum2 += mapq2
        editsum2 += edit2
        if ALN and UN:
            numaln += 1
            numaln_un += 1
        elif ALN:
            numaln += 1
        elif UN:
            numunaln += 1

    #final analysis
    out.close()
    out = open('per-read-stats.txt','w')
    out.write( 'num_0_ctg\t'+str(num_0_ctg)+'\n' )
    out.write( 'num_1_ctg\t'+str(num_1_ctg)+'\n' )
    out.write( 'num_multi_ctg\t'+str(num_multi_ctg)+'\n' )
    out.write( 'readlensum\t'+str(readlensum)+'\n' )
    out.write( 'alnlensum\t'+str(alnlensum)+'\n' )
    out.write( 'mapqsum\t'+str(mapqsum)+'\n' )
    out.write( 'editsum\t'+str(editsum)+'\n' )
    out.write( 'mapqsum2\t'+str(mapqsum2)+'\n' )
    out.write( 'editsum2\t'+str(editsum2)+'\n' )
    out.write( 'numaln\t'+str(numaln)+'\n' )
    out.write( 'numaln_un\t'+str(numaln_un)+'\n' )
    out.write( 'numunaln\t'+str(numunaln)+'\n' )
    out.close()
" 


