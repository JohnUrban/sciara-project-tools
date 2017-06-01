#!/bin/bash

source ~/software/bwa/source.sh
BAM=${TYPE}.bam

echo MTHREADS $MTHREADS
echo TYPE $TYPE
echo BAM $BAM
echo BWAIDX $BWAIDX 
echo FASTQ $FASTQ

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
totalaln = 0
totalunaln = 0
for line in sys.stdin:
    line = line.strip().split()
    name = line[0]
    names.add(name)
    if line[1] != '4':
        ##rlen = sum([int(e) for e in re.findall('(\d+)[MIS=X]', line[5])])
        rlen = sum([int(e) for e in re.findall('(\d+)[MISH=X]', line[5])])
        readlengths[name].append( rlen  )
        totalaln += 1
        alnlengths[name].append( sum([int(e) for e in re.findall('(\d+)[MI=X]', line[5])]) )
        editdists[name].append( int(line[11].split(':')[2])  )
        alnscores[name].append( int(line[13].split(':')[2])  )
        mapqs[name].append( int(line[4])  )
        contigs[name].add( line[2]  )
    else:
        readlengths[name].append( len(line[9]) ) 
        totalunaln += 1
        UN_alnscores[name].append( 0 )
        UN_mapqs[name].append( 0 )

#initate stats vars
summapq_old = 0
num_0_ctg = 0
num_1_ctg = 0
num_multi_ctg = 0
readlensum = 0
alnlensum = 0
alnscoresum = 0
alnscoresum2 = 0
mapqsum = 0
mapqsum2 = 0
editsum = 0
editsumwithunmapped = 0
numaln = 0
numaln_un = 0
numunaln = 0
unmappedreadlensum = 0
matchsum = 0
matchsum2 = 0
totaleditsum = 0
unalnportionsum = 0
sumalnscore_per_aln = 0
#initiate per-read outfile
out = open('$PRE-per-read.txt','w')

for name in list(names):
    ALN = name in alnlengths
    UN = name in UN_alnscores
    readlen = readlengths[name][0]
    if ALN:
        nctg = len( contigs[name] )
        aln_n = len( alnlengths[name] )
        alnlen = sum( alnlengths[name] )
        edit = sum( editdists[name] )
        alnscore = np.mean( alnscores[name] )
        mapq = np.mean( mapqs[name] )
        unmappedreadlen = 0
        if UN:
            editun = 0
            aln_n2 = aln_n + len( UN_alnscores[name] )
            alnscore2 = np.mean( alnscores[name]+UN_alnscores[name] )
            mapq2 = np.mean( mapqs[name]+UN_mapqs[name] )
        else:
            editun = 0
            aln_n2 = aln_n
            alnscore2 = alnscore
            mapq2 = mapq
    elif UN:
        nctg = 0
        alnlen = 0
        aln_n = 0
        edit = 0
        editun = readlen
        alnscore = 0
        mapq = 0
        aln_n2 = 0
        alnscore2 = 0
        mapq2 = 0
        unmappedreadlen = readlen
    unalnportion = readlen - alnlen
    totaledit = edit + unalnportion
    ##match = readlen - unmappedreadlen - edit ## error b/c it counts unaligned portions as matches
    match = alnlen - edit
    match2 = readlen - totaledit ## same as alnlen-edit --> readlen-totaledit = readlen-(edit+unalnportion) = readlen-(edit+readlean-alnlen) = readlen-edit-readlen+alnlen = alnlen+(readlen-readlen-edit) = alnlen-edit
    out.write( ('\t').join([str(e) for e in [name, nctg, aln_n, readlen, alnlen, unalnportion, edit, totaledit, editun, alnscore, mapq, aln_n2, alnscore2, mapq2, match, match2]]) + '\n' )

    #additional analysis
    summapq_old += sum( mapqs[name]+UN_mapqs[name] )
    sumalnscore_per_aln += sum( alnscores[name]  )
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
    editsumwithunmapped += edit + editun
    alnscoresum += alnscore
    mapqsum2 += mapq2
    alnscoresum2 += alnscore2
    unmappedreadlensum += unmappedreadlen
    matchsum += match
    matchsum2 += match2
    totaleditsum += totaledit
    unalnportionsum += unalnportion
    if ALN and UN:
        numaln += 1
        numaln_un += 1
    elif ALN:
        numaln += 1
    elif UN:
        numunaln += 1
#close per-read outfile
out.close()


#final analysis - open stats out file and write
#totalaln = 0
#totalunaln = 0
numentries = totalaln + totalunaln
numuniqentries = numaln + numunaln
totalnumaln = totalaln
totalnumuniqaln = numaln
##summapq_old = summapq_old
avgmapq_old = summapq_old / float( totalaln )
pctaln_old = totalnumuniqaln / float( numuniqentries )
pctunaln_old = (numuniqentries - totalnumuniqaln) / float( numuniqentries )
alnratio = totalnumaln / float( totalnumuniqaln )

avgalnlen_allaln = float(alnlensum) / totalnumaln
avgalnlen_allalnread = float(alnlensum) / totalnumuniqaln
avgalnlen_allread = float(alnlensum) / numuniqentries

avgmapq = mapqsum/float(numaln)
avgmapq2 = mapqsum2/float(numaln)
##avgmapq_with_unmap = ((mapqsum+mapqsum2)/2.0)/float(numaln+numunaln)
avgmapq_with_unmap = mapqsum/float(numaln+numunaln)
avgalnscore = alnscoresum/float(numaln)
avgalnscore2 = alnscoresum2/float(numaln)
##avgalnscore_with_unmap = ((alnscoresum+alnscoresum2)/2.0)/float(numaln+numunaln)
avgalnscore_with_unmap = alnscoresum/float(numaln+numunaln)
avgalnscore_per_aln = sumalnscore_per_aln / float( totalaln )


avgedit = editsum/float(alnlensum)
#editalnratio = editsumwithunmapped / float(alnlensum)
#editmatchratio = editsumwithunmapped / (float(alnlensum) - editsum)

##matches = readlensum - unmappedreadlensum - editsum ## erroneously counts unaligned portions of aligned reads as matches.... (i.e. readlen - alnlen gives unaln portion of read)
## All below agree -- only need 1 for printing -- e.g. matches and pctmatches
matches = alnlensum - editsum   
pctmatches = float(matches)/readlensum
matches2 = readlensum - totaleditsum
pctmatches2 = float(matches2)/readlensum
pctmatchsum = float(matchsum)/readlensum
pctmatchsum2 = float(matchsum2)/readlensum

pctalnlen = float(alnlensum)/readlensum

#totaleditsum-unmappedreadlensum
#readlensum-unmappedreadlensum
pctedit_allreads = float(totaleditsum)/(readlensum)
pctmatch_allreads = float(matches)/(readlensum)
pctedit_alnreads = float(totaleditsum-unmappedreadlensum)/(readlensum-unmappedreadlensum)
pctmatch_alnreads = float(matches)/(readlensum-unmappedreadlensum)
pctedit_alns = editsum/float(alnlensum)
pctmatch_alns = float(matches)/float(alnlensum)


pct_0_ctg = num_0_ctg/float(num_0_ctg + num_1_ctg + num_multi_ctg)
pct_1_ctg = num_1_ctg/float(num_0_ctg + num_1_ctg + num_multi_ctg)
pct_multi_ctg = num_multi_ctg/float(num_0_ctg + num_1_ctg + num_multi_ctg)

out = open('$PRE-per-read-stats.txt','w')
out.write( 'numentries\t'+str(numentries)+'\n' )
out.write( 'numuniqentries\t'+str(numuniqentries)+'\n' )
out.write( 'totalnumaln\t'+str(totalnumaln)+'\n' )
out.write( 'totalnumuniqaln\t'+str(totalnumuniqaln)+'\n' )
##out.write( 'numaln\t'+str(numaln)+'\n' ) ## redundant with totalnumuniqaln
out.write( 'numaln_with_un_should_be_0\t'+str(numaln_un)+'\n' )
out.write( 'num_unaln\t'+str(numunaln)+'\n' )
out.write( 'num_0_ctg\t'+str(num_0_ctg)+'\n' ) ## should be same as num_unaln
out.write( 'num_1_ctg\t'+str(num_1_ctg)+'\n' ) ## should be same as totalnumuniqaln
out.write( 'num_multi_ctg\t'+str(num_multi_ctg)+'\n' )
out.write( 'readlensum\t'+str(readlensum)+'\n' )
out.write( 'alnlensum\t'+str(alnlensum)+'\n' )
out.write( 'unalnportionsum\t'+str(unalnportionsum)+'\n' )
out.write( 'unmappedreadlensum\t'+str(unmappedreadlensum)+'\n' )
out.write( 'editsum\t'+str(editsum)+'\n' )
##out.write( 'editsumwithunmapped\t'+str(editsumwithunmapped)+'\n' ) ## captures unaligned reads as part of sum, but not unaligned portion of reads
out.write( 'totaleditsum\t'+str(totaleditsum)+'\n' )
out.write( 'sum_matches\t'+str(matches)+'\n' )
##out.write( 'pctmatches\t'+str(pctmatches)+'\n' ) ## redundant with pctmatch_allreads below

## below matchsum/matches2 stuff is redundant with matches (above)
##out.write( 'matchsum\t'+str(matchsum)+'\n' ) 
##out.write( 'matchsum2\t'+str(matchsum2)+'\n' )
##out.write( 'pctmatchsum\t'+str(pctmatchsum)+'\n' )
##out.write( 'pctmatchsum2\t'+str(pctmatchsum2)+'\n' )
##out.write( 'matches2\t'+str(matches2)+'\n' )
##out.write( 'pctmatches2\t'+str(pctmatches2)+'\n' )


out.write( 'pctalnlen\t'+str(pctalnlen)+'\n' )
out.write( 'pctedit_allreads\t'+str(pctedit_allreads)+'\n' )
out.write( 'pctmatch_allreads\t'+str(pctmatch_allreads)+'\n' )
out.write( 'pctedit_alnreads\t'+str(pctedit_alnreads)+'\n' )
out.write( 'pctmatch_alnreads\t'+str(pctmatch_alnreads)+'\n' )
out.write( 'pctedit_alns\t'+str(pctedit_alns)+'\n' )
out.write( 'pctmatch_alns\t'+str(pctmatch_alns)+'\n' )


out.write( 'pct_reads_that_aln\t'+str(pctaln_old)+'\n' )
out.write( 'alnratio\t'+str(alnratio)+'\n' )
out.write( 'pct_0_ctg\t'+str(pct_0_ctg)+'\n' ) ## should be same as pct_reads_that_dont_align
out.write( 'pct_1_ctg\t'+str(pct_1_ctg)+'\n' ) ## should be <= pct_reads_that_aln
out.write( 'pct_multi_ctg\t'+str(pct_multi_ctg)+'\n' ) ## should be <= pct_reads_that_aln

out.write( 'avgalnlen_allaln\t'+str(avgalnlen_allaln)+'\n' )
out.write( 'avgalnlen_allalnread\t'+str(avgalnlen_allalnread)+'\n' )
out.write( 'avgalnlen_allread\t'+str(avgalnlen_allread)+'\n' )

out.write( 'summapq_per_aln\t'+str(summapq_old)+'\n' )
out.write( 'avgmapq_per_aln\t'+str(avgmapq_old)+'\n' )
out.write( 'mapqsum_per_read\t'+str(mapqsum)+'\n' )
out.write( 'avgmapq_per_read\t'+str(avgmapq)+'\n' )
out.write( 'avgmapq_with_unmap\t'+str(avgmapq_with_unmap)+'\n' )


out.write( 'alnscoresum_per_aln\t'+str(sumalnscore_per_aln)+'\n' )
out.write( 'avgalnscore_per_aln\t'+str(avgalnscore_per_aln)+'\n' )
out.write( 'alnscoresum_aln_reads\t'+str(alnscoresum)+'\n' )
out.write( 'avgalnscore_aln_reads\t'+str(avgalnscore)+'\n' )
out.write( 'avgalnscore_allreads\t'+str(avgalnscore_with_unmap)+'\n' )

##out.write( 'avgedit\t'+str(avgedit)+'\n' ) ## redundant with pctedit_alns below
##out.write( 'editalnratio\t'+str(editalnratio)+'\n' ) ## not sure this gives any additional useful info 
##out.write( 'editmatchratio\t'+str(editmatchratio)+'\n' ) ## not sure this gives any additional useful info 

out.write( 'mapqsum2\t'+str(mapqsum2)+'\n' )
out.write( 'avgmapq2\t'+str(avgmapq2)+'\n' )
out.write( 'alnscoresum2\t'+str(alnscoresum2)+'\n' )
out.write( 'avgalnscore2\t'+str(avgalnscore2)+'\n' )


out.write( '\t'+str()+'\n' )
out.close()


###################################

out = open('$PRE-per-read-stats.simple.txt','w')

out.write( 'alnlensum\t'+str(alnlensum)+'\n' )
out.write( 'sum_matches\t'+str(matches)+'\n' )

out.write( 'pctalnlen\t'+str(pctalnlen)+'\n' )
out.write( 'pctmatch_allreads\t'+str(pctmatch_allreads)+'\n' )
out.write( 'pctmatch_alnreads\t'+str(pctmatch_alnreads)+'\n' )
out.write( 'pctmatch_alns\t'+str(pctmatch_alns)+'\n' )

out.write( 'pct_reads_that_aln\t'+str(pctaln_old)+'\n' )
out.write( 'alnratio\t'+str(alnratio)+'\n' )
out.write( 'pct_multi_ctg\t'+str(pct_multi_ctg)+'\n' ) ## should be <= pct_reads_that_aln

out.write( 'avgalnlen_allaln\t'+str(avgalnlen_allaln)+'\n' )
out.write( 'avgalnlen_allalnread\t'+str(avgalnlen_allalnread)+'\n' )
out.write( 'avgalnlen_allread\t'+str(avgalnlen_allread)+'\n' )

out.write( 'avgmapq_per_aln\t'+str(avgmapq_old)+'\n' )
out.write( 'avgmapq_per_read\t'+str(avgmapq)+'\n' )
out.write( 'avgmapq_with_unmap\t'+str(avgmapq_with_unmap)+'\n' )

out.write( 'avgalnscore_per_aln\t'+str(avgalnscore_per_aln)+'\n' )
out.write( 'avgalnscore_aln_reads\t'+str(avgalnscore)+'\n' )
out.write( 'avgalnscore_allreads\t'+str(avgalnscore_with_unmap)+'\n' )

out.close()

" 
