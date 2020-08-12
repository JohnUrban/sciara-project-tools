#!/usr/bin/env python2.7
import sys, argparse, pybedtools, scipy.stats
from collections import defaultdict
import numpy as np

import rpy2.robjects as robjects
ksmooth = robjects.r['ksmooth']
intvec = robjects.IntVector
fltvec = robjects.FloatVector


parser = argparse.ArgumentParser(description="""

 Given gc_fe_table.txt
 Calculate median FE for each GC content 0,1,2,......,100

 Typical procedure beforehand:
 bowtie2 -p $P --very-sensitive -N 1 -x hg19 -U treat.fastq | samtools view -bSh -F 4 -q 0 | samtools sort --threads $P -T TREAT | samtools rmdup -s - treat.bam
 genomeCoverageBed -bg -i treat.bam -g hg19.genome (-scale factor) > treat.bedGraph

 bowtie2 -p $P --very-sensitive -N 1 -x hg19 -U control.fastq | samtools view -bSh -F 4 -q 0 | samtools sort --threads $P -T TREAT | samtools rmdup -s - control.bam
 genomeCoverageBed -bg -i control.bam -g hg19.genome (-scale factor) > control.bedGraph

 bedtools makewindows -g hg19.genome -w 100 -s 100 | bedtools intersect -v -a - -b excludedregions.bed | awk '$3-$2 == 100 {OFS="\t"; print $1,$2,$3}' | bedtools nuc -fi hg19.fa -bed - | awk '!/^#/ {OFS="\t"; print $1,$2,$3,$5}' > gc.bedGraph
 bedtools unionbedg -i gc.bedGraph treat.bedGraph control.bedGraph | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$5/$6}' > gc_fe_table.txt


 NOTE:
 Pipeline for correcting signal that also has CNVs:
     1. GC median
     2. CNV correction given GC median
     3. GC median update
     4. Subtraction or Z-score correction
     5. Peak calling via height and length cutoffs or BCP analysis. Optional p-value/q-value tracks, peak-callng, etc.
    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('--table', '-i', '-f',
                   type= str,
                   help='''Path to gc_fe_table.txt.''')


parser.add_argument('--sigcol',
                   type=int, default=7,
                   help='''Assuming gc_fe_table.txt, default is 7. Change if necessary.''')

parser.add_argument('--gccol', 
                   type=int, default=4,
                   help='''Assuming gc_fe_table.txt, default is 7. Change if necessary.''')

parser.add_argument('--scale', 
                   type=float, default=100,
                   help='''Multiplies GC column by this number. Default = 100. Assumes GC are proportions output by BEDtools.
Also automatically rounds to nearest integer.''')


parser.add_argument('--mean', '-M', 
                   action='store_true', default=False,
                   help='''Subtract GC mean from signal instead of GC median.''')

parser.add_argument('--chrcol', '-c', 
                   type=int, default=1,
                   help='''Column chr name found in. Default = 1''')
parser.add_argument('--startcol', '-S', 
                   type=int, default=2,
                   help='''Column start coordinate found in. Default = 2''')
parser.add_argument('--endcol', '-E', 
                   type=int, default=3,
                   help='''Column end coordinate found in. Default = 3''')

parser.add_argument('--sigcoltocorrect',
                   type=int, default=False,
                   help='''Assuming gc_fe_table.txt, default is to use same column in --sigcol -- i.e. correct the signal used to get parameters. Change if necessary.
                        i.e. if want to use the params from signal in one col to correct signal in another.''')

parser.add_argument('--control',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the GC_median over each bin given the bin's GC content.
This can be used (1) in a browser, (2) as part of command-line subtraction, fold-enrcichment and other operations w/ paste/awk and a signal file,
(3) down stream significance tests outside of this program, etc etc.''')



parser.add_argument('--mad',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the median absolute deviation from each GC_median over each bin given the bin's GC content.
This can be used (1) in a browser, (2) as part of command-line operations w/ paste/awk and a signal and/or median file,
(3) down stream significance tests outside of this program, etc etc.''')


parser.add_argument('--zscore',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the zscore over each bin.
The z-score is obtained by (bin_value - gc_mean)/gc_std_dev given the bin's GC content.''')

parser.add_argument('--robust_zscore',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the zscore over each bin.
The z-score is obtained by (bin_value - gc_median)/(c*gc_mad) given the bin's GC content.
When MedianAD is 0, c*gc_mad = 1.253314*MeanAD
When MedianAD is not 0, c*gc_mad = 1.4826*MedianAD
According to Wikipedia -- this should use 1.4826 
https://en.wikipedia.org/wiki/Median_absolute_deviation
According to IBM - 1.486:
See: https://www.ibm.com/support/knowledgecenter/en/SSWLVY_1.0.0/com.ibm.spss.analyticcatalyst.help/analytic_catalyst/modified_z.html
or 1.483 according to another source (which is 1.4826 rounded)...
or 1.4296 according to another (1.43 rounded)...
''')

parser.add_argument('--madunits',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the zscore-like madunits over each bin.
The z-score-like madunits are obtained by (bin_value - gc_median)/gc_mad given the bin's GC content.
When MAD is 0, it attempts to replace MAD with the MeanAD from the median.
When both MAD and MeanAD are 0, it attempts to replace MAD with 0.67449*std_dev.
When all are 0 (should not happen) - it replaces MAD with 1.
This also differs from the robust_zscore because it does not scale the MAD.
The robust_zscore does scale MAD:  (x-MED)/scaled_MAD.
''')


parser.add_argument('--medfe',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the fold enrichment of bin/med_gc over each bin given the bin's GC content.
This can be used with smoothing to get long-range mean relative copy numbers.
Long-range RCNs over each bin can be used to scale the bin's value (bin/bin_RCN) before subtracting GC_med or getting z-score in subsequent steps.
When median is 0, it is changed to 1 for the denominator.
This may or may not give the desired effect.''')

parser.add_argument('--meanfe',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the fold enrichment of bin/mean_gc over each bin given the bin's GC content.
This can be used with smoothing to get long-range mean relative copy numbers.
Long-range RCNs over each bin can be used to scale the bin's value (bin/bin_RCN) before subtracting GC_med or getting z-score in subsequent steps.
When mean is 0, it is changed to 1 for the denominator.
This may or may not give the desired effect.''')


parser.add_argument('--subtractmed',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the median-subtracted scores of bin-med_gc over each bin given the bin's GC content.''')

parser.add_argument('--subtractmean',
                   type= str, default=False,
                   help='''Give an output prefix to also output a bedGraph that contains the mean-subtracted scores of bin-med_gc over each bin given the bin's GC content.''')


parser.add_argument('--dist',
                   type= str, default=False,
                   help='''Give an output prefix to also output a 2-column table that includes GC content value and all values encountered (used to calculate median, mean, stdev, mad).

This can be used (1) to re-create the stats dict made here w/o seeing the original file again,
(2) to plot distributions of values seen for each GC...., etc...''')

parser.add_argument('--bdg', '-b',
                   type= str, default=False,
                   help='''OUTDATED OPTION. THERE ARE MANY POSSIBLE BEDGRAPHS TO OUTPUT NOW. MOREOVER THIS IS ACTUALLY JUST A SUBTRACTION BDG.
Give an output prefix to also output a bedGraph of corrected signal.
(ONLY KEPT OPTION TO PREVENT BREAKAGE OF OLDER PIPELINES)...''')


args = parser.parse_args()

def name_bdg(bdg):
    if not bdg.endswith('.bedGraph'):
        return bdg + '.bedGraph'
    return bdg


## DECLARE SOME VALUES
gccol = args.gccol-1
sigcol = args.sigcol-1
chrcol = args.chrcol-1
startcol = args.startcol-1
endcol = args.endcol-1
if not args.sigcoltocorrect:
    args.sigcoltocorrect = args.sigcol
sigcoltocorrect = args.sigcoltocorrect-1

gcsub = 0
if args.mean:
    gcsub = 1

## COLLECT LIST OF SIGNAL VALUES FOR EACH GC VALUE
## initialize
gc2felist = defaultdict(list)
for i in range(0,101,1):
    gc2felist[i]
with open(args.table) as table:
    for row in table:
        row = row.strip().split()
        gc2felist[int(100.0*float(row[gccol]))].append( float(row[sigcol]) )




        
## COMPUTE SIGNAL STATS FOR EACH GC VALUE
## IF OPTED FOR, OUTPUT THE SIGNAL DISTRIBUTIONS FOR EACH GC VALUE
if args.dist:
    distout = open(args.dist+'.gc_signal_distributions.txt', 'w')
# initialize dictionary
statdict = defaultdict(list)
mad_sd_ratio = []
for i in range(0,101,1):
    n = len(gc2felist[i])
    if n > 0:
        median = np.median(gc2felist[i])
        mean =  np.mean(gc2felist[i])
        binsum = np.sum(gc2felist[i])
        if n == 1: ## stdev will be NAN and MAD will be 0 
            ## This is just a band-aid for now
            ## Moreover, it is not guaranteed to work 
            try:
                vals = np.array(gc2felist[i-1] + gc2felist[i] + gc2felist[i+1])
                std = np.std(vals, ddof=1)
                mad = np.median( np.absolute(vals - median) )
                mad2 = np.mean( np.absolute(vals - median) )
            except: ## when i=0 or i=N, this will crash.
                std = 1
                mad = 1
                mad2 = 1
        else:
            std = np.std(gc2felist[i], ddof=1) ## ddof=1 gives same result as R
            mad = np.median( np.absolute(gc2felist[i] - median) ) ## median absolute deviation from median
            mad2 = np.mean( np.absolute(gc2felist[i] - median) ) ## mean absolute deviation from median

        statdict[i] += [median, mean, std, mad, mad2, n, binsum]
        if args.dist:
            out = [str(i), (',').join([str(e) for e in gc2felist[i]])]
            outmsg = ('\t').join(out)
            distout.write(outmsg + '\n')

if args.dist:
    distout.close()



## RETURN STATS TABLE
## PRINTED TO STDOUT FOR NOW -- THIS WAS AT FIRST THE MAIN OUTPUT OF THIS PROGRAM
## OUTPUT SIGNAL STATS FOR EACH GC VALUE
for i in range(0,101,1):
    if len(gc2felist[i]) != 0:
        out = [i] + statdict[i]
        print ("\t").join([str(e) for e in out])




## 2018-01-04 - I would like to re-visit if having bedGraph generation in this script is desirable.
## Seems to me that (1) it needs to go back through the table file anyway, and (2) loading the stat information back into python in a different script is not costly...
## Keeping it below is starting to make a mess - especially as I change how I'd like to proceed.

## WRITE EXTRA BEDGRAPHS IF DESIRED
if args.bdg or args.control or args.mad or args.zscore or args.robust_zscore or args.madunits or args.medfe or args.meanfe or args.subtractmed or args.subtractmean:
    if args.bdg: ## WRITE CORRECTED SIGNAL FILE IF DESIRED
        sigbdg = open(name_bdg(args.bdg),'w')
    if args.control: ## WRITE MEDIAN CONTROL BEDGRAPH IF DESIRED (mean printed out if --mean specified)
        conbdg = open(name_bdg(args.control),'w')
    if args.mad: ## WRITE MAD BEDGRAPH IF DESIRED
        madbdg = open(name_bdg(args.mad),'w')
    if args.zscore:
        zbdg = open(name_bdg(args.zscore), 'w')
    if args.robust_zscore:
        rzbdg = open(name_bdg(args.robust_zscore), 'w')
    if args.madunits:
        madunitbdg = open(name_bdg(args.madunits), 'w')
    if args.medfe:
        febdg = open(name_bdg(args.medfe), 'w')
    if args.meanfe:
        mufebdg = open(name_bdg(args.meanfe), 'w')
    if args.subtractmed:
        submedbdg = open(name_bdg(args.subtractmed), 'w')
    if args.subtractmean:
        submubdg = open(name_bdg(args.subtractmean), 'w')
    with open(args.table) as table:
        for row in table:
            row = row.strip().split()
            gc = int(100.0*float(row[gccol]))
            if args.bdg:
                sig = float(row[sigcoltocorrect])
                newsig = sig - statdict[gc][gcsub]
                out = [row[chrcol], row[startcol], row[endcol], newsig]
                outmsg = ("\t").join([str(e) for e in out])
                sigbdg.write( outmsg + '\n' )
            if args.control:
                curr_gc_med = statdict[gc][gcsub]
                out = [row[chrcol], row[startcol], row[endcol], curr_gc_med]
                outmsg = ("\t").join([str(e) for e in out])
                conbdg.write( outmsg + '\n' )
            if args.mad:
                curr_gc_mad = statdict[gc][3]
                out = [row[chrcol], row[startcol], row[endcol], curr_gc_mad]
                outmsg = ("\t").join([str(e) for e in out])
                madbdg.write( outmsg + '\n' )
            if args.zscore: #(bin - median)
                sig = float(row[sigcoltocorrect])
                mu = statdict[gc][1]
                std = statdict[gc][2]
                zscore = (sig - mu) / std
                out = [row[chrcol], row[startcol], row[endcol], zscore]
                outmsg = ("\t").join([str(e) for e in out])
                zbdg.write( outmsg + '\n' )
            if args.robust_zscore: #(bin - median)
                sig = float(row[sigcoltocorrect])
                med = statdict[gc][0]
                mad = statdict[gc][3]
                mad2 = statdict[gc][4]
                if mad == 0:
                    denom = 1.253314*mad2
                else:
                    denom = 1.4826*mad
                zscore = (sig - med) / denom
                out = [row[chrcol], row[startcol], row[endcol], zscore]
                outmsg = ("\t").join([str(e) for e in out])
                rzbdg.write( outmsg + '\n' )
            if args.madunits: #(bin - median)
                sig = float(row[sigcoltocorrect])
                med = statdict[gc][0]
                mad = statdict[gc][3]
                mad2 = statdict[gc][4]
                std = statdict[gc][2]
                if mad != 0:
                    denom = mad
                elif mad == 0 and mad2 != 0:
                    denom = mad2
                elif mad == 0 and mad2 == 0 and std > 0:
                    denom = 0.67449 * std
                else:
                    denom = 1
                zscore = (sig - med) / denom
                out = [row[chrcol], row[startcol], row[endcol], zscore]
                outmsg = ("\t").join([str(e) for e in out])
                madunitbdg.write( outmsg + '\n' )
            if args.medfe:
                sig = float(row[sigcoltocorrect])
                med = statdict[gc][0]
                if med == 0:
                    med = 1
                fe = sig / med
                out = [row[chrcol], row[startcol], row[endcol], fe]
                outmsg = ("\t").join([str(e) for e in out])
                febdg.write( outmsg + '\n' )
            if args.meanfe:
                sig = float(row[sigcoltocorrect])
                mean = statdict[gc][1]
                if mean == 0:
                    mean = 1
                fe = sig / mean
                out = [row[chrcol], row[startcol], row[endcol], fe]
                outmsg = ("\t").join([str(e) for e in out])
                mufebdg.write( outmsg + '\n' )
            if args.subtractmed:
                sig = float(row[sigcoltocorrect])
                newsig = sig - statdict[gc][0]
                out = [row[chrcol], row[startcol], row[endcol], newsig]
                outmsg = ("\t").join([str(e) for e in out])
                submedbdg.write( outmsg + '\n' )
            if args.subtractmean:
                sig = float(row[sigcoltocorrect])
                newsig = sig - statdict[gc][1]
                out = [row[chrcol], row[startcol], row[endcol], newsig]
                outmsg = ("\t").join([str(e) for e in out])
                submubdg.write( outmsg + '\n' )
    if args.bdg:
        sigbdg.close()
    if args.control:
        conbdg.close()
    if args.mad:
        madbdg.close()
    if args.zscore:
        zbdg.close()
    if args.robust_zscore:
        rzbdg.close()
    if args.madunits:
        madunitbdg.close()
    if args.medfe:
        febdg.close()
    if args.meanfe:
        mufebdg.close()
    if args.subtractmed:
        submedbdg.close()
    if args.subtractmean:
        submeanbdg.close()

## ANOTHER WAY TO CORRECT: might be to do Signal * Overall_median/GC_median
