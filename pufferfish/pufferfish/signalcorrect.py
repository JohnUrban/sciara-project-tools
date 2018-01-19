#!/usr/bin/env python2.7
import sys, argparse, pybedtools, scipy.stats
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description="""

 Given gcmedian.txt and signal file,
 Generate updated signal file corrected for gcfemedian...

 NOTE: if gcfmedian obtained on FE can only use on FE signal.
 If gcmedian done on SPMR, can only use on SPMR signal.
 Etc.

 NOTE2: signal bins being corrected need to be same bins GC analyzed in.
    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--signal', '-a', '-i', '-f',
                   type= str,
                   help='''Path to signal bedGraph that has cols: chr, start, end, GCprop, signal-to-correct.Can tell algo what cols to look for.''')
parser.add_argument('--scalefactors', '-b', 
                   type=str, 
                   help='''Path to gc median parameter table.''')

parser.add_argument('--gccol', '-g', 
                   type=int, default=6,
                   help='''Column GC proportions found in. Default = 6''')
parser.add_argument('--signalcol', '-s', 
                   type=int, default=4,
                   help='''Column signal found in. Default = 4''')
parser.add_argument('--chrcol', '-c', 
                   type=int, default=1,
                   help='''Column chr name found in. Default = 1''')
parser.add_argument('--startcol', '-S', 
                   type=int, default=2,
                   help='''Column start coordinate found in. Default = 2''')
parser.add_argument('--endcol', '-E', 
                   type=int, default=3,
                   help='''Column end coordinate found in. Default = 3''')
parser.add_argument('--mean', '-M', 
                   action='store_true', default=False,
                   help='''Subtract GC mean from signal instead of GC median.''')

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


args = parser.parse_args()
## could also use NS-seq medians to correct its own signal...
## could do that all in the first step -- either to treatment signal, FE signal, or both...
## if do raw signa or raw spmr -- will still need to normalize for copy/number -- and/or do local lambda peak calling

def name_bdg(bdg):
    if not bdg.endswith('.bedGraph'):
        return bdg + '.bedGraph'
    return bdg

gccol = args.gccol-1
sigcoltocorrect = args.signalcol-1
chrcol = args.chrcol-1
startcol = args.startcol-1
endcol = args.endcol-1

gcsub = 0
if args.mean:
    gcsub = 1


## READ IN GC STATS INFORMATION
statdict = defaultdict(list)
with open(args.scalefactors) as table:
    for row in table:
        row = row.strip().split()
        statdict[int(row[0])] = [float(e) for e in row]
        



## WRITE BEDGRAPHS

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
    
with open(args.signal) as table:
    for row in table:
        row = row.strip().split()
        gc = int(100.0*float(row[gccol]))
##        if args.bdg:
##            sig = float(row[sigcoltocorrect])
##            newsig = sig - statdict[gc][gcsub]
##            out = [row[chrcol], row[startcol], row[endcol], newsig]
##            outmsg = ("\t").join([str(e) for e in out])
##            sigbdg.write( outmsg + '\n' )
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


### THIS WILL BE DELETED EVENTUALLY
## IT IS HERE FOR POSTERITY SO OLDER PIPELINES DO NOT GO BELLY UP
if not (args.control or args.mad or args.zscore or args.robust_zscore or args.madunits or args.medfe or args.meanfe or args.subtractmed or args.subtractmean):
    ## If bedGraph type not specified... then it must be an older pipeline OR user wants subtracmed bedGraph
    with open(args.signal) as table:
        for row in table:
            row = row.strip().split()
            gc = int(100.0*float(row[gccol]))
            sig = float(row[sigcoltocorrect])
            newsig = sig - statdict[gc][gcsub]
            out = [row[chrcol], row[startcol], row[endcol], newsig]
            print ("\t").join([str(e) for e in out])
            

