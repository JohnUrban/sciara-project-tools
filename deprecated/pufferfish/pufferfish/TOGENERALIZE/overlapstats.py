#!/usr/bin/env python2.7
import sys, argparse, pybedtools, scipy.stats
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description="""

# Given:
# a set of origin intervals (known positives)
# a set of peak intervals (called positives)
# genome file (with chr lengths)
# determine:
# Sensitivity : goal = 1
# FNR : goal = 0
# PPV : goal = 1
# FDR : goal = 0
# SP : goal = 1
# FPR : goal = 0
# NPV : goal = 1
# FOR : goal = 0
# F : goal = 1
# Accuracy: Goal = 1

# a lot more...

    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('--knowns', '-k', '-a',
                   type= str,
                   help='''Path to BED file of coordinates for knowns/groundtruth (e.g. origins). Or to file A in BEDtools sense.''')
parser.add_argument('--peaks', '-p', '-b',
                   type= str,
                   help='''Path to BED file of coordinates called as positives (e.g. peaks). Or to file B in BEDtools sense.''')
parser.add_argument('--genome', '-g',
                   type= str,
                   help='''Path to genome file with chr lengths needed by BEDtools.''')
parser.add_argument('--slop', '-s',
                   type= int, default=False,
                   help='''Optional. When integer provided, peak coordinates are extended that much in both directions before all analyses.''')
parser.add_argument('--shuffle', '-S',
                   action='store_true', default=False,
                   help='''Optional. Shuffle peaks before analysis. This occurs after slopping if that is used. noOverlapping always set to True. Need to set chrom if you want it.''')
parser.add_argument('--chrom', '-c',
                   action='store_true', default=False,
                   help='''Optional. When/If shuffling, keep peaks on same chromosome they start on.''')
##parser.add_argument('--exclude', '-e',
##                   type= int, default=False,
##                   help='''Optional. When/If shuffling, prove excluded regions file.''')

args = parser.parse_args()


def get_sum_lengths(intervals):
    return sum([interval.length for interval in intervals])

def get_sum_overlaps(intervals):
    return sum([int(interval[-1]) for interval in intervals])

def get_percentile_col(intervals, col=-1, q=range(10,100,10)):
    return np.percentile([int(interval[col]) for interval in intervals], q)

def get_last_percentile_with_given_value(intervals, col=-1, q=range(0,101,1), value=0.0):
    pctl = get_percentile_col(intervals, col, q)
    if pctl[0] > value:
        return -1
    for i in range(1,len(q)):
        if pctl[i] > 0:
            return q[i-1]
    return q[i]

def get_mean_col(intervals, col=-1):
    return np.mean([int(interval[col]) for interval in intervals])

def basicBinomialTestForOverlapSignificance(obsNumOverlap, numPeaksA, meanPeakWidthA, numPeaksB, meanPeakWidthB, genomeSize, numConnectedComponents, minOverlap=1):
  totalPositions = genomeSize - (numConnectedComponents*(meanPeakWidthA-1))
  numSuccessfulPositions = (meanPeakWidthA + meanPeakWidthB - minOverlap)*numPeaksB
  numSuccessfulPositions = min(numSuccessfulPositions, totalPositions)
  obsProportion = obsNumOverlap/float(numPeaksA)
  expProportion = numSuccessfulPositions/float(totalPositions)
#   expNumPeaks = expProportion*numPeaksA
  pValue = sum( scipy.stats.binom.pmf(range(obsNumOverlap, numPeaksA+1), numPeaksA, expProportion) )
      #q=obsNumOverlap, size=numPeaksA, prob=expProportion, lower=FALSE) 
  return {'obs':obsProportion, 'exp':expProportion, 'p':pValue}




#0.	num_genome_bp = total genome bp
chrsizes = {}
num_genome_bp = 0
with open(args.genome) as G:
    for chrom in G:
        chrom = chrom.strip().split()
        chrsizes[chrom[0]] = int(chrom[1])
        num_genome_bp += chrsizes[chrom[0]]

## initialize
scorenames = ['SN','FNR','PPV','FDR','SP','FPR', 'NPV', 'FOR', 'F_score','Accuracy', 'Benchmark_Accuracy', 'SN_Exp_Num_Knowns_Ovlp_Peaks', 'SN_Obs_Num_Knowns_Ovlp_Peaks', 'SN_Pvalue_Num_Knowns_Ovlp_Peaks', 'PPV_Exp_Num_Peaks_Ovlp_Knowns', 'PPV_Obs_Num_Peaks_Ovlp_Knowns', 'PPV_Pvalue_Num_Peaks_Ovlp_Knowns', 'Exp_F_score_overlap_method', 'Obs_F_score_overlap_method']
scores = {}
for name in scorenames:
    scores[name] = 0

altnames = {'SN':'P(B|A)','FNR':'P(Bc|A)','PPV':'P(A|B)','FDR':'P(Ac|B)','SP':'P(Bc|Ac)','FPR':'P(B|Ac)', 'NPV':'P(Ac|Bc)', 'FOR':'P(A|Bc)', 'F_score':'2*(SN*PPV)/(SN+PPV)','Accuracy':'P(A,B)+P(Ac,Bc)', 'Benchmark_Accuracy':'NoPositiveTests',  'SN_Exp_Num_Knowns_Ovlp_Peaks':'E(B|A)', 'SN_Obs_Num_Knowns_Ovlp_Peaks':'O(B|A)', 'SN_Pvalue_Num_Knowns_Ovlp_Peaks':'P(O>=O(B|A)|E(B|A))', 'PPV_Exp_Num_Peaks_Ovlp_Knowns':'E(A|B)', 'PPV_Obs_Num_Peaks_Ovlp_Knowns':'O(A|B)', 'PPV_Pvalue_Num_Peaks_Ovlp_Knowns':'P(O>=O(A|B)|E(A|B))',  'Exp_F_score_overlap_method':'2*(SN_exp*PPV_exp)/(SN_exp+PPV_exp)', 'Obs_F_score_overlap_method':'2*(SN_obs*PPV_obs)/(SN_obs+PPV_obs)'}

    
# 1.	Ori_bp = all bp inside origin coordinates
## when finding Ori regions, merge any overlapping origins so as to not count their intersecting lengths twice
oris = pybedtools.BedTool(args.knowns).merge().slop(b=0, g=args.genome)
num_ori_bp = get_sum_lengths(oris)
pA = float(num_ori_bp) / num_genome_bp

    
#2.	Non_Ori_bp = all bp outside of origin coordinates
ori_negative = oris.complement(g=args.genome)
num_nonori_bp = get_sum_lengths(ori_negative)
pAc = float(num_nonori_bp) / num_genome_bp

#3.	Peak_bp = all bp inside peaks
peaks = pybedtools.BedTool(args.peaks).merge().slop(b=0, g=args.genome)
if args.slop:
    peaks = peaks.slop(b=args.slop, g=args.genome).merge()
if args.shuffle:
    peaks = peaks.shuffle(g=args.genome, chrom=args.chrom, noOverlapping=True).sort().merge()
    
num_peak_bp = get_sum_lengths(peaks)
pB = float(num_peak_bp) / num_genome_bp


#4.	Non_Peak_bp = all bp outside of peak coordinates
peak_negative = peaks.complement(g=args.genome)
num_nonpeak_bp = get_sum_lengths(peak_negative)
pBc = float(num_nonpeak_bp) / num_genome_bp

#5.	SN = P(B|A) = Prob(Peak_bp | Ori_bp) = numPeakBpInsideOriCoords / Ori_bp = % Ori_bp labeled as a peak = TruePositiveRate TPR
peak_AND_ori = oris.intersect(peaks, wo=True)
num_bp_labeled_peak_AND_ori = get_sum_overlaps(peak_AND_ori)
SN = float(num_bp_labeled_peak_AND_ori) / num_ori_bp
scores['SN'] = SN

#6.	FNR = P(Bc|A) = Prob(Non_peak_bp | Ori_bp) = 1-SN = numNonPeakBpInsideOriCoords / Ori_bp = % Ori_bp NOT labeled as a peak
FNR = 1 - SN
scores['FNR'] = FNR

#7.	PPV = P(A|B) = Prob(Ori_bp | Peak_bp) = numOriBpInsidePeakCoords / Peak_bp = % Peak_bp labeled as origin
try:
    PPV = float(num_bp_labeled_peak_AND_ori) / num_peak_bp
except ZeroDivisionError:
    PPV = 0
scores['PPV'] = PPV

#8. FDR = P(Ac|B) = Prob(Non_Ori_bp | Peak_bp) = 1 - PPV = numNonOriBpInsidePeakCoords/Peak_bp = % Peak_bp not labeled as origin
FDR = 1 - PPV
scores['FDR'] = FDR

#9.	SP = P(Bc|Ac) = Prob(Non_Peak_bp | Non_Ori_bp) = numNonPeakBpInsideNonOriCoords / Non_Origin_bp = % Non_Ori_bp NOT labeled as Peak_bp
nonpeak_AND_nonorigin = ori_negative.intersect(peak_negative, wo=True)
num_bp_labeled_nonpeak_AND_nonori = get_sum_overlaps(nonpeak_AND_nonorigin)
SP = float(num_bp_labeled_nonpeak_AND_nonori) / num_nonori_bp
scores['SP'] = SP

#10. FPR = P(B|Ac) = Prob(Peak_bp | Non_Ori_bp) = 1-SP = numPeakBpInsideNonOriCoords / Non_Ori_bp = % Non_Ori_bp labeled as a peak
FPR = 1-SP
scores['FPR'] = FPR

#11.	NPV = P(Ac|Bc) = Prob(Non_Ori_bp | Non_Peak_bp) = numNonOriBpInsideNonPeakCoords / Non_Peak_bp = % Non_Peak_bp labeled as Ori_bp
try:
    NPV = float(num_bp_labeled_nonpeak_AND_nonori) / num_nonpeak_bp
except ZeroDivisionError:
    NPV = 0
scores['NPV'] = NPV

#12. FOR = P(A|Bc) = Prob(Ori_bp | Non_Peak_bp) = numOriBpInsideNonPeakCoords / Non_Peak_bp = % Non_Peak_bp labeled as Ori_bp
FOR =  1-NPV
scores['FOR'] = FOR

#13. Calculate F_score = 2.0*(SN*PPV)/(SN+PPV) = 2.0 * (P(B|A) * P(A|B)) / ((P(B|A) + P(A|B))
try:
    F = 2.0*(SN*PPV)/(SN+PPV)
except ZeroDivisionError:
    F = 0
scores['F_score'] = F

#14. Calculate accuracy = (TP+TN)/Total
## Accuracy = (TP + TN)/(TP+TN+FP+FN) = (Number correct)/Number all
## TP = True Pos = num_bp_labeled_peak_AND_ori
## TN = true neg = num_bp_labeled_nonpeak_AND_nonori
## If prevalence is known -- a priori probability of known ... num_ori_bp/num_genome_bp
## Accuracy = (sensitivity) (prevalence) + (specificity) (1 - prevalence) = P(B|A)*P(A) + P(Bc|Ac)*P(Ac) = P(A,B)+P(Ac,Bc)
## Total = num_genome_bp
## print SN*(float(num_ori_bp)/num_genome_bp) + SP*(1 - float(num_ori_bp)/num_genome_bp)
try:
    Accuracy = (num_bp_labeled_peak_AND_ori + num_bp_labeled_nonpeak_AND_nonori) / float(num_genome_bp)
except ZeroDivisionError:
    Accuracy = 0
scores['Accuracy'] = Accuracy

#15. Benchmark accuracy is the accuracy we would get if we just labeled all bases as nonorigin
# "nonorigin' is chosen here b/c most of the genome is non-origin -- so one can be pretty confident any
# base chosen at random is non-origin
## This means it is 1 - num_ori_bp/Total
Bench_Acc = 1 - float(num_ori_bp)/num_genome_bp
scores['Benchmark_Accuracy'] = Bench_Acc

#16. calculate exp and obs number of peaks that overlap origins and vice versa - and other fixed params
numPeaks = len(peaks)
try:
    meanPeakWidth = num_peak_bp / float(numPeaks)
except ZeroDivisionError:
    meanPeakWidth = 0
numOris = len(oris)
meanOriWidth = num_ori_bp / float(numOris)
peaksOverlapOris = peaks.intersect(oris, u=True)
numPeaksOverlapOris = len(peaksOverlapOris)
orisOverlapPeaks = oris.intersect(peaks, u=True)
numOrisOverlapPeaks = len(orisOverlapPeaks)
genomeSize = num_genome_bp
numConnectedComponents = len(chrsizes.keys())

#17. calculate bionomial p-value for number of peaks overlapping origins
obsNumOverlap = numOrisOverlapPeaks
numPeaksA = numOris
meanPeakWidthA = meanOriWidth
numPeaksB = numPeaks
meanPeakWidthB = meanPeakWidth
ori_ovlp_peak = basicBinomialTestForOverlapSignificance(obsNumOverlap, numPeaksA, meanPeakWidthA, numPeaksB, meanPeakWidthB, genomeSize, numConnectedComponents, minOverlap=1)
scores['SN_Exp_Num_Knowns_Ovlp_Peaks'] = ori_ovlp_peak['exp']
scores['SN_Obs_Num_Knowns_Ovlp_Peaks'] = ori_ovlp_peak['obs']
scores['SN_Pvalue_Num_Knowns_Ovlp_Peaks'] = ori_ovlp_peak['p']



#18. calculate bionomial p-value for number of peaks overlapping origins
obsNumOverlap = numPeaksOverlapOris
numPeaksA = numPeaks
meanPeakWidthA = meanPeakWidth
numPeaksB = numOris
meanPeakWidthB = meanOriWidth
peak_ovlp_ori = basicBinomialTestForOverlapSignificance(obsNumOverlap, numPeaksA, meanPeakWidthA, numPeaksB, meanPeakWidthB, genomeSize, numConnectedComponents, minOverlap=1)
scores['PPV_Exp_Num_Peaks_Ovlp_Knowns'] = peak_ovlp_ori['exp']
scores['PPV_Obs_Num_Peaks_Ovlp_Knowns'] = peak_ovlp_ori['obs']
scores['PPV_Pvalue_Num_Peaks_Ovlp_Knowns'] = peak_ovlp_ori['p']

#19. F_score of the overlap method -- = 2.0*(SN*PPV)/(SN+PPV) = 2.0 * (P(B|A) * P(A|B)) / ((P(B|A) + P(A|B))
try:
    F2_E = 2.0 * (ori_ovlp_peak['exp'] * peak_ovlp_ori['exp']) / (ori_ovlp_peak['exp'] + peak_ovlp_ori['exp'])
except ZeroDivisionError:
    F2_E = 0
scores['Exp_F_score_overlap_method'] = F2_E
try:
    F2_O = 2.0 * (ori_ovlp_peak['obs'] * peak_ovlp_ori['obs']) / (ori_ovlp_peak['obs'] + peak_ovlp_ori['obs'])
except ZeroDivisionError:
    F2_O = 0
scores['Obs_F_score_overlap_method'] = F2_O

#20. 

## one could also derive the expected numbers for all the other things in this analysis as well...
## e.g. expNumSucPositions is actually the exp conditionals -- i.e. expected conditional A|B -- or B|A depending on which way analysis is going..
## e.g. since expNumSucPositions assumes independence, pB = B|A  and pA = A|B
## e.g. since conditionals and priors are known, can find joint  A,B = p(B|A)p(A) = p(A|B)p(B)) = p(A)p(B)
## e.g. since joint and priors are known, can get exp union = p(A)+p(B)-p(A)p(B)
## e.g. since priors are known, can get complements: Ac = 1-A  and Bc = 1-B
## e.g. since conditionals are known, can get expected complements = Ac|B = 1-A|B; Bc|A = 1-B|A
## e.g. since Ac|B and Bc|A are known as well as priors, can get: B|Ac = p(Ac|B)p(B)/p(Ac); and A|Bc = p(Bc|A)p(A)/p(Bc)
##   --> since independence assumed B|Ac = p(Ac)p(B)/p(Ac) = p(B)
##   --> since independence assumed A|Bc = p(Bc)p(A)/p(Bc) = p(A)
## e.g. since independence assumed more generally -- all p(X|Y) equates to p(X) for expectation
## OBSERVED  -- I will have to figure this out later...
## e.g. we can observe the conditionals A|B and B|A -- which I take by proportion of peaks... but it could be... 
## e.g.    observed joint (intersect) could be..... taking peaks from both sets that overlap, merging, and summing the length...
## e.g. Union is length of all merged....
## e.g. expected ~A OR ~B = 1-expNumSucPositions = A+B-(A,B)
## SN = P(B|A) already used -- but this means we can give FNR=1-SN = P(Bc|A)
## PPV = P(A|B) already used -- but this means we cn give FDR=1-PPV = P(Ac|B)
## That leaves P(A|Bc), P(B|Ac), P(Ac|Bc) and P(Bc|Ac)
## P(Bc|Ac) = P(Ac|Bc)P(Bc)/P(Ac) = P(Ac,Bc)/P(Ac)
## P(B|A) = P(A|B)P(B)/P(A)
## P(A) = P(B|A)/(P(A|B)P(B)) = P(B|A)/P(A,B)
## problem here is... p(A)=p(A|A) only if p(A) independent of p(A)... we know they are not b/c p(A|A) = 1 != p(A)
## --------> A,B independent if P(A|B) = P(A), P(B|A) = P(B), and P(A,B) = P(A)P(B)
## if A and B are knowable, so are observed p(Ac) and p(Bc)
## Do I have to assume p(A) is same for this as it was for the labeling/single_bp approach? i.e. sum_A_lengths/genome
##
## if p(A) known -- then p(Ac) known -- same for B
## P(Ac,Bc) = 1-AUB
## AUB = p(A)+p(B)-p(A,B)
## p(A,B) = p(A|B)p(B) = p(B|A)p(A)




## JACCARD
orijaccard = oris.jaccard(peaks)['jaccard']
peaksjaccard = peaks.jaccard(oris)['jaccard']

orisfisher = oris.fisher(peaks, g=args.genome).right_tail
peaksfisher = peaks.fisher(oris, g=args.genome).right_tail

# Mean distances
# Mean distance from Ori to closest peak -- In general I thought this was too influenced by outliers
##closestori = get_mean_col(peaks.closest(oris,d=True))
# Mean distance from Ori to closest peak -- In general I thought this was too influenced by outliers
##closestpeak = get_mean_col(oris.closest(peaks,d=True))

# Various percentiles of distances - found that most were distances of 0 assuming decent %overlap
##print get_percentile_col(peaks.closest(oris,d=True))
##print get_percentile_col(oris.closest(peaks,d=True))

## So I tried getting the last percentile where 0 was the distances to closest....
## but this just essentially gives lower res pct overlap...
## i.e. the actual percentile that is the last one with 0 IS THE PCT OVERLAP from binomial analysis
##print get_last_percentile_with_given_value(peaks.closest(oris,d=True))
##print get_last_percentile_with_given_value(oris.closest(peaks,d=True))
##print get_last_percentile_with_given_value(peaks.closest(oris,d=True),value=100)
##print get_last_percentile_with_given_value(oris.closest(peaks,d=True),value=1000)


## OUTPUT
##TODO
# P(A|B) = P(A,B)/P(B); P(B|A) = P(A,B)/P(A)
## A,B independent if P(A|B) = P(A), P(B|A) = P(B), and P(A,B) = P(A)P(B)
print "Prevalence\tP(A)\t" + str(pA) + "\tExp_PPV_and_Exp_FOR"
print "PositiveRate\tP(B)\t" + str(pB) + "\tExp_SN_and_Exp_FPR"
print "Product\tP(A)P(B)\t" + str(pA*pB) + "\tExp_joint"
print "Conditional_P(A|B)\tP(B|A)P(A)/P(B)\t" + str(SN*pA/pB)
print "Conditional_P(B|A)\tP(A|B)P(B)/P(A)\t" + str(PPV*pB/pA)
print "Joint_P(A,B)\tP(B|A)P(A)\t" + str(SN*pA)
print "Joint_P(A,B)\tP(A|B)P(B)\t" + str(PPV*pB)
print "Union\tP(A)+P(B)-P(A,B)\t"+ str(pA+pB-SN*pA)
print "Exp_Union\tP(A)+P(B)-P(A)p(B)\t"+ str(pA+pB-(pA*pB))
print "ExpJaccard\tP(A)P(B)/(P(A)+P(B)-P(A)P(B))\t" + str((pA*pB/(pA+pB-(pA*pB))))
print "ObsJaccard\tP(A,B)/P(AUB)\t" + str((SN*pA) / (pA+pB-SN*pA))
print "Exp_FDR\tP(Ac)\t" + str(pAc) + "\t" + str("Exp_NPV")
print "Exp_SP\tP(Bc)\t" + str(pBc) + "\t" + str("Exp_FNR")
print "Exp_F_score\t2*(Exp_SN*Exp_PPV)/(Exp_SN+Exp_PPV)\t" + str(2.0*(pA*pB)/(pA+pB))
## Accuracy = (sensitivity) (prevalence) + (specificity) (1 - prevalence) = P(B|A)*P(A) + P(Bc|Ac)*P(Ac) = P(A,B)+P(Ac,Bc) = pA*pB+pAc*pBc
print "Exp_Accuracy\tp(A)*p(B)+p(Ac)*p(Bc)\t" + str(pA*pB+pAc*pBc)
outmsg = ("\n").join([name+"\t"+altnames[name]+"\t"+str(scores[name]) for name in scorenames])
print outmsg



##print "P(Ac) P(Bc) P(A,Bc) ... etc...."

print "BEDtoolsKnownJaccard\tIntersect/Union\t" + str(orijaccard)
print "BEDtoolsPeakJaccard\tIntersect/Union\t" + str(peaksjaccard)
print "ProbJaccard\tP(A,b)/P(AUB)\t" + str((SN*pA) / (pA+pB-SN*pA))
print "Fisher_AB\tExactTestPvalue\t" + str(orisfisher)
print "Fisher_BA\tExactTestPvalue\t" + str(peaksfisher)
print "ROC\tx,y\t" + str(FPR)+","+str(SN)


#NOTE: ROC curves are x vs y = FPR vs Sensitivity (TPR)
## see: http://www.lexjansen.com/nesug/nesug10/hl/hl07.pdf


## TODO:
## Given optional bedGraph of signal from B, determine if signal more enriched in A coordinates than:
## 1. B coordinate avg
## 2. genomic bp avg
## 3. shuffle avg....
## Same for signal from A -- if pertinent....
## Can do parametric and non-parametric tests....
## Also find out if true positives have stronger enrichments than false positives...
## Find out if false negatives have higher enrichment than true negatives...
## Find out if true-positives are enriched in top ranks of peaks sorted by p, q, FE, etc...
## Can also get union set of coordinates to ensure all known regions are in ranked list....
