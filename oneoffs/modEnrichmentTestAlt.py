#!/usr/bin/env python2.7
import sys, os, argparse, pybedtools, scipy.stats, string, random
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description="""
#
# Given:
# a set of modification sites
# a set of target regions
# genome file (with chr lengths)

# A target base or dinuc

# Assumes following available:
- BEDtools

# See modEnrichmentTest.py help message for more information.
Basically, that is testing whether target regions have a higher modification rate.
E.g. What is Pr(mod | region)? Or What is the percent of bases modified in the region?
That is compare to Pr( mod | genome ): or percent of bases mod in genome.

For example:
1% of bases are modified in the genome.
2% are modified in the region.

This script asks the for the opposite:
What is Pr(region | mod)? Or What is the percent of all mods labeled as region?
That is compared to Pr( region | genome ): or percent of genome labeled as region.

For example:
40% of the genome is labeled as in repeats.
55% of modifications are labeled as in repeats.

So one is looking at whether modification rate goes up.
The other is looking at whether "Region rate" goes up.
They should be related, since:
Pr( mod | region ) = Pr(region|mod)*Pr(mod)/Pr(region) = Pr(mod,region)/Pr(region)
Pr( region | mod ) = Pr(mod|region)*Pr(region)/Pr(mod) = Pr(mod,region)/Pr(mod)

i.e. Both are the joint prob scaled by something.....





    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--modBED', '-a',
                   type= str,
                   help='''Path to BED file of coordinates for DNA mods.''')
parser.add_argument('--regionBED', '-b',
                   type= str,
                   help='''Path to BED file of target regions. Or to file B in BEDtools sense.''')
##parser.add_argument('--fasta', '-F',
##                   type= str,
##                   help='''Path to fasta file with sequences in it (e.g. genome sequence).''')
##parser.add_argument('--background', '-bf',
##                   type= str, default=False,
##                   help='''Path to faCount file from fasta to skip faCount step .''')
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
parser.add_argument('--noOverlapping', '-N',
                   action='store_true', default=False,
                   help='''Optional. When/If shuffling, use BEDTools' noOverlapping argument.''')

parser.add_argument('--complement', '-C',
                   action='store_true', default=False,
                   help='''Optional. Analyze complement regions as well.''')

##parser.add_argument('--perm', '-P',
##                   type=int, default=False,
##                   help='''Add this number of permutation tests at the end.''')

##parser.add_argument('--exclude', '-e',
##                   type= int, default=False,
##                   help='''Optional. When/If shuffling, prove excluded regions file.''')

parser.add_argument('--altoutput', '-A',
                   action='store_true', default=False,
                   help='''Optional.....''')

parser.add_argument('--table', '-T',
                   type=str, default=False,
                   help='''Give table output with following columns:

NOTE:
User must provide comma-separated string for what to put in the Feature and Modification columns. E.g. Repeats,5mC or Exons,6mA''')


parser.add_argument('--tableheader', '-H',
                   action='store_true', default=False,
                   help='''Optional.....''')

args = parser.parse_args()







def get_random_barcode(k=20):
    return ''.join(random.choice(string.ascii_lowercase+string.ascii_uppercase+string.digits) for x in range(k))



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
  expNumPeaks = expProportion*numPeaksA
  pValue = sum( scipy.stats.binom.pmf(range(obsNumOverlap, numPeaksA+1), numPeaksA, expProportion) )
      #q=obsNumOverlap, size=numPeaksA, prob=expProportion, lower=FALSE) 
  return {'obs':obsProportion, 'exp':expProportion, 'p':pValue, 'expnum':expNumPeaks}



##def readFaCount(name, rm=True):
##    with open(name) as fh:
##        keys = fh.readline().strip().split()[1:]
##        values = fh.readline().strip().split()[1:]
##        counts = {keys[i]:float(values[i]) for i in range(len(keys))}
##        values = fh.readline().strip().split()[1:]
##        prcnts = {keys[i]:float(values[i]) for i in range(len(keys))}
##
##        if rm:
##            os.system('rm '+name)
##            #pass
##        return {'counts':counts, 'prcnts':prcnts}
##
##def faCount(fasta):
##    
##    name = 'tmp-modEnrich.'+get_random_barcode()+'.txt'
##    cmd = 'faCount -strands -dinuc -summary ' + fasta + ' > ' + name
##    success = os.system(cmd)
##    if success == 0:
##        return readFaCount(name)
##    else:
##        assert 'faCount Error'

##
##def faCountTarget(fasta, B):
##    fa = 'tmp-modEnrich.'+get_random_barcode()+'.fasta'
##    B.sequence(fasta, fo=fa)
##    counts = faCount(fa)
##    os.system('rm '+fa)
##    return counts
##       




#1. Use faCount to get base counts in entire sequence file
##sys.stderr.write('1. Bg FaCount \n')
##if args.background:
##    refCounts = readFaCount(args.background, rm=False)
##else:
##    refCounts = faCount(args.fasta)

# 2.	Load in mod BED
sys.stderr.write('2. Load A\n')

A = pybedtools.BedTool(args.modBED) # Don't want to merge modifications here... certainly not across strands....merge().slop(b=0, g=args.genome)

# 3. Load in region BED
sys.stderr.write('3. Load B\n')
B = pybedtools.BedTool(args.regionBED).merge().slop(b=0, g=args.genome)
if args.slop:
    B = B.slop(b=args.slop, g=args.genome).merge()
if args.shuffle:
    B = B.shuffle(g=args.genome, chrom=args.chrom, noOverlapping=args.noOverlapping).sort().merge()

##Bc = B.complement(g=args.genome)

# 4. Use BEDtools to extract target region sequences and faCount to get base counts in target regions
##sys.stderr.write('4. Region FaCount A\n')
##regionCounts = faCountTarget(args.fasta, B)

# 5. Find the mods that overlap regions
sys.stderr.write('5. A in B\n')
modsInRegions = A.intersect(B, u=True)

##modsInComplement = A.intersect(Bc, u=True)

#6. For the given target base, calculate the expected value, observed value, and binomial stats
sys.stderr.write('6. Binomial stats... ')

G = sum([int(line.strip().split()[1]) for line in open(args.genome)]) * 2 ## x2 b/c mods can be on either strand
M = len(A) ## could also be get_sum_lengths(A)
R = get_sum_lengths(B) * 2  ## x2 b/c mods can be on either strand
##Rc = get_sum_lengths(Bc)
M_R = len(modsInRegions)
##M_Rc = len(modsInComplement)
pR = float(R)/G
E_M_R = pR*M
pMR = float(M_R)/M
##pMC = float(M_C)/M

#sys.stderr.write('pValR... ')
#pValueR = sum( scipy.stats.binom.pmf(range(observedInRegions, int(positionsInRegions+1)), positionsInRegions, expectedPropInRegions) )
sys.stderr.write('pValR... ')
pValueR = scipy.stats.binom_test(M_R, M, pR, 'greater')
#sys.stderr.write('pValL... ')
#pValueL = sum( scipy.stats.binom.pmf(range(0,observedInRegions), positionsInRegions, expectedPropInRegions) )
sys.stderr.write('pValL... ')
pValueL = scipy.stats.binom_test(M_R, M, pR, 'less' )
##sys.stderr.write('pVal2... \n')
##pValue2 = scipy.stats.binom_test((M_R, M, pR)


if not args.table:
    print 'G', G
    print 'N mods', M
    print 'Region length', R
    print 'Observed N mods in region', M_R
    print 'pRegion', pR
    print 'Expected N mods in region', E_M_R
    print 'p(Region | mod )', pMR
    print 'pValueR', pValueR
    print 'pValueL', pValueL





if args.table:
    header = ['Feature', 'Modification', 'G', 'SumRegionLength', 'pRegion', 'N_Mods_in_region', 'Expected_N_mods_in_region', 'p(region|mod)', 'Binomial P-value For Enrichment', 'Binomial P-value for Depletion']
    ftr, mod = args.table.strip().split(',')
    out = [ftr, mod, int(G), int(R), pR, int(M_R), int(E_M_R), pMR, pValueR, pValueL]

    if args.tableheader:
        print '\t'.join(header)
    #OUTPUT    
    print '\t'.join([str(e) for e in out])

