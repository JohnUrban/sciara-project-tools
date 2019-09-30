#!/usr/bin/env python2.7
import sys, os, argparse, pybedtools, scipy.stats
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description="""

# Given:
# a set of modification sites
# a set of target regions
# genome file (with chr lengths)

# A target base or dinuc

# Assumes following available:
- BEDtools
- faCount from Kent Utilities (UCSC)



# There could simply be more Cs in regions that explain the enrichment of modCs there...
# So need to count the number of Cs in the Genome on both strands
#  and the number of modified Cs in the entire genome. .. expected proportion = nmod/N
# Then the number of Cs in target regions on both strands


The perm tests are such:
1. Shuffle the regions
2. Get observed proportion of target bases that are modified in region
3. Repeat a bunch
4. Count how many times the randomly observed proportions are > then the real observation and divide by the number of tests.

NOTE: The permutation tests take a lot longer, but are more robust.
They are necessary as the number of features in A and B become large as th

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--targetbase', '-B',
                   type= str, required=True,
                   help='''Provide either: A, C, G, T.''')
parser.add_argument('--modBED', '-a',
                   type= str,
                   help='''Path to BED file of coordinates for DNA mods.''')
parser.add_argument('--regionBED', '-b',
                   type= str,
                   help='''Path to BED file of target regions. Or to file B in BEDtools sense.''')
parser.add_argument('--fasta', '-F',
                   type= str,
                   help='''Path to fasta file with sequences in it (e.g. genome sequence).''')
parser.add_argument('--background', '-bf',
                   type= str, default=False,
                   help='''Path to faCount file from fasta to skip faCount step .''')
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

parser.add_argument('--perm', '-P',
                   type=int, default=False,
                   help='''Add this number of permutation tests at the end.''')

##parser.add_argument('--exclude', '-e',
##                   type= int, default=False,
##                   help='''Optional. When/If shuffling, prove excluded regions file.''')

parser.add_argument('--altoutput', '-A',
                   action='store_true', default=False,
                   help='''Optional.....''')

parser.add_argument('--table', '-T',
                   type=str, default=False,
                   help='''Give table output with following columns:
Feature	Modification	Total Target Base	Total Observed Target Base Modifications	Expected Proportion Modified	Number Target Base in given regions	Expected Number Modified	Observed Number modified	Observed Proportion	Binomial P-value For Enrichment	Binomial P-value for Depletion	Permutation P-value for Enrichment	Permutation P-value for Depletion.

NOTE:
User must provide comma-separated string for what to put in the Feature and Modification columns. E.g. Repeats,5mC or Exons,6mA''')


parser.add_argument('--tableheader', '-H',
                   action='store_true', default=False,
                   help='''Optional.....''')

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
  expNumPeaks = expProportion*numPeaksA
  pValue = sum( scipy.stats.binom.pmf(range(obsNumOverlap, numPeaksA+1), numPeaksA, expProportion) )
      #q=obsNumOverlap, size=numPeaksA, prob=expProportion, lower=FALSE) 
  return {'obs':obsProportion, 'exp':expProportion, 'p':pValue, 'expnum':expNumPeaks}



def readFaCount(name, rm=True):
    with open(name) as fh:
        keys = fh.readline().strip().split()[1:]
        values = fh.readline().strip().split()[1:]
        counts = {keys[i]:float(values[i]) for i in range(len(keys))}
        values = fh.readline().strip().split()[1:]
        prcnts = {keys[i]:float(values[i]) for i in range(len(keys))}

        if rm:
            os.system('rm '+name)
        return {'counts':counts, 'prcnts':prcnts}

def faCount(fasta):
    
    name = 'tmp.'+str(np.random.randint(10000,99999))+'.txt'
    cmd = 'faCount -strands -dinuc -summary ' + fasta + ' > ' + name
    success = os.system(cmd)
    if success == 0:
        return readFaCount(name)
    else:
        assert 'faCount Error'


def faCountTarget(fasta, B):
    fa = 'tmp.'+str(np.random.randint(10000,99999))+'.fasta'
    B.sequence(fasta, fo=fa)
    counts = faCount(fa)
    os.system('rm '+fa)
    return counts
       




#1. Use faCount to get base counts in entire sequence file
sys.stderr.write('1. Bg FaCount \n')
if args.background:
    refCounts = readFaCount(args.background, rm=False)
else:
    refCounts = faCount(args.fasta)

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
    

# 4. Use BEDtools to extract target region sequences and faCount to get base counts in target regions
sys.stderr.write('4. Region FaCount A\n')
regionCounts = faCountTarget(args.fasta, B)

# 5. Find the mods that overlap regions
sys.stderr.write('5. A in B\n')
modsInRegions = A.intersect(B, u=True)



#6. For the given target base, calculate the expected value, observed value, and binomial stats
sys.stderr.write('6. Binomial stats... ')
totalPositions = refCounts['counts'][args.targetbase]
totalObserved = len(A) ## could also be get_sum_lengths(A)
positionsInRegions = regionCounts['counts'][args.targetbase]
observedInRegions = len(modsInRegions)
expectedPropInRegions = float(totalObserved/totalPositions)
expectedNumInRegions = expectedPropInRegions*positionsInRegions
observedPropInRegions = float(observedInRegions)/positionsInRegions
#sys.stderr.write('pValR... ')
#pValueR = sum( scipy.stats.binom.pmf(range(observedInRegions, int(positionsInRegions+1)), positionsInRegions, expectedPropInRegions) )
sys.stderr.write('pValR... ')
pValueR = scipy.stats.binom_test(observedInRegions, positionsInRegions, expectedPropInRegions, 'greater')
#sys.stderr.write('pValL... ')
#pValueL = sum( scipy.stats.binom.pmf(range(0,observedInRegions), positionsInRegions, expectedPropInRegions) )
sys.stderr.write('pValL... ')
pValueL = scipy.stats.binom_test(observedInRegions, positionsInRegions, expectedPropInRegions, 'less' )
sys.stderr.write('pVal2... \n')
pValue2 = scipy.stats.binom_test(observedInRegions, positionsInRegions, expectedPropInRegions )


if not args.table:
    print 'totalPositions', totalPositions
    print 'totalObserved', totalObserved
    print 'positionsInRegions', positionsInRegions
    print 'observedInRegions', observedInRegions
    print 'expectedPropInRegions', expectedPropInRegions
    print 'expectedNumInRegions', expectedNumInRegions
    print 'observedPropInRegions', observedPropInRegions
    print 'pValueR', pValueR
    print 'pValueL', pValueL
    print 'pValue-two-side', pValue2

if args.perm:
    randomObs = []
    sys.stderr.write('Permutations: ')
    for i in range(args.perm):
        sys.stderr.write('.')
        B = B.shuffle(g=args.genome, chrom=args.chrom, noOverlapping=args.noOverlapping).sort().merge()
        regionCountsD = faCountTarget(args.fasta, B)
        regionCounts = regionCountsD['counts'][args.targetbase]
        modsInRegions = A.intersect(B, u=True)
        modCounts = len(modsInRegions)
        randomObs.append( float(modCounts)/regionCounts )

    sys.stderr.write('\n')
    randomObs = np.array(randomObs)
    L = float(len(randomObs))
    NR = sum(randomObs > observedPropInRegions)
    NL = sum(randomObs < observedPropInRegions)
    if not args.table:
        print "NL NR L NL/L NR/L"
        print NL, NR, L, NL/L, NR/L
        print randomObs



if args.table:
    header = ['Feature', 'Target Base', 'Modification', 'Total Target Base', 'Total Observed Target Base Modifications', 'Expected Proportion Modified', 'Number Target Base in given regions', 'Expected Number Modified', 'Observed Number modified', 'Observed Proportion', 'Binomial P-value For Enrichment', 'Binomial P-value for Depletion']
    ftr, mod = args.table.strip().split(',')
    out = [ftr, args.targetbase, mod, int(totalPositions), int(totalObserved), expectedPropInRegions, int(positionsInRegions), expectedNumInRegions, int(observedInRegions), observedPropInRegions, pValueR, pValueL]

    if args.perm:
        header += ['Permutation P-value for Enrichment', 'Permutation P-value for Depletion']
        out += [NR/L, NL/L]
    if args.tableheader:
        print '\t'.join(header)
    #OUTPUT    
    print '\t'.join([str(e) for e in out])
