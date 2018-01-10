#!/usr/bin/env python2.7
import sys, argparse, pybedtools, scipy.stats
from collections import defaultdict
import numpy as np
from CovBedClass import *
from pk2txt import bdgmsg, newmsg



parser = argparse.ArgumentParser(description="""

 Given:
     (i) HMM STATE bedGraph from pufferfish
     (ii) Target signal bedGraph to normalize
 Find mean of each state.
 Return (and/or):
     (i) Normalized target bedGraph -- where each target value is divided by the mean of its HMM state
         The desired outcome is that the target signal becomes centered on CN=1.
         Downstream steps can then make the assumption of CN=1 for all bins.
     (ii) Mean level bedGraph where each bin is assigned the mean of its state.
          Can be used to:
              - visually compare to original signal
              - use awk on target bdg vs mean level bdg downstream to normalize, signif test, or other


    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--signal', '-i', '-f',
                   type= str,
                   help='''Path to signal bedGraph that usually has cols: chr, start, end, signal-to-correct.
Can tell algo what col to look at.''')

parser.add_argument('--states', '-i2', '-f2', 
                   type=str, 
                   help='''Path to bedGraph that usually has cols: chr, start, end, HMM state.
Can tell algo what col to look at.
chr/start/end should be identical in values and sort order as --signal.''')

parser.add_argument('--signalcol', '-s', 
                   type=int, default=4,
                   help='''1-based column that signal found in. Default = 4''')
parser.add_argument('--statecol', '-S',
                   type=int, default=4,
                   help='''1-based column that signal found in. Default = 4''')

parser.add_argument('--chrcol', 
                   type=int, default=1,
                   help='''1-based column that chr/seq name found in. Default = 1''')
parser.add_argument('--startcol',
                   type=int, default=2,
                   help='''1-based column that start coordinate found in. Default = 2''')
parser.add_argument('--endcol', 
                   type=int, default=3,
                   help='''1-based column that end coordinate found in. Default = 3''')

parser.add_argument('--levels',
                   type= str,
                   help='''By default, the CN normalized bedGraph is written to stdout.
Using this flag and providing a file name tells the program to also write
a bedGraph to that file name for the means of the states over each bin. 
''')

parser.add_argument('--levels_only',
                   action='store_true', default=False,
                   help='''By default, the CN normalized bedGraph is written to stdout.
Using this flag tells the program to only return the levels bedGraph (to stdout by default if --levels not used). 
''')

parser.add_argument('--normbdg',
                   type= str, default=False,
                   help='''By default, the CN normalized bedGraph is written to stdout.
This redirects it into a filename provided. 
''')

parser.add_argument('-c', '--collapsed', action='store_true', default=False,
                           help='''Return collapsed variable-step bedGraph instead of expanded single-step bedGraph.
This is often a much smaller file.''')

parser.add_argument('-q', '--quiet', action='store_true', default=False,
                           help='''QUIET.''')


args = parser.parse_args()



sigcol = args.signalcol-1
statecol = args.statecol-1
chrcol = args.chrcol-1
startcol = args.startcol-1
endcol = args.endcol-1









##def run(parser, args):

if not args.quiet:
    sys.stderr.write(str(datetime.datetime.now()) +": ..Loading files...\n")

signal = CovBed(args.signal)
states = CovBed(args.states, count_only=True)

## FIRST GET STATE MEANS
if not args.quiet:
    sys.stderr.write(str(datetime.datetime.now()) +": ..Learning state means...\n")

## MEANS: Sum up data over each state and tally the number of times the state was observed.
statesum = defaultdict(float) ## Sum of signal over each state
stateobs = defaultdict(float) ## Number times state is observed.
for chrom in states.chromosomes:
    numbins = len(states.count[chrom])
    for i in range(numbins):
        state = states.count[chrom][i]
        emission = signal.count[chrom][i]
        statesum[state] += emission
        stateobs[state] += 1.0
## MEANS: Divide Sum of data over each state by number of times the state was observed.
statemeans = defaultdict(float)
for state, n in stateobs.iteritems():
    statemeans[state] = statesum[state] / float(n)

## NORMALIZE BY STATE MEANS
if not args.quiet:
    sys.stderr.write(str(datetime.datetime.now()) +": ..Normalizing emissions to state means...\n")

## NORM: Create dictionary that contains chroms w/ lists of signal/statemean
normsig = defaultdict(list)
levels = defaultdict(list)
for chrom in signal.chromosomes:
    numbins = len(signal.count[chrom])
    for i in range(numbins):
        state = states.count[chrom][i]
        emission = signal.count[chrom][i]
        statemean = statemeans[state]
        levels[chrom].append( statemean )
        normsig[chrom].append( emission / statemean )

## OUTPUT:
if not args.levels_only:
    if args.normbdg:
        normsigout = open(args.normbdg, 'w')
    else:
        normsigout = sys.stdout
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..Writing normalized signal bedGraph...\n")
    normsigout.write(signal.get_bdg(normsig, args.collapsed))
    if args.normbdg:
        normsigout.close()


if args.levels or args.levels_only:
    if args.levels:
        levelsout = open(args.levels, 'w')
    elif args.levels_only: ## "levels_only and levels" already taken care of by first condition
        levelsout = sys.stdout
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..Writing state level means bedGraph...\n")
    levelsout.write(signal.get_bdg(levels, args.collapsed))
    if args.levels:
        levelsout.close()










