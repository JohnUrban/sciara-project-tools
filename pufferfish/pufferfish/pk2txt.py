import sys, os, gzip
import datetime
from CovBedClass import *
import cPickle as pickle
from loaddata import load_pickle

def newmsg(x):
    sys.stderr.write(str(datetime.datetime.now()) +": .. " + x + "...\n")
    
def bdgmsg(x, collapsed=False):
    if not collapsed:
        return sys.stderr.write(str(datetime.datetime.now()) +": ..printing " + x + " as expanded single-step bedGraph...\n")
    else:
        return sys.stderr.write(str(datetime.datetime.now()) +": ..printing " + x + " as collapsed var-step bedGraph...\n")
    
def run(parser, args):
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..Un-pickling...\n")
    f = load_pickle(args.inpickle)

    if args.counts:
        if not args.quiet:
            bdgmsg(x="(normalized) counts from all files")
        o = open(args.prefix+".counts.bedGraph",'w')
        o.write(f.get_counts_bdg())
        o.close()

    if args.filtered_counts:
        if not args.quiet:
            bdgmsg(x="counts from filtered contigs from all files")
        o = open(args.prefix+".filtered_counts.bedGraph",'w')
        o.write(f.get_filtered_contigs_bdg())
        o.close()
        
    if args.viterbi_expanded:
        if not args.quiet:
            bdgmsg(x="viterbi state path")
        o = open(args.prefix+".viterbi_expanded.bedGraph",'w')
        o.write(f.get_cor_state_bdg())
        o.close()
        
    if args.viterbi_collapsed:
        if not args.quiet:
            bdgmsg(x="viterbi state path", collapsed=True)
        o = open(args.prefix+".viterbi_collapsed.bedGraph",'w')
        o.write(f.get_cor_state_bdg(collapsed=True))
        o.close()
        
    if args.correlations_expanded:
        if not args.quiet:
            bdgmsg(x="bin/stage correlation scores")
        o = open(args.prefix+".corscore_expanded.bedGraph",'w')
        o.write(f.get_corscore_bdg())
        o.close()
        
    if args.correlations_collapsed:
        if not args.quiet:
            bdgmsg(x="bin/stage correlation scores", collapsed=True)
        o = open(args.prefix+".corscore_collapsed.bedGraph",'w')
        o.write(f.get_corscore_bdg(collapsed=True))
        o.close()
        
    if args.smoothed_correlations_expanded:
        if not args.quiet:
            bdgmsg(x="smoothed bin/stage correlation scores")
        o = open(args.prefix+".smoothed_corscore_expanded.bedGraph",'w')
        o.write(f.get_smooth_corscore_bdg())
        o.close()
        
    if args.smoothed_correlations_collapsed:
        if not args.quiet:
            bdgmsg(x="smoothed bin/stage correlation scores", collapsed=True)
        o = open(args.prefix+".smoothed_corscore_collapsed.bedGraph",'w')
        o.write(f.get_smooth_corscore_bdg(collapsed=True))
        o.close()

    if args.slopes_expanded:
        if not args.quiet:
            bdgmsg(x="bin/stage slopes")
        o = open(args.prefix+".slope_expanded.bedGraph",'w')
        o.write(f.get_slope_bdg())
        o.close()
        
    if args.slopes_collapsed:
        if not args.quiet:
            bdgmsg(x="bin/stage slopes", collapsed=True)
        o = open(args.prefix+".slope_collapsed.bedGraph",'w')
        o.write(f.get_slope_bdg(collapsed=True))
        o.close()
    


