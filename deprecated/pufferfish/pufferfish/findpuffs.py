import sys, os, gzip
import datetime
from CovBedClass import *
import cPickle as pickle
from loaddata import load_pickle



def run(parser, args):
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..initializing...\n")
        
    if args.input:
        if not args.quiet:
            sys.stderr.write(str(datetime.datetime.now()) +": ..Loading files...\n")
        f = MultiCovBed(args.input)
    elif args.inpickle:
        if not args.quiet:
            sys.stderr.write(str(datetime.datetime.now()) +": ..Un-pickling...\n")
        f = load_pickle(args.inpickle)

    f.filter_null_contigs(args.max_empty_bin_pct, args.max_offending_samples_pct) ## filter bad contigs

    if args.smoothbeforemedian:
        f.ksmooth_counts(bw=args.smoothbeforemedian) ## smooth before finding median

##    if args.mediannorm:
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..median normalizing...\n")
    f.median_normalize_data()

    if args.smoothbeforecor:
        f.ksmooth_counts(bw=args.smoothbeforecor) ## smooth counts before getting cor in bins - prob should not do both b4median and b4cor

##    if args.calcbincors:
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..calc corscores...\n")
    f.cor_score()

##    if args.smoothcors:
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..ksmoothing corscores...\n")
    f.ksmooth_corscores(bw=args.corsmoothbandwidth)


    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..getting viterbi state paths...\n")
    f.get_cor_states()


    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..getting slopes in bins...\n")
    f.find_slopes()
    
    if not args.quiet:
        sys.stderr.write(str(datetime.datetime.now()) +": ..Pickling...\n")
    if not args.outpickle.endswith('.gz'):
        args.outpickle += '.gz'
    if os.path.exists(args.outpickle):
        os.remove(args.outpickle)
    with gzip.open(args.outpickle,'wb') as pkfile:
        pickle.dump(f, pkfile)


##sys.stderr.write(str(datetime.datetime.now()) +": ..printing meta-bdg with median norm scores from all files...\n")
##print f
##print f.get_corscor_bdg()
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing ksmoothed corscores...\n")
##print f.get_smooth_corscor_bdg()
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing viterbi state paths...\n")
##print f.get_cor_state_bdg()
##sys.stderr.write(str(datetime.datetime.now()) +": ..analyzing 0-state bins...\n")
##f.analyze_state0_bins()
##print f.state0_medians



