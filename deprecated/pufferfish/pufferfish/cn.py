import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg
from normalize import protocol1, protocol2, protocol3, protocol4, protocol5, protocol6, normalize
from hmm_fxns_for_R import *



## TODO:
## ALLOW OPTION TO USE KMEANS CLUSTERING TO FIND MEAN OF K STATES (w/ or w/o normalization)
## --kmeans 5 -- this will over-ride the normal emissions params given or used and will learn them...
## might also want to give option to randomly sample from data if kmeans takes to long on all data..

def run(parser, args):
##    if not args.quiet:
##        newmsg("loading late stage file")
##    late = CovBed(args.latestage)
##    if args.earlystage:
##        if not args.quiet:
##            newmsg("loading early stage file")
##        early = CovBed(args.earlystage)
##    else:
##        early = False
##    ## todo -- add filtering out 0 contigs option...
##    
##    ##only protocol1 for now
##    if not args.quiet:
##        if args.earlystage:
##            emsg = ' with additional early stage normalization'
##        else:
##            emsg = ' without additional early stage normalization'
##        if args.protocol1:
##            newmsg("following normalization protocol 1"+emsg)
##        elif args.protocol2:
##            newmsg("following normalization protocol 2"+emsg)
##        elif args.protocol3:
##            newmsg("following normalization protocol 3"+emsg)
##        elif args.protocol4:
##            newmsg("following normalization protocol 4"+emsg)
##    if args.protocol1:
##        late = protocol1(late, early, args.pseudo)
##    elif args.protocol2:
##        late = protocol2(late, early, args.bandwidth, args.pseudo)
##    elif args.protocol3:
##        late = protocol3(late, early, args.bandwidth, args.pseudo)
##    elif args.protocol4:
##        late = protocol4(late, early, args.bandwidth, args.pseudo)
##
    ## TODO - just change to 1 argument: --protocol -- with options [1,2,3,4]

    
    if args.protocol1:
        protocol=1
    elif args.protocol2:
        protocol=2
    elif args.protocol3:
        protocol=3
    elif args.protocol4:
        protocol=4
    elif args.protocol5:
        protocol=5
    elif args.protocol6:
        protocol=6
    elif args.protocol7:
        protocol=7
    
    late = normalize(latestage=args.latestage, protocol=protocol, earlystage=args.earlystage, pseudo=args.pseudo, bandwidth=args.bandwidth, quiet=args.quiet, impute=args.impute, replace=args.replace, replace_with=args.replace_with, replace_this=args.replace_this)

    if args.counts:
        if not args.quiet:
            bdgmsg("final normalized late stage counts", args.collapsed)
        o = open(args.counts + ".bedGraph", 'w')
        o.write(late.get_bdg(late.count, args.collapsed))
        o.close()
        
    if not args.quiet:
        newmsg("Constructing probability matrices...")


    ## CONSTRUCT EMISSIONS/TRANSITION/INITIAL PROBABILITY MATRIXES FOR R
    if args.kmeans is not None:
        data = []
        for chrom in late.count:
            data += list(late.count[chrom])
            
        eprobs, tprobs, iprobs = help_get_state_emissions_from_kmeans(data, args.kmeans)
    else:
        eprobs, tprobs, iprobs = help_get_prob_matrices_from_params(args.mu, args.sigma, args.mu_scale, args.leave_special_state, args.leave_other, args.special_idx, args.init_special, args.initialprobs)

    if not args.quiet:
        newmsg("\nTransition Probs:\n"+str(tprobs))
        newmsg("\nEmission Probs:\n"+str(eprobs))
        newmsg("\nInitial Probs:\n"+str(iprobs))
        
    if not args.quiet:
        newmsg("finding state path")

    ## HIDDEN MARKOV MODEL: Find most probable path through states

    statepath = hmmR(late, args.path, args.emodel, eprobs=eprobs, tprobs=tprobs, iprobs=iprobs)

    ## If opted for, get the state means/levels
    if args.levels:
        newmsg("getting levels")
        levels = get_levels(statepath)

    ##
    if not args.quiet:
            bdgmsg("state path", args.collapsed)
            
    if args.levels:
        sys.stdout.write(late.get_bdg(levels, args.collapsed))
    else:
        sys.stdout.write(late.get_bdg(statepath, args.collapsed))
        
