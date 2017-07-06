import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg
from normalize import protocol1, protocol2, protocol3, protocol4, protocol5, protocol6, normalize

def hmm7(late, path, emodel):
    states = {}
    if path == 'viterbi':
        for chrom in late.chromosomes:
##            sys.stderr.write( chrom + "\n" )
            if len(late.count[chrom]) > 1:
                v = puffR.viterbi_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
                states[chrom] = list(v[0])
            else:
                states[chrom] = [0] ## if only 1 bin, assign a non-state
    elif path == 'posterior':
        for chrom in late.chromosomes:
            f = puffR.forward_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
            b = puffR.backward_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
            p = posterior(f[0], b[0], [1,2,3,4,5,6,7])
            states[chrom] = list(p[0])            
    return states

##def ORI_bed(late):
##    ## BED with chr

def get_levels(states):
    levels = {}
    for chrom in states.keys():
        levels[chrom] = 2**(np.array(states[chrom])-1)
    return levels


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
    
    late = normalize(latestage=args.latestage, protocol=protocol, earlystage=args.earlystage, pseudo=args.pseudo, bandwidth=args.bandwidth, quiet=args.quiet)

    if args.counts:
        if not args.quiet:
            bdgmsg("final normalized late stage counts", args.collapsed)
        o = open(args.counts + ".bedGraph", 'w')
        o.write(late.get_bdg(late.count, args.collapsed))
        o.close()
        
    if not args.quiet:
        newmsg("finding state path")
    states = hmm7(late, args.path, args.emodel)
    if args.levels:
        newmsg("getting levels")
        levels = get_levels(states)

    if not args.quiet:
            bdgmsg("state path", args.collapsed)
    if args.levels:
        sys.stdout.write(late.get_bdg(levels, args.collapsed))
    else:
        sys.stdout.write(late.get_bdg(states, args.collapsed))
        
