import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg
from normalize import protocol1, protocol2, protocol3, protocol4, protocol5, protocol6, normalize
from hmm_fxns_for_R import *




def run(parser, args):

    bedgraph = CovBed(args.bedgraph)

    if not args.quiet:
        newmsg("finding state path")


    ## CONSTRUCT EMISSIONS PROBABILITY MATRIX FOR R


    eprobs, nstates = help_get_emission_probs(args.mu, args.sigma, args.mu_scale)


    ## CONSTRUCT TRANSITIONS PROBABILITY MATRIX FOR R


    tprobs = help_get_transition_probs(args.leave_special_state, args.leave_other, args.special_idx, nstates)

    
    ## CONSTRUCT INITIAL PROBABILITY MATRIX FOR R


    iprobs = help_get_initial_probs(nstates, args.special_idx, args.init_special, args.initialprobs)


    ## HIDDEN MARKOV MODEL: Find most probable path through states

    statepath, emitted_data = generate_hmmR(bedgraph, args.emodel, eprobs=eprobs, tprobs=tprobs, iprobs=iprobs)


    ##
    if not args.quiet:
            bdgmsg("state path", False)
            

    sys.stdout.write(bedgraph.expanded_bdg_two_cols(emitted_data, statepath))

        
