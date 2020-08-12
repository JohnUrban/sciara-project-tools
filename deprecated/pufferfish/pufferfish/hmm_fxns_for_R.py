from CovBedClass import *
from collections import defaultdict

## PREVIOUSLY THIS FUNCTION USED VARIABLES DEFINED INSIDE R.
## NEWER VERSION BELOW USES ONLY VARIABLES DEFINED IN PYTHON FIRST.
##def hmm7(late, path, emodel):
##    states = {}
##    if path == 'viterbi':
##        for chrom in late.chromosomes:
####            sys.stderr.write( chrom + "\n" )
##            if len(late.count[chrom]) > 1:
##                v = puffR.viterbi_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
##                states[chrom] = list(v[0])
##            else:
##                states[chrom] = [0] ## if only 1 bin, assign a non-state
##    elif path == 'posterior':
##        for chrom in late.chromosomes:
##            f = puffR.forward_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
##            b = puffR.backward_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
##            p = posterior(f[0], b[0], [1,2,3,4,5,6,7])
##            states[chrom] = list(p[0])            
##    return states


## NEWER VERSION OF HMM7 BELOW USES ONLY VARIABLES DEFINED IN PYTHON -- SEE PROBABILITY MATRICES BELOW.
## IT ALSO ALLOWS SOME TOGGLING OF HMM PARAMETERS
def hmmR(late, path, emodel, eprobs=None, tprobs=None, iprobs=None):
    if eprobs is None:
        eprobs = emissions7()
    if tprobs is None:
        tprobs = transitions7()
    if iprobs is None:
        iprobs = initial7()
    states = intvec(range(1,len(iprobs)+1))
    statepath = {}
    if path == 'viterbi':
        for chrom in late.chromosomes:
##            sys.stderr.write( chrom + "\n" )
            if len(late.count[chrom]) > 1:
                v = puffR.viterbi_puff(emissions = eprobs, transitions = tprobs, initial = iprobs, states = states, emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
                statepath[chrom] = list(v[0])
            else:
                statepath[chrom] = [0] ## if only 1 bin, assign a non-state
    elif path == 'posterior':
        for chrom in late.chromosomes:
            f = puffR.forward_puff(emissions = eprobs, transitions = tprobs, initial = iprobs, states = states, emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
            b = puffR.backward_puff(emissions = eprobs, transitions = tprobs, initial = iprobs, states = states, emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
            p = posterior(f[0], b[0], [1,2,3,4,5,6,7])
            statepath[chrom] = list(p[0])            
    return statepath

def generate_hmmR(late, emodel, eprobs=None, tprobs=None, iprobs=None):
    states = intvec(range(1,len(iprobs)+1))
    statepath = {}
    emitted_data = {}
    for chrom in late.chromosomes:
        statepathlen = len(late.count[chrom])
        if statepathlen > 1:
            ans = puffR.generate(emissions = eprobs, transitions = tprobs, initial = iprobs, states = states, statepathlen = statepathlen, emodel = emodel)
            statepath[chrom] = list(ans[0])
            emitted_data[chrom] = list(ans[1])
        else:
            pass ## do not use that chrom for now
##            statepath[chrom] = [0] ## if only 1 bin, assign a non-state
            
    return statepath, emitted_data


## FXNS TO HELP GENERATE PROB MATRICES FOR R
def get_emission_probs(mu=[1,2,4,8,16,32,64],sig=None):
    mu = fltvec(mu)
    if sig is None:
        sig = fltvec([e**0.5 for e in mu])
    else:
        sig = fltvec(sig)
    eprobmat = matrixr(mu+sig, nrow=2, byrow=True)
    return eprobmat


##def get_transition_probs(nstates=7, changestate=0.001, samestate=0.999):
##    t = np.zeros([nstates,nstates], dtype=float)
##    for i in range(nstates):
##      t[i,] = changestate #0.000001
##      t[i,i] = samestate #0.999999
##      t[i,] = t[i,]/sum(t[i,])
##    rowsvec = []
##    for i in range(nstates):
##        rowsvec += list(t[i,])
##    rowsvecr = fltvec(rowsvec)
##    return matrixr(rowsvecr, nrow=nstates, byrow=True)

def get_transition_probs(nstates=7, special_state_idx=0, leave_special=0.001, leave_non_to_special=0.001, leave_non_to_othernon=0.001):
    stay_special = 1 - (leave_special * (nstates-1))
    stay_non = 1 - leave_non_to_special - (leave_non_to_othernon * (nstates-2))
    nonspecial = range(nstates)
    nonspecial.pop(special_state_idx)
    t = np.zeros([nstates,nstates], dtype=float)
    
    for i in range(nstates):
        ## mark whole row as non to other non
        t[i,] = leave_non_to_othernon
        ## mark cell to special
        t[i, special_state_idx] = leave_non_to_special
        ## mark all diagonal as non to self
        t[i,i] = stay_non

    ## mark entire row of special to other states as leave_special
    t[special_state_idx,] = leave_special
    ## mark special-to-self
    t[special_state_idx, special_state_idx] = stay_special
    ## ensure all rows sum to 1
    for i in range(nstates):
      t[i,] = t[i,]/sum(t[i,])
    ##convert to R
    rowsvec = []
    for i in range(nstates):
        rowsvec += list(t[i,])
    rowsvecr = fltvec(rowsvec)
    return matrixr(rowsvecr, nrow=nstates, byrow=True)


def get_initial_probs(nstates=7, special_state_idx=0, special_state=0.997, other_states=None, inits=None):
    if inits is None:
        if other_states is None:
            ## Each other state is given uniform prob 
            other_states = (1.0-special_state)/float(nstates-1)
        inits = [other_states]*nstates
        inits[special_state_idx] = special_state
    ##    return matrixr( fltvec([0.997] + [0.0005]*6), nrow=1 )
    ## ENSURE it sums to 1
    inits = list( np.array(inits)/np.array(inits).sum() )
    return matrixr( fltvec(inits), nrow=1 )












## THESE ARE THE PROBABILITY MATRICES USED IN THE 7-STATE DNA PUFF FINDING MODEL

def emissions7():
    return get_emission_probs(mu=[1,2,4,8,16,32,64],sig=None)

def transitions7():
    return get_transition_probs(nstates=7, statechange=0.001, samestate=0.999)

def initial7():
    return get_initial_probs(nstates=7, special_state_idx=0, special_state=0.997)




##def ORI_bed(late):
##    ## BED with chr

def get_levels(states):
    levels = {}
    for chrom in states.keys():
        levels[chrom] = 2**(np.array(states[chrom])-1)
    return levels








### These help clean up the sub-program files
def help_get_emission_probs(mu, sigma=None, mu_scale=None):
    '''mu is a string w/ comma-separated state means
        sigma is a string with comma-sep stat stdevs
        RETURNS: e_prob matrix and nstates'''
    ## emissions: determine state means
    e_mu = [float(e) for e in mu.strip().split(',')]
    ## emissions: determine number of states from state means
    nstates = len(e_mu)
    ## emissions: determine state sigmas
    if sigma is not None:
        e_sig = [float(e) for e in sigma.strip().split(',')]
    elif mu_scale is not None: 
        e_sig = [e*mu_scale for e in e_mu]
    else:
        e_sig = [e**0.5 for e in e_mu]
    assert len(e_sig) == nstates
    ## emissions: get R object
    eprobs = get_emission_probs(mu=e_mu, sig=e_sig)
    return eprobs, nstates


def help_get_transition_probs(leave_special_state, leave_other, special_state_idx, nstates):
    '''
        leave_special_state is probability of leaving special state (e.g. CN=1) - can make it same as others.
        leave_other is probability of leaving a non-special state
        special_state_idx is the 0-based idx of where to find special state params
        nstates is number of states in model
    '''
    if leave_special_state > 1:
        leave_special_state = 1.0/leave_special_state
    if leave_other is None:
        leave_non_to_special = leave_special_state
        leave_non_to_other = leave_special_state
    else:
        leave_other = [float(e) for e in leave_other.split(',')]
        if len(leave_other) == 1:
            if leave_other[0] > 1:
                lp = 1.0/leave_other[0]
            else:
                lp = leave_other[0]
            leave_non_to_special = lp
            leave_non_to_other = lp
        elif len(leave_other) > 1:
            if leave_other[0] > 1:
                leave_non_to_special = 1.0/leave_other[0]
            else:
                leave_non_to_special = leave_other[0]
            if leave_other[1] > 1:
                leave_non_to_other = 1.0/leave_other[1]
            else:
                leave_non_to_other = leave_other[1]

    tprobs = get_transition_probs(nstates=nstates, special_state_idx=special_state_idx, leave_special=leave_special_state, leave_non_to_special= leave_non_to_special, leave_non_to_othernon=leave_non_to_other)
    return tprobs


def help_get_initial_probs(nstates, special_state_idx, init_special, initialprobs=None):
    if initialprobs is None:
        iprobs = get_initial_probs(nstates=nstates, special_state_idx=special_state_idx, special_state=init_special)
    else:
        inits = [float(e) for e in initialprobs.strip().split(',')]
        iprobs = get_initial_probs(inits=inits)
    assert len(iprobs) == nstates
    return iprobs


def help_get_prob_matrices_from_params(mu, sigma, mu_scale, leave_special_state, leave_other, special_idx, init_special, initialprobs):
    ## CONSTRUCT EMISSIONS PROBABILITY MATRIX FOR R
    eprobs, nstates = help_get_emission_probs(mu, sigma, mu_scale)

    ## CONSTRUCT TRANSITIONS PROBABILITY MATRIX FOR R
    tprobs = help_get_transition_probs(leave_special_state, leave_other, special_idx, nstates)
    
    ## CONSTRUCT INITIAL PROBABILITY MATRIX FOR R
    iprobs = help_get_initial_probs(nstates, special_idx, init_special, initialprobs)

    return eprobs, tprobs, iprobs



def k_state_means(data, k):
    return kmeans(x=data, centers=k)

def get_k_mean_centers(km):
    ## km is kmeans object
    return list(km[1])

def get_k_mean_sigmas(km):
    ## km is kmeans object
    d = defaultdict(int)
    for k in km[0]:
        d[k-1] += 1 ## k-1 to put it in python indexing
    sigs = []
    for k in sorted(d.keys()):
        sig = (km[3][k]/d[k])**0.5
        sigs.append( sig )
    return sigs

def get_k_mean_cluster_lengths(km):
    d = defaultdict(list)
    last_c = km[0][0]
    L = 1
    for i in range(1,len(km[0])):
        this_c = km[0][i]
        if this_c == last_c:
            L+=1
        else:
            d[last_c].append(L)
            last_c = this_c
            L = 1
    d[this_c] = L
    return d ## all lengths found

def get_k_mean_cluster_mean_lengths(km_lengths):
    means = []
    for c in sorted(km_lengths.keys()):
        means.append( np.array(km_lengths[c]).mean() )
    return means

def get_k_mean_cluster_trans_probs(kmmeans):
    leave = 1/np.array(kmmeans, dtype=float)
    return list(1-(leave)), list(leave/(len(leave)-1))

def get_k_mean_cluster_initial_probs(km_lengths):
    sums = []
    for c in sorted(km_lengths.keys()):
        sums.append( np.array(km_lengths[c]).sum() )
    return list( np.array(sums, dtype=float)/np.array(sums,dtype=float).sum() )

## this method is not yet implemented in full work flow
def get_k_mean_cluster_trans_probs_matrix(km):
    t = np.ones([len(km[1]), len(km[1])]) ## all start w/ pseudo count of 1
    for i in range(1,len(km[0])):
        c_from = km[0][i-1] - 1 #for python array indexing
        c_to = km[0][i] - 1 #for python array indexing
        t[c_from, c_to] += 1
    for i in range(len(km[1])):
        t[i,] = t[i,]/t[i,].sum()
    return t



## TODO: update this to allow better transition estimates. e.g. can just make i,j matrix of all i-->j
def get_state_emissions_from_kmeans(data, k):
    km = k_state_means(data, k)
    mu = get_k_mean_centers(km)
    sig = get_k_mean_sigmas(km)
    km_lengths = get_k_mean_cluster_lengths(km)
    init_probs = get_k_mean_cluster_initial_probs(km_lengths)
    mean_lengths = get_k_mean_cluster_mean_lengths(km_lengths) #in order of clustnums like other2
    stay_probs, leave_probs = get_k_mean_cluster_trans_probs(mean_lengths)
    sorted_params = sorted( zip( mu, sig , init_probs, stay_probs, leave_probs) )
    return [mu for mu,sig,I,S,L in sorted_params], [sig for mu,sig,I,S,L in sorted_params], [I for mu,sig,I,S,L in sorted_params], [S for mu,sig,I,S,L in sorted_params], [L for mu,sig,I,S,L in sorted_params]

def help_get_state_emissions_from_kmeans(data, k):
    e_mu, e_sig, inits, stay, leave = get_state_emissions_from_kmeans(data, k)
    eprobs = get_emission_probs(mu=e_mu, sig=e_sig)
    rowsvec = []
    for i in range(k):
        row = [leave[i]]*k
        row[i] = stay[i]
        rowsvec += row
    rowsvecr = fltvec(rowsvec)
    tprobs = matrixr(rowsvecr, nrow=k, byrow=True)
    iprobs = matrixr( fltvec(inits), nrow=1 )
    return eprobs, tprobs, iprobs
                                 
                            



## FUNCTIONS NOT BEING USED
def window_medians(x, flanksize=2, includeflanks=True):
    ''' Given vector of FEs and flank size, find median FE in window centered over each bin by entending flanksize in each direction.
    TODO: For edge cases, take median in windows of flanksize*2 by adding/subtracting on each side as nec...'''
    meds = []
    lenx = len(x)
    start = flanksize
    end = lenx-flanksize
    Pos = start
    ans = np.median( x[Pos-flanksize:Pos+1+flanksize] )
    if includeflanks:
        # Append ans for all flank positions before Start
        meds += [ans]*flanksize
    meds.append(ans)
    #Iter from start+1
    for Pos in range(start+1, end):
        meds.append( np.median( x[Pos-flanksize:Pos+1+flanksize] ) )
    if includeflanks:
        # Append ans for all flank positions before Start
        meds += [meds[-1]]*flanksize
    return np.array(meds)
    

def window_sums(x, flanksize=2, returnmeans=False, includeflanks=True):
    ''' Given vector and flank size, get sum or mean of bin and flanks to each side.
        Note that this is sliding window mean-smoothing with uniform weighting.
        Can use ksmooth or loess in R to weight closer bins higher.'''
    sums = []
    lenx = len(x)
    start = flanksize
    end = lenx-flanksize 
    windowsize = float(flanksize*2 + 1) ## 2 flanks + Pos it is centered on
    Pos = start
    #First window
    ans = np.sum( x[Pos-flanksize:Pos+flanksize+1] )
    if includeflanks:
        # Append ans for all flank positions before Start
        sums += [ans]*flanksize
    # Append ans for Start
    sums.append( ans ) 

    #Iterate
    for Pos in range(start+1, end):
        ## Subtract first element of last window
        ans -= x[Pos-1-flanksize]
        ## Add last element of current window
        ans += x[Pos+flanksize]
        ## Append ans to sums
        sums.append( ans )

    if includeflanks:
        # Append last answer for all flank positions after End
        sums += [sums[-1]] * flanksize
        
    ## convert sums to np
    sums = np.array(sums)
    ## If means desired
    if returnmeans:
        return sums/windowsize
    return sums


def ksmooth_counts(self, bw=10000):
    for chrom in self.chromosomes:
        x = self.start[chrom]
        y = self.count[chrom]
        k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
        self.count[chrom] = np.array(k[1])

## HMM followed by window_modes of states could help...
