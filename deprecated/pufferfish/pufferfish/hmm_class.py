import itertools
import numpy as np
import time
from collections import defaultdict
from scipy.stats import norm, expon, poisson, geom

def get_normal_fxn(self):
    return norm(self.emission_probs[0,:], self.emission_probs[1,:])

def get_expon_fxn(self):
    return

def get_pois_fxn(self):
    return

def get_geom_fxn(self):
    return


puff.exponential <- function(x, emissions){
  dexp(x = x, rate = emissions[1, ], log = FALSE)
}

puff.poisson <- function(x, emissions){
  dpois(x = round(x), lambda = emissions[1, ], log = FALSE)
}

puff.geometric <- function(x, emissions){
  dgeom(x = round(x), prob = emissions[1, ], log = FALSE)
}



class HMM(object):

    def __init__(self, states=None, emissions=None):
        # data
        self.states = states
        self.statepath = None
        self.initial_probs = None
        self.transition_probs = None
        self.emission_probs = None
        self.emissions = None
        self.num_states = None
        self.num_emits = None

        ## FOR NOW - DO NOT AUTOMATICALLY GENERATE AT RANDOM...
##        if emissions is None: # initialize with random values
##            if self.states is None:
##                self.states = STATES_FIVEMERS
##            self.randomize_emission_probs()
##            self.randomize_transition_probs()
##            self.randomize_initial_probs()
##            self.randomize_statepath()
##            self.randomize_emissions()
##        else:
##            self.emissions = emissions
            

    def add_initial_probs(self, initial_probs):
        ## should be 1 x nstates np.array
        assert len(initial_probs) == len(self.states)
        self.initial_probs = initial_probs
        

    def add_transition_probs(self, transition_probs):
        ## must be nstate x nstate np.array
        self.transition_probs = transition_probs

    def add_emission_probs(self, emission_probs):
        ## should be be 2 x nstate np.array (if normal -- need means and stdevs)
        self.emission_probs = emission_probs

    def add_states(self, states):
        self.states = states

    def get_initial_probs(self):
        return self.initial_probs

    def get_transition_probs(self):
        return self.transition_probs

    def get_emission_probs(self):
        return self.emission_probs

    def randomize_initial_probs(self, uniform=True):
        nstates = len(self.states)
        if uniform:
            self.initial_probs = [1.0 / nstates] * nstates
        else:
            self.initial_probs = np.random.poisson(lam=10.0, size=nstates)
            self.initial_probs /= sum(self.initial_probs)

    def randomize_transition_probs(self, allow_gaps=True):
        """
        If allow_gaps = False, assumes each k-mer has another kmer overlapped by k-1
        """
        k = len(self.states[0])
        if k > 2:
            nstates = len(self.states)
            self.transition_probs = np.zeros([nstates, nstates])
            # make prefix-suffix dict -- overlaps of k-1 and k-2
            prefix = defaultdict(list)
            for i in xrange(nstates):
                pref = self.states[i][:k-1]
                prefix[pref].append(i)
                pref = self.states[i][:k-2]
                prefix[pref].append(i)
            # create transition probs -- can soft code the sampling parameters later if want
            for i in xrange(nstates):
                ## overlap by k-1 (move = 1)
                current_suffix = self.states[i][1:]
                poisson = np.random.poisson(lam=365.0, size=k-1)
                for t, j in enumerate(prefix[current_suffix]):
                    self.transition_probs[i,j] = poisson[t]
                if allow_gaps:
                    ## overlap by k-2 (move = 2) -- add additional counts
                    current_suffix = self.states[i][2:]
                    poisson = np.random.poisson(lam=4.0, size=(k-1)**2)
                    for t, j in enumerate(prefix[current_suffix]):
                        self.transition_probs[i,j] += poisson[t]
                    ## stay in place: add additional probability to staying in place (move = 0)
                    current_suffix = self.states[i]
                    poisson = np.random.poisson(lam=20.0, size=1)
                    self.transition_probs[i,i] += poisson
                ## normalize all counts by sum to create probs that sum to 1
                self.transition_probs[i,:] /= sum(self.transition_probs[i,:])

    def randomize_emission_probs(self, mu_mean=65.5, mu_sd=6.5, sigma_mean=1.1, sigma_sd=0.4):
        ## generates n state mus from normal dist using mu_mean and mu_sd
        ## generates n state sigmas from normal dist using sigma_mean and sigma_dist
        ## result is n pairs of state mus and sigmas.
        nstates = len(self.states)
        self.emission_probs = np.zeros([2, nstates])
        for i in xrange(nstates):
            self.emission_probs[0,i] = np.random.normal(mu_mean, mu_sd)
            self.emission_probs[1,i] = abs(np.random.normal(sigma_mean, sigma_sd))

    def randomize_statepath(self, length=10):
        nstates = len(self.states)
        self.statepath = [np.random.choice(nstates, p=self.initial_probs)]
        for _ in xrange(length-1):
            self.statepath.append(
                    np.random.choice(
                            nstates,
                            p=self.transition_probs[self.statepath[-1]]))


    def randomize_emissions(self):
        means = self.emission_probs[0, self.statepath]
        stdevs = self.emission_probs[1, self.statepath]
        self.emissions = np.random.normal(means, stdevs)


    def get_num_states(self):
        if self.num_states == None:
            self.num_states = len(self.states)
        return self.num_states

    def get_num_emits(self):
        if self.num_emits == None:
            self.num_emits = len(self.emissions)
        return self.num_emits
    
    def forward(self):
        nstates = self.get_num_states()
        nemits = self.get_num_emits()
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        Forward = np.zeros([nstates, nemits])
        scalefactors = np.zeros([2, nemits])
        # initial
        Forward[:, 0] = np.multiply(self.initial_probs, ep.pdf(self.emissions[0]))
        # scale to prevent underflow -- keep track of scaling
        scalefactors[0,0] = sum(Forward[:,0])
        scalefactors[1,0] = np.log(scalefactors[0,0])
        Forward[:,0] /= scalefactors[0,0]
        # iterate
        for k in xrange(1, nemits):
            emit = ep.pdf(self.emissions[k])
            Forward[:,k] = np.multiply(emit, np.dot(Forward[:,k-1], self.transition_probs))
            scalefactors[0,k] = sum(Forward[:,k])
            scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k-1]
            Forward[:,k] /= scalefactors[0,k]
        self.Forward = Forward
        self.Forward_scalefactors = scalefactors

    def get_forward_matrix(self):
        try:
            return self.Forward
        except:
            self.forward()
            return self.Forward

    def backward(self):
        nstates = self.get_num_states()
        nemits = self.get_num_emits()
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        Backward = np.zeros([nstates, nemits])
        scalefactors = np.zeros([2, nemits])
        end = nemits - 1
        # initial
        Backward[:, end] = 1
        # scale to prevent underflow -- keep track of scaling
        scalefactors[0,end] = sum(Backward[:,end])
        scalefactors[1,end] = np.log(scalefactors[0,end])
        Backward[:,end] /= scalefactors[0,end]
        # iterate
        for k in xrange(end-1, -1, -1):
            emit = ep.pdf(self.emissions[k+1])
            a = np.multiply(Backward[:,k+1], emit).transpose()
            Backward[:,k] = np.dot(self.transition_probs, a).transpose()
            scalefactors[0,k] = sum(Backward[:,k])
            scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k+1]
            Backward[:,k] /= scalefactors[0,k]
        self.Backward = Backward
        self.Backwars_scalefactors = scalefactors

    def get_backward_matrix(self):
        try:
            return self.Backward
        except:
            self.backward()
            return self.Backward

    def posterior_decoding(self):
        ##F and B are scaled long seq matrices -- the scales are scalefactors that come with them out of long fxns
        ## ensure F and B available
        self.get_forward_matrix()
        self.get_backward_matrix()
        nstates = self.get_num_states()
        nemits = np.shape(self.Forward)[1]
        self.posterior_path = np.zeros(nemits, dtype=int)
        for i in xrange(nemits):
            fb = self.Forward[:,i] * self.Backward[:,i]
            self.posterior_path[i] = int(fb.argmax())
        return self.posterior_path

    def prob_data(self, nemits=None):
        #ensure forward available
        self.get_forward_matrix()
        if nemits == None:
            end = np.shape(self.Forward)[1]-1
        else:
            end = nemits-1
        return sum(self.Forward[:,end])*np.exp(self.Forward_scalefactors[1,end])

##    @profile

    def viterbi(self):
        np.seterr(divide='ignore')
        nstates = self.get_num_states()
        nemits = self.get_num_emits()
        initial_probs = np.log(self.initial_probs)
        tran_probs = np.log(self.transition_probs)
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        pointer = np.zeros([nemits, nstates])
        Viterbi = np.zeros([nstates, nemits])  
        ## need to add log_probs instead of multiply probs to prevent underflow
        Viterbi[:,0] = self.initial_probs + ep.logpdf(self.emissions[0])
        pointer[0,:] = 1
        for j in range(1,nemits):
            selection = Viterbi[:,j-1] + tran_probs.transpose() 
            maxstates = np.apply_along_axis(max_and_index, 1, selection)
            Viterbi[:,j] = ep.logpdf(self.emissions[j]) + maxstates[:,1]
            pointer[j,:] = maxstates[:,0]
        end = nemits - 1
        #path init
        viterbi_path = np.zeros(nemits).astype(int)
        viterbi_path[end] = Viterbi[:,end].argmax()
        #prob
        viterbi_prob = Viterbi[viterbi_path[end], end]
        #path iter
        for j in range(end,0,-1):
            viterbi_path[j-1] = pointer[j,viterbi_path[j]]
        return viterbi_path, viterbi_prob


    def compare_statepath(self, dst, src=None):
        if src is None:
            src = self.statepath
        ident = sum(a==b for a,b in itertools.izip(dst, src))
        edit_dist = len(src) - ident
        return edit_dist, ident, 100.0*ident/len(src)


    def baumwelch(self):
        pass


##    def randomize_emissions_twoemits(self):
##        pass

##    def viterbi2(self):
##        np.seterr(divide='ignore')
##        num_states = self.get_num_states()
##        num_emits = self.get_num_emits()
##        initial_probs = np.log(self.initial_probs)
##        tran_probs = np.log(self.transition_probs)
##        ep1 = norm(emission_probs[0,:], emission_probs[1,:])
##        ep2 = norm(emission_probs2[0,:], emission_probs2[1,:])
##        pointer = np.zeros([num_emits, num_states])
##        Viterbi = np.zeros([num_states, num_emits])  
##        ## need to add log_probs instead of multiply probs to prevent underflow
##        Viterbi[:,0] = initial_probs + ep1.logpdf(emitted_data[0]) + ep2.logpdf(emitted_data2[0])
##        pointer[0,:] = 1
##        for j in range(1,num_emits):
##            selection = Viterbi[:,j-1] + tran_probs.transpose() 
##            maxstates = np.apply_along_axis(max_and_index, 1, selection)
##            Viterbi[:,j] = ep1.logpdf(emitted_data[j]) + ep2.logpdf(emitted_data2[j]) + maxstates[:,1]
##            pointer[j,:] = maxstates[:,0]
##        end = num_emits - 1
##        #path init
##        viterbi_path = np.zeros(num_emits).astype(int)
##        viterbi_path[end] = Viterbi[:,end].argmax()
##        #prob
##        viterbi_prob = Viterbi[viterbi_path[end], end]
##        #path iter
##        for j in range(end,0,-1):
##            viterbi_path[j-1] = pointer[j,viterbi_path[j]]
##        return viterbi_path, viterbi_prob
##


#### FUNCTION
def get_emiss_probs_from_model(model, twoemits=False):
    ''' model is object returned from get_stored_model() in model_tools '''
    states = sorted(model[1].keys())
    num_states = len(states)
    t_emissions = np.zeros([2,num_states])
    c_emissions = np.zeros([2,num_states])
    if twoemits:
        t_emissions2 = np.zeros([2,num_states])
        c_emissions2 = np.zeros([2,num_states])
    for i in range(num_states):
        t_emissions[0,i] = model[1][states[i]][0]
        t_emissions[1,i] = model[1][states[i]][1]
        c_emissions[0,i] = model[2][states[i]][0]
        c_emissions[1,i] = model[2][states[i]][1]
        if twoemits:
            t_emissions2[0,i] = model[1][states[i]][2]
            t_emissions2[1,i] = model[1][states[i]][3]
            c_emissions2[0,i] = model[2][states[i]][2]
            c_emissions2[1,i] = model[2][states[i]][3]
    if twoemits:
        return t_emissions, c_emissions, t_emissions2, c_emissions2
    return t_emissions, c_emissions
