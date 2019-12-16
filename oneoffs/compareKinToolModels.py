#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict, Counter
#from Bio import SeqIO
#from numpy.random import choice
import scipy
from scipy.stats import chisquare

parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    Read in a test and background kinTools kmer model,
    Output a report on each grouping.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--test', '-t', type=str, required=True, help='''Path to KinTools model to test.''')
parser.add_argument('--background', '-b', type=str, required=True, help='''Path to KinTools background model to use.''')
parser.add_argument('--klist', '-k', type=str, default='1,2,3,4,5,6,7', help='''Comma-separated list of kmer-sizes (integers) to build model for. Default: 1,2,3,4,5,6,7 ''')

parser.add_argument('--pseudo', '-c', type=int, default=0, help='''Add a pseudo count. Default: 0. ''')

args = parser.parse_args()






class MidBaseModel(object):
    def __init__(self, fh, kmers=[1,2,3,4,5,6,7]):
        self.fh = fh
        self.model = {k:{i:defaultdict(int) for i in range(k)} for k in kmers}
        self.__initialize_model()
        self.kmers = kmers
        
    def __initialize_model(self):
        with open(self.fh) as f:
            for line in f:
                self.__add(line)

                
    def __add(self, line):
        '''seq is string of upper-case'''
        k, i, seq, count = self.__process(line)
        try:
            self.model[k][i][seq] = count
        except:
            pass ## Fail silently on k values not anticipated...
        

    def __process(self, line):
        k, i, seq, count = line.strip().split('\t')
        k = int(k)
        i = int(i)-1
        count = int(count)
        return k, i, seq, count


    def __str__(self):
        #return str(self.model)
        out = ''
        gate = set('ACGT') ## does not allow kmers with N or any other non-canonical
        for k in self.kmers:
            for i in sorted(self.model[k].keys()):
                for seq in sorted(self.model[k][i].keys()):
                    chars = set(seq)
                    if not chars.difference(gate):
                        ## Add 1 to i so it is 1-based
                        out += '\t'.join([str(e) for e in [k, i+1, seq, self.model[k][i][seq]]]) + '\n'
        return out.strip()

    def update_model(self, other):
        for k in other.model.keys():
            if k not in self.model.keys():
                self.model[k] = defaultdict(int)
            for i in other.model[k].keys():
                if i not in self.model[k].keys():
                    self.model[k][i] = 0
                    
    def add_pseudo(self, pseudo=1):
        for k in self.kmers:
            for i in self.model[k].keys():
                for seq in self.model[k][i].keys():
                    self.model[k][i][seq] += pseudo

class CompareModels(object):
    def __init__(self, test_fh, bg_fh, kmers=[1,2,3,4,5,6,7], pseudo=0):
        self.kmers = kmers
        self.test =  MidBaseModel(fh = args.test, kmers=kmers)
        self.bg = MidBaseModel(fh = args.background, kmers=kmers)
        self.test.update_model(self.bg)
        self.bg.update_model(self.test)
        self.pseudo = pseudo
        ## Pseudo will only be used for "FE" calculation in __str__ and only if a 0 count found in denom
        #if pseudo > 0:
        #    self.test.add_pseudo(1)
        #    self.bg.add_pseudo(1)
        
    def __str__(self):
        out = ''
        for k in self.kmers:
            for i in sorted(self.test.model[k].keys()):
                keys = sorted(self.test.model[k][i].keys())
                #test_sum = float(sum(self.test.model[k][i].values()))
                #bg_sum = float(sum(self.bg.model[k][i].values()))
                obs_test_counts = scipy.array([self.test.model[k][i][seq] for seq in keys])
                bg_counts = scipy.array([self.bg.model[k][i][seq] for seq in keys])
                test_sum = float(obs_test_counts.sum())
                bg_sum = float(bg_counts.sum())
                test_probs = obs_test_counts / test_sum
                bg_probs = bg_counts / bg_sum
                exp_test_counts = test_sum * bg_probs
                chisquare_stat, chisquare_p = chisquare(obs_test_counts, f_exp=exp_test_counts)
                residuals = (obs_test_counts - exp_test_counts) / scipy.sqrt(exp_test_counts)
                if len(keys) > 2:
                    mu = scipy.mean(residuals)
                    sd = scipy.std(residuals, ddof=1)
                    med = scipy.median(residuals)
                    absdiffs = scipy.absolute(residuals - med)
                    mad = scipy.median(absdiffs)
                    if mad == 0:
                        mad = scipy.mean(absdiffs)
                        if mad == 0: #still
                            mad = sd/1.4862
                    z = (residuals - mu)/sd
                    rz = (residuals - med)/mad
                else:
                    z = "."
                    rz = "."
                for j in range(len(keys)):
                    seq = keys[j]
                    test_count = obs_test_counts[j]
                    test_prob = test_probs[j]
                    #bg_count = bg_counts[j]
                    bg_prob = bg_probs[j]
                    exp_count = exp_test_counts[j]
                    res = residuals[j]
                    if exp_count == 0: ## PSEUDO USED HERE ONLY WHEN NEEDED
                        fe = test_count / float(exp_count + self.pseudo)
                    else:
                        fe = test_count / float(exp_count)
                        # I do not include bg count
                    out += '\t'.join([str(e) for e in [k, i+1, seq, test_prob, bg_prob, test_count, exp_count, res, fe, chisquare_p, chisquare_stat, z[j], rz[j]] ]) + '\n'
        return out.strip()

    
#### EXECUTE

try:
    ## Initialize models
    testing = CompareModels(test_fh = args.test, bg_fh = args.background, kmers=[int(e) for e in args.klist.split(',')], pseudo=args.pseudo)

    ## RETURN MODEL
    print testing

except IOError:
    pass ## Fail silently here (e.g. when piped into less or head)




quit()

