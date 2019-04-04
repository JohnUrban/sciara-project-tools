import os
from collections import defaultdict
import numpy as np
##np.seterr(divide='raise', invalid='raise')
np.seterr(divide='ignore', invalid='raise')
##np.seterr(divide='ignore', invalid='ignore')
import rpy2.robjects as robjects
ksmooth = robjects.r['ksmooth']
kmeans = robjects.r['kmeans']
intvec = robjects.IntVector
fltvec = robjects.FloatVector
matrixr = robjects.r.matrix
r = robjects.r
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as stap
from puffR import *
puffR = stap(puffRstring, 'puffR')
import sys, datetime
from scipy.stats import spearmanr

def s50(counts, stages, x=[50]):
        """
        counts,stages,x lists
        Returns sX for all x for a list of numbers "counts".
        Default: 50
        Assumes all values in list x are between 0 and 100.
        Interpretation: Returns stage at which >= 50% of reads is hit.
        """
        n = len(counts)
        x_to_stage = {e:0 for e in x}
        count_sum = sum(counts)
        total = 0
        i=0
        for e in sorted(x):
                target = count_sum*e/100.0
                while total < target and counts:
                        total += counts[i]
                        lastcount = counts[i]
                        laststage = stages[i]
                        i+=1
                try:
                    x_to_stage[e] = laststage
                except UnboundLocalError:
                    x_to_stage[e] = "."
        return x_to_stage


TMP_DIR = ".pufferfish_tmp_dir"
TMP_DIR = TMP_DIR[1:]

class CovBed(object):
    def __init__(self, covbedfile, count_only=False, replace=False, replace_with='0', replace_this='.'):
        ## "replace" means if you see the "replace_this" character in the count column, make it "replace_with"
        ## Made to deal with "." --> 0 by default when replace used.
        self.fopen = False
        self.connection = None
        self.file = covbedfile
        self.start = {}
        self.end = {}
        self.count = {}
        self.chromosomes = set([])
        self.median = None
        self.mean = None
        self.sd = None
        self.count_only=count_only ## When False, start/end dicts are initialized, but remain empty: useful when comparing 2 bedgraphs of identical coords. See also MultiCovBed (though that currently requires a "stage file")
        self._extract_data(replace, replace_with, replace_this)
        
    def open(self):
        if self.fopen:
            self.close()        
        self.connection = open(self.file, 'r')
        self.fopen = True


    def close(self):
        self.connection.close()

    def _add_chromosome(self, chrom):
        self.chromosomes.add(chrom)
        self.start[chrom] = []
        self.end[chrom] = []
        self.count[chrom] = []

    def _update_data(self, chrom, start, end, count):
        if chrom not in self.chromosomes:
            self._add_chromosome(chrom)
        if not self.count_only: ## This allows the start/end dicts to be initialized, but remain empty
                self.start[chrom].append(start)
                self.end[chrom].append(end)
        self.count[chrom].append(count)

    def _finalize_data(self):
        ## convert set to list
        self.chromosomes = sorted(list(self.chromosomes))
        #convert lists to np arrays
        for chrom in self.chromosomes:
            self.start[chrom] = np.array(self.start[chrom])
            self.end[chrom] = np.array(self.end[chrom])
            self.count[chrom] = np.array(self.count[chrom])

                
    def _extract_data(self, replace=False, replace_with='0', replace_this='.'):
        self.open()
        for line in self.connection:
            chrom, start, end, count = line.strip().split()
            if replace and count == replace_this:
                count = replace_with
            self._update_data(chrom, int(start), int(end), float(count))
            ### JAN 9, 2018 -- I changed above int(float(count)) to float(count)
            ###         At this point, I don't know what it might break...
            ##          But now that I am using this more generally for signal data - not just counts - I need it as float
        self._finalize_data()
        self.close()


    def get_mean(self):
        if self.mean is None:
            counts = np.concatenate(self.count.values())
            self.mean = float(np.mean(counts))
            self.sd = float(np.std(counts))
        return self.mean

    def get_sd(self):
        if self.sd is None:
            counts = np.concatenate(self.count.values())
            self.mean = float(np.mean(counts))
            self.sd = float(np.std(counts,ddof=1))
        return self.sd

    def _get_median(self):
        counts = np.concatenate(self.count.values())
        self.median = float(np.median(counts))
##        counts = []
##        for chrom in self.chromosomes:
##            counts.append(self.count[chrom])
##        self.median = float(np.median(counts))
        
    def get_median(self):
        if self.median is None:
            self._get_median()
        return self.median

    def median_normalize_x(self, x):
        #x is np.array
        return x/self.get_median()

    def median_normalize_data(self):
        for chrom in self.chromosomes:
            self.count[chrom] = self.median_normalize_x(self.count[chrom])
            
    def expanded_bdg(self, bdg):
        ##bdg is just what should be in the 4th column
        string = ''
        for chrom in self.chromosomes:
            for i in range(len(self.start[chrom])):
                string += ('\t').join([chrom, str(self.start[chrom][i]),  str(self.end[chrom][i]),  str(bdg[chrom][i])]) + "\n"
        return string

    def expanded_bdg_two_cols(self, bdg1, bdg2):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(bdg1[chrom][i]), str(bdg2[chrom][i])]) + "\n"
        return string
    
    def collapsed_bdg(self, bdg):
        ##bdg is just what should be in the 4th column
        string = ''
        for chrom in self.chromosomes:
            if len(self.start[chrom]) > 1:
                #init
                start = self.start[chrom][0]
                value = bdg[chrom][0]
                for i in range(1, len(self.start[chrom]) ):
                    if bdg[chrom][i] != value:
                        string += ('\t').join([chrom, str(start), str(self.end[chrom][i-1]), str(value)]) + "\n"
                        start = self.start[chrom][i]
                        value = bdg[chrom][i]
                ##finish chrom
                string += ('\t').join([chrom, str(start), str(self.end[chrom][i]), str(value)]) + "\n"
            else: #only 1 bin (very tiny contig)
                string += ('\t').join([chrom, str(self.end[chrom][0]), str(self.end[chrom][0]), str(bdg[chrom][0])]) + "\n"
        return string


    def get_bdg(self, bdg, collapsed=False):
        if not collapsed:
            return self.expanded_bdg(bdg)
        else:
            return self.collapsed_bdg(bdg)

    def filtered_bdg(self, relation = ">", value = 0, bdg=None):
        ##bdg is just what should be in the 4th column
        ## for this might typically be self.count
        if bdg is None:
            bdg = self.count
            string = ''
        if relation == "gt":
            keep = lambda x: x > value
        elif relation == "ge":
            keep = lambda x: x >= value
        elif relation == "lt":
            keep = lambda x: x < value
        elif relation == "le":
            keep = lambda x: x <= value
        elif relation == "eq":
            keep = lambda x: x == value
        elif relation == "ne":
            keep = lambda x: x != value
        for chrom in self.chromosomes:
            for i in range(len(self.start[chrom])):
                if keep(bdg[chrom][i]):
                    string += ('\t').join([chrom, str(self.start[chrom][i]),  str(self.end[chrom][i]),  str(bdg[chrom][i])]) + "\n"
        return string
    
    def __str__(self):
        return self.get_bdg(self.count)
    
    def get_chromosomes(self):
        return self.chromosomes

    def get_start_dict(self):
        return self.start
    def get_end_dict(self):
        return self.end
    def get_count_dict(self):
        return self.count

    def ksmooth_counts(self, bw=10000):
        for chrom in self.chromosomes:
            x = self.start[chrom]
            y = self.count[chrom]
            k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
            self.count[chrom] = np.array(k[1])
##            sys.stderr.write(chrom+"\n")
##            sys.stderr.write(str(x)+"\n")
##            sys.stderr.write(str(y)+"\n")
##            sys.stderr.write(str(np.array(k[1]))+"\n\n\n")
####            quit()

    def normalize_to_other(self, other, pseudocount=0.01):
        #other is another CovBed object with same bins from same genome
        for chrom in self.chromosomes:
            self.count[chrom] = (np.array(self.count[chrom])+pseudocount)/(np.array(other.count[chrom])+pseudocount)

    def impute_zeros(self, bw):
        '''When requiring mapq to be stringent, it leaves stretches of 0 that could benefit from being imputed.
        The 0 score often causes a state change in HMMs and can also lead to very inflated scores in FE after pseudocount added (if numerator is non-zero - e.g. 10/0.1 = 100)
        So this smooths the counts... and only uses the resulting smoothed values to substitute 0s.
        This means in very long 0 regions (e.g. contigs with no coverage), the score will remain 0 as desired.'''
        for chrom in self.chromosomes:
            x = self.start[chrom]
            y = self.count[chrom]
            k = puffR.impute_zeros(x = fltvec(x), y = fltvec(y), bw = bw)
            self.count[chrom] = np.array(k)
        
        










class MultiCovBed(object):
    ## assumes all covbeds in list have same exact elements in 1st 3 colums and are all identically ordered
    ## this should be true if all are outputs from getcov using the same genome file (just different bams)
    def __init__(self, stagefile):
        # 'stagefile' is a FOFN-like file with 2 columns: 1=stage_integer,2=bincounts_filepaths
        ##  where stage integer is time points -- e.g. 1,2,3,4,5
        ##  stage_ints can be used on replicates (i.e. can use same int >1 time)
        self.stagefile = stagefile
        self.nfiles = 0
        self.nbins = 0
        self.covbeds = {}
        self.stages = {}
        self.stagelist = []
        self.files = {}
        self.count = {}
        self.median = {}
        self.start = {}
        self.end = {}
        self.chromosomes = set([])
        self._parse_stagefile()
        self._extract_data()
        self.corscores = {}
        self.rksmooth = None
        self.smooth_corscores = {}
        self.cor_states = {}
        self.s50 = {}
        self.a50 = {}
        self.filtered = {'start':{}, 'end':{}, 'count':{k:{} for k in range(self.nfiles)}}
        self.filteredchromosomes = []
        self.ntest = None


    ######################
    ##  INITIALIZATION
    ######################
    def _parse_stagefile(self):
        i = 0
        with open(self.stagefile) as f:
            for line in f:
                stage, fname = line.strip().split()
                self.stages[i] = int(stage)
                self.files[i] = fname
                self.stagelist.append(int(stage))
                i+=1
        self.nfiles = len(self.files.keys())
            
        
    def _extract_data(self):
        for i in sorted(self.files.keys()):
            self.add_covbed(findex = i)

    def add_covbed(self, findex):
        covbed = CovBed(self.files[findex])
        if not self.chromosomes:
            self._initialize(covbed)
        self.count[findex] = covbed.get_count_dict()
            

    def _initialize(self, covbed):
        self.start = covbed.get_start_dict()
        self.end = covbed.get_end_dict()
        self.chromosomes = sorted(list(covbed.get_chromosomes()))


    ######################
    ## Operations
    ######################
    def _get_median(self, findex):
        counts = np.concatenate(self.count[findex].values())
        self.median[findex] = float(np.median(counts))
        
    def get_median(self, findex, refresh=False):
        if refresh:
            self._get_median(findex)
        try:
            return self.median[findex]
        except KeyError as e:
            self._get_median(findex)
            return self.median[findex]

    def _refresh_medians(self):
        ## will re-calculate medians every time called
        for findex in range(self.nfiles):
            self._get_median(findex)

    def normalize(self, x, denom):
        #x is np.array
        # denom is float
        return x/denom

    def normalize_findex_by_x(self,findex, x):
        for chrom in self.chromosomes:
            self.count[findex][chrom] = self.normalize(self.count[findex][chrom], x)

    def nomalize_data_by_xdict(self,xdict):
        for findex in range(self.nfiles):
            self.normalize_findex_by_x(findex, xdict[findex])

    def median_normalize_findex(self,findex, refresh=False):
        self.normalize_findex_by_x(findex, x = self.get_median(findex, refresh))

    def median_normalize_data(self, refresh=False):
        for findex in range(self.nfiles):
            self.median_normalize_findex(findex, refresh)

    def ksmooth_counts(self, bw=10000):
        for chrom in self.chromosomes:
            for findex in range(self.nfiles):
                x = self.start[chrom]
                y = self.count[findex][chrom]
                k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
                self.count[findex][chrom] = np.array(k[1])

    def expanded_bdg(self, bdg):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(bdg[chrom][i])]) + "\n"
        return string

    def expanded_bdg_two_cols(self, bdg1, bdg2):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(bdg1[chrom][i]), str(bdg2[chrom][i])]) + "\n"
        return string
    
    def collapsed_bdg(self, bdg):
        string = ''
        for chrom in self.chromosomes:
            #init
            start = self.start[chrom][0]
            value = bdg[chrom][0]
            for i in range(1, len(self.start[chrom]) ):
                if bdg[chrom][i] != value:
                    string += ('\t').join([chrom, str(start), str(self.end[chrom][i-1]), str(value)]) + "\n"
                    start = self.start[chrom][i]
                    value = bdg[chrom][i]
            ##finish chrom
            string += ('\t').join([chrom, str(start), str(self.end[chrom][i]), str(value)]) + "\n"
        return string

    def get_bdg(self, bdg, collapsed=False):
        if not collapsed:
            return self.expanded_bdg(bdg)
        else:
            return self.collapsed_bdg(bdg)
        
    def find_slopes(self, stagelist=''):
        if not stagelist:
            stagelist = self.stagelist
        self.slopes = {}
        for chrom in self.chromosomes:
            self.slopes[chrom] = []
            for i in range(len(self.start[chrom])):
                counts = [self.count[j][chrom][i] for j in range(self.nfiles)]
                try:
                    slope = np.polyfit(x = stagelist, y = counts, deg = 1)[0]
                except FloatingPointError:
                    slope = 0
                self.slopes[chrom].append(slope)

    def get_slope_bdg(self, collapsed=False):
        ##assumes self.slopes already present
        return self.get_bdg(self.slopes, collapsed)


                
    def cor_score(self, stagelist=''):
        if not stagelist:
            stagelist = self.stagelist
        for chrom in self.chromosomes:
            self.corscores[chrom] = []
            for i in range(len(self.start[chrom])):
                counts = [self.count[j][chrom][i] for j in range(self.nfiles)]
                try:
##                    score = np.nan_to_num( np.corrcoef(x = [stagelist, counts])[0,1] )
                    score = np.corrcoef(x = [stagelist, counts])[0,1] 
                except FloatingPointError:
                    score = 0
                self.corscores[chrom].append(score)


    def get_corscore_bdg(self, collapsed=False):
        if not self.corscores:
            self.cor_score()
        return self.get_bdg(self.corscores, collapsed)

    def ksmooth_corscores(self, bw=10000):
        if not self.corscores:
            self.cor_score()
        for chrom in self.chromosomes:
            x = self.start[chrom]
            y = self.corscores[chrom]
            k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
            self.smooth_corscores[chrom] = np.array(k[1])


    def get_smooth_corscore_bdg(self, collapsed=False):
        return self.get_bdg(self.smooth_corscores, collapsed)


    def get_cor_states(self, smoothed=False, emodel="normal"):
        if not self.corscores:
            self.cor_score()
        if smoothed and not self.smooth_corscores:
            sys.stderr.write("Smoothing cor scores with default bandwidth")
            self.ksmooth_corscores()
        if smoothed:
            scores = self.smooth_corscores
        else:
            scores = self.corscores
        for chrom in self.chromosomes:
##            sys.stderr.write( chrom + "\n" )
            v = puffR.viterbi_puff(emissions = puffR.emissions, transitions = puffR.transitions, initial = puffR.initial, states = intvec([1,2,3]), emitted_data = fltvec(scores[chrom]), emodel = emodel, logprobs=False)
##            f = puffR.forward_puff(emissions = puffR.emissions, transitions = puffR.transitions, initial = puffR.initial, states = intvec([-1,0,1]), emitted_data = fltvec(scores[chrom]), emodel = emodel, logprobs=False)
##            b = puffR.backward_puff(emissions = puffR.emissions, transitions = puffR.transitions, initial = puffR.initial, states = intvec([-1,0,1]), emitted_data = fltvec(scores[chrom]), emodel = emodel, logprobs=False)
            
            self.cor_states[chrom] = np.array(list(v[0]))


    def get_cor_state_bdg(self, collapsed=False):
        return self.get_bdg(self.cor_states, collapsed)

    

    
    def calc_s50(self, x=[50], stagelist=''):
        if not stagelist:
            stagelist = self.stagelist
        for chrom in self.chromosomes:
            self.s50[chrom] = []
            for i in range(len(self.start[chrom])):
                counts = [self.count[j][chrom][i] for j in range(self.nfiles)]
                ans = s50(counts,stagelist,x=x)
                ans = ans[x[0]]
                self.s50[chrom].append(ans)

                
    def get_s50_bdg(self):
        if not self.s50:
            self.calc_s50()
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                if self.s50[chrom][i] != ".":
                    string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(self.s50[chrom][i])]) + "\n"
        return string

    def analyze_state0_bins(self):
        ## assumes median normalized (or otherwise)
        if not self.corscores:
            self.cor_score()
        if not self.cor_states:
            self.get_cor_states()
        counts = {k:[] for k in range(self.nfiles)}
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                if self.cor_states[chrom][i] == 2:#states are 1,2,3 for -1,0,1
                    for j in range(self.nfiles):
                        counts[j].append( self.count[j][chrom][i] )
        ## mean would be skewed above 1 since even if the count of numbers < 1 equals the count of numbers > 1, the mean will be >1.
        total_median = [np.median(counts.values())]
        medians = {findex:np.median(counts[findex]) for findex in range(self.nfiles)}
        self.state0_medians = {0:medians,1:total_median}

    def _need_more_cycles(self):
        for m in self.state0_medians[0].values():
            if m != 1:
                return True
        return False

    def _normalize_data_by_state0_median(self):
        for findex in range(self.nfiles):
            if self.state0_medians[0][findex] != 1:
                self.normalize_findex_by_x(findex, self.state0_medians[0][findex])

    def find_cn1(self, n_iter=10, max_empty_bin_pct=0.4, max_offending_samples_pct=0.4, verbose=True):
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - getting medians before and after filtering..\n")
        self._refresh_medians() ## get medians for first time
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - medians before filtering were %s...\n" % (str(self.median)))
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - filtering..\n")
        self.filter_null_contigs(max_empty_bin_pct, max_offending_samples_pct) ## filter bad contigs
        self._refresh_medians() ## get medians after filtration
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - medians after filtering were %s...\n" % (str(self.median)))
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - median normalizing..\n")
        self.median_normalize_data() # median norm
        self._refresh_medians() ## get medians after normalization
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - medians (of all bins) after normalizing should be 1 and were %s...\n" % (str(self.median)))
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - getting bin correlations \n")
        self.cor_score() # give correlation score to each bin
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - getting state path \n")
        self.get_cor_states() # get state path through bins
        self.analyze_state0_bins() # get median of 0-correlation (most likely cn=1) bins
        for i in range(n_iter):
            if verbose:
                sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 iter %d starting 0-corr medians: %s...\n" % (i+1, str(self.state0_medians)))
            if self._need_more_cycles(): # if 0-corr median from any sample != 1, re-normalize by 0-corr median, est new correlations, get updated statepath, get new 0corr medians
                self._normalize_data_by_state0_median()
                self.cor_score()
                self.get_cor_states()
                self.analyze_state0_bins()
            else:
                if verbose:
                    sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 iter %d canceled - 0-corr medians all = 1...\n" % (i+1))
                break
        
    def filter_null_contigs(self, max_empty_bin_pct=0.4, max_offending_samples_pct=0.4):
        ##default: if >40% of samples have >40% empty bins on a chrom, remove it
        sample_thres = max_offending_samples_pct*self.nfiles
        for chrom in self.chromosomes:
            nbins = len(self.start[chrom])
            bin_thres = max_empty_bin_pct*nbins
            n_bad_samp = 0
            for findex in range(self.nfiles):
                n_bad_bins = len(self.count[findex][chrom][self.count[findex][chrom] == 0])
                if n_bad_bins > bin_thres:
                    n_bad_samp += 1
            if n_bad_samp > sample_thres:
                self.filtered['start'][chrom] = self.start.pop(chrom)
                self.filtered['end'][chrom] = self.end.pop(chrom)
                for findex in range(self.nfiles):
                    counts = (self.count[findex]).pop(chrom)
                    self.filtered['count'][findex][chrom] = counts
                self.filteredchromosomes.append(chrom)
                self.chromosomes.remove(chrom)


            
    def pct_state(self, state=3):
        total = 0
        nstate = 0
        for chrom in self.chromosomes:
            total += len(self.cor_states[chrom])
            nstate += sum(self.cor_states[chrom][self.cor_states[chrom] == state])
        return 100.0*nstate/total

    def n_state(self, state=3):
        nstate = 0
        for chrom in self.chromosomes:
            nstate += sum(self.cor_states[chrom][self.cor_states[chrom] == state])
        return nstate
    
    def eFDR1(self, stage=3, n_iter=10):
        ## shuffles stages -- not best way since many permutations uphold correlations...
        ## assumes find_cn1 already run
        if self.ntest is None:
            self.ntest = self.n_state(3)
        controls = []
        for i in range(n_iter):
            stagelist = self.stagelist[:] #clone it
            np.random.shuffle(stagelist)
            sys.stderr.write(str(stagelist)+"\n") ## PRINT TEMP
            self.cor_score(stagelist)
            self.get_cor_states()
            controls.append(self.n_state(3))
        return self.ntest, controls, 100.0*np.array(controls)/self.ntest

    def eFDR2(self):
        ## 
        pass

    def pval(self):
        # for each bin, shuffle scores, take cor, store -- get 1000 of these, p ~ n_gt_cor/N
        # BH correct p-values - only keep bins with q < 0.1
        ## OR -- take state+ bins, combine counts in all for each stage, do this
        ## that would be less tests and it would not treat all bins independently - would treat regions independently
        ## hmmm.... maybe not worth doing this at all...
        pass

    def discretize_cor_values(self):
        self.dcorscores = {}
        for chrom in self.chromosomes:
            self.dcorscores[chrom] = map(round, np.array(self.corscores[chrom])*4+5) #9-sided dice emission symbols 1-9 (for Rindexing) where 1-4 neg, 5=0, 6-9pos
    ######################
    ## Printing etc
    ######################
    def get_counts_bdg(self):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i])] + [str(self.count[j][chrom][i]) for j in range(self.nfiles)]) + "\n"
        return string

    def get_filtered_contigs_bdg(self):
        string = ''
        for chrom in self.filteredchromosomes:
            for i in range(len( self.filtered['start'][chrom] )):
                string += ('\t').join([chrom, str(self.filtered['start'][chrom][i]), str(self.filtered['end'][chrom][i])] + [str(self.filtered['count'][j][chrom][i]) for j in range(self.nfiles)]) + "\n"
        return string
    
    def __str__(self):
        return self.get_counts_bdg()

    
            









## DOES NOT BEHAVE IN USEFUL WAY - i.e. results not more useful than s50 (probably less so).
##    def calc_a50(self,x=[50]):
##        ## ASSUMES COUNTS ARE MEDIAN NORMALIZED
##        ## for all cn=1 areas, s50=3 -- i.e. 50% of normalized read counts seen halfway through
##        ## a50 is trying to mask cn=1 to highlight when 50% of amplification is done
##        ## it does this by subtracting 1 from the median normalized counts (and taking max of that or 0 to avoid negatives)
##        ##  --> and later ignoring anything that is "."
##        ## this is experimental - not totally sure it will give what I want
##        for chrom in self.chromosomes:
##            self.a50[chrom] = []
##            for i in range(len(self.start[chrom])):
##                counts = [max([self.count[j][chrom][i]-1,0]) for j in range(self.nfiles)]
##                ans = s50(counts,self.stagelist,x=x)
##                ans = ans[x[0]]
##                self.a50[chrom].append(ans)
##                
##    def get_a50_bdg(self):
##        if not self.a50:
##            self.calc_a50()
##        string = ''
##        for chrom in self.chromosomes:
##            for i in range(len( self.start[chrom] )):
##                if self.a50[chrom][i] != ".":
##                    string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(self.a50[chrom][i])]) + "\n"
####                if i > 20: break
####            if i > 20: break
##        return string
