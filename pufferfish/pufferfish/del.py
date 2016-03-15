## getcov for all -- saved in bincounts files
## corscores
## take in FOFN-like file of 2 columns: 1=stage_integer,2=bincounts_filepaths
##  where stage integer is time points -- e.g. 1,2,3,4,5
##  stage_ints can be used on replicates (i.e. can use same int >1 time)
## python parsers those
##  builds covBed representations
##  builds multiCovBed representation


## TODO:
## getcov needs samtools option to filter mappings based on mapq
## Need to remove chromosomes where >X% of bins have 0 reads in >Y% of samples
##      These will be ignored, and assumed to be issues with assembly (or L)
## Need to iterate:
##  median norm
##  get cor scores
##  get state path
##  take median from 0 states in path
##  go back to median norm step.
##  break iteration when condition_X met
##      --- change in nbins=0state < x
##          -- doesnt necessarily solve problem of too low median causing general ..
##      --- abs(mean_0state_bins - 1) < 1e-2 --- i.e median/median = 1, using good median means mean(X/median) is close to 1. diff between median used to normalize

import sys, os
import datetime
from CovBedClass import *
import cPickle as pickle
##f = CovBed(sys.argv[1])
##
##print f


## OPEN PICKLE INSTEAD OF REDOING ALL CALCS EACH TIME WHILE TESTING
##sys.stderr.write(str(datetime.datetime.now()) +": ..Un-Pickling...\n")
##with open('data.pkl','rb') as pkfile:
##    f = pickle.load(pkfile)
    
##sys.stderr.write(str(datetime.datetime.now()) +": ..initializing...\n")
##f = MultiCovBed(sys.argv[1])
##sys.stderr.write(str(datetime.datetime.now()) +": ..median normalizing...\n")
##f.median_normalize_data()
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing meta-bdg with median norm scores from all files...\n")
##print f
##print f.get_corscor_bdg()



##sys.stderr.write(str(datetime.datetime.now()) +": ..calc corscores...\n")
##f.cor_score()
##sys.stderr.write(str(datetime.datetime.now()) +": ..ksmoothing corscores...\n")
##f.ksmooth_corscores(bw=15000)
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing ksmoothed corscores...\n")
##print f.get_smooth_corscor_bdg()
##sys.stderr.write(str(datetime.datetime.now()) +": ..getting viterbi state paths...\n")
##f.get_cor_states()
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing viterbi state paths...\n")
##print f.get_cor_state_bdg()
##sys.stderr.write(str(datetime.datetime.now()) +": ..analyzing 0-state bins...\n")
##f.analyze_state0_bins()
##print f.state0_medians

##sys.stderr.write(str(datetime.datetime.now()) +": ..calc s50 values...\n")
##f.calc_s50()
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing s50 values...\n")
##print f.get_s50_bdg()

##sys.stderr.write(str(datetime.datetime.now()) +": ..calc a50 values...\n")
##f.calc_a50()
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing a50 values...\n")
##print f.get_a50_bdg()



##print f


## PICKLE f AFTER DOING ALL OPERATIONS SO DONT NEED TO DO THEM OVER AND OVER WHILE TESTING
##sys.stderr.write(str(datetime.datetime.now()) +": ..Pickling...\n")
##with open('data.pkl','wb') as pkfile:
##    pickle.dump(f, pkfile)


sys.stderr.write(str(datetime.datetime.now()) +": ..Un-Pickling...\n")
with open('data.pk2','rb') as pkfile:
    g = pickle.load(pkfile)

sys.stderr.write( str(g.files) +"\n")

g.discretize_cor_values()
for chrom in g.chromosomes:
    print g.dcorscores[chrom][:20]
    print g.corscores[chrom][:20]
    print min(g.dcorscores[chrom]), max(g.dcorscores[chrom])
    print min(g.corscores[chrom]), max(g.corscores[chrom])
    print
quit()
##sys.stderr.write(str(datetime.datetime.now()) +": ..initializing...\n")
##g = MultiCovBed(sys.argv[1])
##sys.stderr.write( str(g.files) +"\n")
##with open('data.pk3','rb') as pkfile:
##    g = pickle.load(pkfile)
##if os.path.exists("data.pk3"):
##    os.remove("data.pk3")
####with open('data.pk3','wb') as pkfile:
##    pickle.dump(g, pkfile)
##sys.stderr.write(str(datetime.datetime.now()) +": ..filtering...\n")
##g.filter_null_contigs(0.4,0.4)
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing filtered contig info...\n")
##print g.get_filtered_contigs_bdg()

##sys.stderr.write(str(datetime.datetime.now()) +": ..finding CN=1...\n")
##g.find_cn1(n_iter=2)

##sys.stderr.write(str(datetime.datetime.now()) +": ..Pickling...\n")
##if os.path.exists("data.pk2"):
##    os.remove("data.pk2")
##with open('data.pk2','wb') as pkfile:
##    pickle.dump(g, pkfile)

sys.stderr.write(str(datetime.datetime.now()) +": ..printing viterbi state paths...\n")
print g.get_cor_state_bdg()

##sys.stderr.write(str(datetime.datetime.now()) +": ..estimating eFDR...\n")
##g.ntest = None
##for i in range(10):
##    ans = g.eFDR(stage=3, n_iter=1)
##    sys.stderr.write(("\n").join([str(e) for e in ans])+"\n")
##sys.stderr.write(str(datetime.datetime.now()) +": ..printing viterbi state paths...\n")
##print g.get_cor_state_bdg()

### thinking through probabilities
## genome size 250-350e6
## 18 puffs
## mean puff size -- assume at most 100 kb mean
## under-rep regions -- if present -- assume only 18 as well
## assume at most, mean size is mean size of puff .. prob smaller though...
## initial prob for -1,0,+1: 0.006, 0.988, 0.006
## transitions:
## there are up to 18 dna puffs and 18 underrep regions -- 36 total.
## haploid autosomes + single X likely ~275 Mb (though asm size ~300-310, and genome with L ~360)
## there is likely on average 7.5 - 8 Mb between centers of puffs+under regions
## --> 275/36 or 300 or 360 / 36 ---
## so when in cn=1 stretch, the prob of seeing a cn+- is 1/8M each
##      and the prob of staying self is 1-x
## since DNA puffs are at least 100 kb long on average, the prob of leaving puff
##  is 1/100k to each state, and 1-x for staying in
##  same transitions can be given to under-rep regions
##
## still to do:
## what is the fraction of reads in amplicons across stages?
## even though Im able to estimate regions that are CN=1
## What fraction of bins are 0-state? I checked --~94.6%
## What fraction of bins are expected to be if puffs are 500k long?
##  --> ~96.5-97.4 (250-350e6 G)
##  if mean size is 100e3: ~99.3-99.5 bins will be 0 state

def fdr(exp_proportion_positive, power, sig_level, ntests):
    exp_proportion_positive = float(exp_proportion_positive)
    TP = ntests*exp_proportion_positive*power
    FN = ntests*exp_proportion_positive*(1-power)
    TN = ntests*(1-exp_proportion_positive)*(1-sig_level)
    FP = ntests*(1-exp_proportion_positive)*sig_level
    FDR = FP/(FP+TP)
    return TP/ntests, FN/ntests, TN/ntests, FP/ntests, FDR

#ntest=607450
#exp prop = 3600.0/607450 = 0.005926413696600543
# power = 0.8
# siglevel = 0.05 --> need 1e-5 for FDR < 0.1%

## can maybe estimate FDR by randomly assigning stages counts
## i.e. just change order of stagelist... how many are positive?
## do that 1000 times... what is average eFDR?
