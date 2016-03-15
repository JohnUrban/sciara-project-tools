import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg

def protocol1(late, early, pseudocount=0.1):
    late.median_normalize_data()
    if early:
        early.median_normalize_data()
        late.normalize_to_other(early, pseudocount)
    return late

def protocol2(late, early, bw=10000, pseudocount=0.1):
    late.ksmooth_counts(bw=bw)
    late.median_normalize_data()
    if early:
        early.ksmooth_counts(bw=bw)
        early.median_normalize_data()
        late.normalize_to_other(early, pseudocount)
    return late

def protocol3(late, early, bw=10000, pseudocount=0.1):
    late.median_normalize_data()
    late.ksmooth_counts(bw=bw)
    if early:
        early.median_normalize_data()
        early.ksmooth_counts(bw=bw)
        late.normalize_to_other(early, pseudocount)
    return late

def protocol4(late, early, bw=10000, pseudocount=0.1):
    late.median_normalize_data()
    if early:
        early.median_normalize_data()
        late.normalize_to_other(early, pseudocount)
    late.ksmooth_counts(bw=bw)
    return late

##TODO - allow imputing values locally when a bin is 0 -- i.e. if surrounded by 2 bins with values >0, take average.
## various 0 spots are causing short state changes in the CN hmm.
## perhaps do ksmoothing with something like 10-50kb bandwidth -- then go back to raw signal, and wherever 0 is, substitute with Ksmoothed value.
## 0 spots lead to 2 problems:
##   1. these become fe=1 when both early and late are 0 -- leads to state drops and even bad smoothing
##   2. when not 0 in late, this becomes (count+0.1)/0.1 in fe step (after pseudocount)
##    --> which could mean a really high bin
## best stage to impute would be prior to ANY normalization
## i.e. impute late and early separately --> then median norm --> then FE --> then smooth if want
## If one does protocol 2 or 3... technically imputing is done just by smoothing...
##   only problem is it still shows a significant drop at 0 spots... whereas would be less the case if pre-imputed even before this smoothing

## NOTES: imputing in this way creates its own problems.
## e.g. if late needs to be imputed, but early does not
##  then late will get value close to other values, but early will remain near 0 -- so you get massive spike
##  would want to change same bins in both samples...
##  didnt seem to change state path at all
def normalize(latestage, protocol=1, earlystage=False, pseudo=0.1, bandwidth=2500, quiet=False, impute=False):
    if not quiet:
        newmsg("loading late stage file")
    late = CovBed(latestage)
    if impute:
        if not quiet:
            newmsg("imputing late stage bins with missing data")
        late.impute_zeros(bw=impute)
            
    if earlystage:
        if not quiet:
            newmsg("loading early stage file")
        early = CovBed(earlystage)
        if impute:
            if not quiet:
                newmsg("imputing early stage bins with missing data")
            early.impute_zeros(bw=impute)
    else:
        early = False
    ## todo -- add filtering out 0 contigs option...
    

    if not quiet:
        if earlystage:
            emsg = ' with additional early stage normalization'
        else:
            emsg = ' without additional early stage normalization'
        if protocol == 1:
            newmsg("following normalization protocol 1"+emsg)
        elif protocol == 2:
            newmsg("following normalization protocol 2"+emsg)
        elif protocol == 3:
            newmsg("following normalization protocol 3"+emsg)
        elif protocol == 4:
            newmsg("following normalization protocol 4"+emsg)
    if protocol == 1:
        late = protocol1(late, early, pseudo)
    elif protocol == 2:
        late = protocol2(late, early, bandwidth, pseudo)
    elif protocol == 3:
        late = protocol3(late, early, bandwidth, pseudo)
    elif protocol == 4:
        late = protocol4(late, early, bandwidth, pseudo)
    return late


def run(parser, args):
    ## TODO - just change to 1 argument: --protocol -- with options [1,2,3,4]
    if args.protocol1:
        protocol=1
    elif args.protocol2:
        protocol=2
    elif args.protocol3:
        protocol=3
    elif args.protocol4:
        protocol=4
    late = normalize(latestage=args.latestage, protocol=protocol, earlystage=args.earlystage, pseudo=args.pseudo, bandwidth=args.bandwidth, quiet=args.quiet, impute=args.impute)
    sys.stdout.write( late.get_bdg(late.count, args.collapsed) )
