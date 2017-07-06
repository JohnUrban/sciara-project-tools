

## Take in:
## - 7state path bdg
## - later stage bdg (and early stage bdg)
## - define puff:
##    default:
# state > 1,
# merge regions within 10kb,
# keep products >= 50kb,
# merge kept products that are within 40kb (keep note of states),
# IF LATEST stage, only keep products who have an internal (within middle 80% of coordinates) max state > 2....
### in other words, things that appear to have amplified 2 or more times at last stage
### for earlier stages, use these coordinates to identify anything inside...
import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg
from normalize import protocol1, protocol2, protocol3, protocol4, normalize
from pybedtools import *
from pybedtools.featurefuncs import greater_than
import cStringIO as StringIO


def find_candidate_regions(states, thresh_state=1, merge1=10e3, minwidth=50e3, merge2=40e3, max_state_thresh=2, internal=0.8):
    '''states is a filename pointing to statepath bdg output of cn'''
    '''works on both collapsed and expanded versions'''
    a = BedTool(states)
    return a.filter(lambda x: float(x[3]) > thresh_state).merge(d=int(merge1),c=4,o='max').filter(greater_than, int(minwidth)).merge(d=int(merge2),c=4,o='max').filter(lambda x: float(x[3]) > max_state_thresh)
##    return a.filter(lambda x: float(x[3]) > thresh_state).merge(d=int(merge1),c=4,o='max').filter(lambda x: len(x) > int(minwidth)).merge(d=int(merge2),c=4,o='max').filter(lambda x: float(x[3]) > min_internal_height)
    # could do filter(lambda x: len(x) > int(minwidth)) or filter(greater_than, int(minwidth)) - where greater_than is from pybedtools.featurefuncs 

def bedtool2bedfile(bedtool):
    '''bedtool is a BedTool object
       this returns list as if reading lines from BED file
    '''
    return str(bedtool).strip().split("\n")

def add_len_to_feature(feature):
    '''Example of appending info to BED lines'''
    return feature.append(str(len(feature)))

def add_len_to_bedtool(bedtool):
    '''Example of using "each" method - similar to map() but returns BedTool object and supports chaining; similar to filter'''
    return bedtool.each(add_len_to_feature)

## For each interval in file A, take overlapping sub-interval in file B with highest score in column 4
# for starters: intersectBed -u -a b.bed -b a.bed | mergeBed -i - -c 2,3,4 -o collapse,collapse,collapse
# a.intersect(b, sortout=True).merge(c=(2,3,4), o='collapse,collapse,collapse')
def max_subinterval(feature):
    '''Assumes 6 columns: chr,start,end,collapsed_starts,collapsed_ends,collapsed_scores
       Returns middle-most interval with max score as chr,start,end
       returns all collapsed starts and ends in cols 4 and 5
       col 6 simply has the max score value.
       This fxn is designed to be used with the .each() method for pyBedTools.
       e.g. a.intersect(b, sortout=True).merge(c=(2,3,4), o='collapse,collapse,collapse').each(max_subinterval)
       It will identify the max score in col 6 and the indexes of where it occurs in the list given in col6.
       It will then return the feature with the coordinates of the middle-most max score in the list.
       But will also return in the coordinates of all in cols 4 (starts) and 5 (ends).
       Col6 is returned as the value of the max score.
       '''
    starts = [int(e) for e in feature[3].split(",")]
    ends = [int(e) for e in feature[4].split(",")]
    scores = [float(e) for e in feature[5].split(",")]
    mscore = max(scores)
    nscores = len(scores)
    mi = [i for i in range(nscores) if scores[i] == mscore]
    nmax = len(mi)
    mid_mi = mi[(nmax+1)//2 - 1]
    feature.start = starts[mid_mi]
    feature.end = ends[mid_mi]
    feature[3] = [(",").join([str(starts[i]) for i in mi]) if nmax > 1 else "-"][0]
    feature[4] = [(",").join([str(ends[i]) for i in mi]) if nmax > 1 else "-"][0]
    feature[5] = str(mscore)
    return feature


def summits(a,b):
    '''Assumes 'a' and 'b' BedTool objects.
       Assumes 'a' is a bedGraph (with scores in col4) and has subintervals of 'b' -- and 'b' contains subintervals of 'a'.
       i.e. assumes 'a' will contain the summits and 'b' has the regions you wish to find summits.
       Returns a BedTool.
       Each interval will describe the coordinates of the middle-most max score in a given interval from 'b'.
       But will also return in the coordinates of all summits in the 'b' interval in cols 4 (starts) and 5 (ends).
       Col6 is returned as the value of the max score (summit score) of that region.'''
    return a.intersect(b, sortout=True).merge(c=(2,3,4), o='collapse,collapse,collapse').each(max_subinterval)

def run(parser, args):
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

    if args.regions:
        # find maximums (summits) within regions given
        regions = BedTool(args.regions)
        
    else:
        # find peak regions by algorithm at top, then summits within them
        ## Read in states bedGraph, identify peaks
##        states = CovBed(args.states)
        regions = find_candidate_regions(args.states, thresh_state=1, merge1=10e3, minwidth=50e3, merge2=40e3, max_state_thresh=2, internal=0.8)

    ##Covert CovBed object to BedTool object
    a = BedTool( StringIO.StringIO( late.get_bdg(bdg=late.count, collapsed=True) ) )
    ans = summits(a = a, b = regions)
    print str(ans).strip()

####print "Loading..."
##states = sys.argv[1]
####print "Finding..."
##peakregions = find_candidate_regions(states, thresh_state=1, merge1=10e3, minwidth=50e3, merge2=40e3, min_internal_height=2, internal=0.8)
####sys.stdout.write( str(peakregions) )
####print peak_regions_bed(peakregions).strip()
##
