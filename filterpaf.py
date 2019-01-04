#!/usr/bin/env python2.7
import sys, pandas
import argparse
import numpy as np
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

 Return PAF lines given rules.


    Note: I've only begun tinkering with this script.

    For sorting, assuming a normal PAF, this unix command will give the same results (or should and has in some tests):
    sort -k1,1 -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n

    For bubbles/redundant contigs, I'd recommend using awk for now.
    Something similar to:
    # create a modified PAF describing putative bubbles/redundant contigs (RCs).
    awk '$1!=$6 && (($4-$3)/$2 >= 0.6 || ($9-$8)/$7 >= 0.6) {OFS="\t"; A=($4-$3)/$2; B=($9-$8)/$7; C=$10/$11; if (A >= 0.6 && B < 0.6) print "query",A,B,C,$0; else if (A < 0.6 && B >= 0.6) print "orig_target",B,A,C,$6,$7,$8,$9,$5,$1,$2,$3,$4,$10,$11,$12,$13,$14,$15,$16; else if (A>=0.6 && B>=0.6 && A>=B) print "both_q",A,B,C,$0; else if (A>=0.6 && B>=0.6 && B>=A) print "both_t",B,A,C,$6,$7,$8,$9,$5,$1,$2,$3,$4,$10,$11,$12,$13,$14,$15,$16; else print  "none",A,B,C,$0}' bubble-search-asm10.paf > filtered-bubble-search-asm10.paf

    # define a set of bubbles/RCs using cutoffs on the putative bubble set
    awk '(($2>=0.8 && $4>=0.8)) && $6/$11 < 0.8' filtered-bubble-search-asm10.paf 


    # Merging -- if bubbles can be query or target, will want to merge twice: once on the normal PAF, once on a re-structured PAF using awk:
    awk 'OFS="\t" {print $6,$7,$8,$9,$5,$1,$2,$3,$4,$10,$11,$12}' normal.paf > restructured.paf
    sort -k1,1 -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n normal.paf > sorted-normal.paf
    sort -k1,1 -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n restructured.paf > sorted-restructured.paf
    filterpaf.py --merge -MR 1e3,1e3 -i sorted-normal.paf > merged-sorted-normal.paf
    filterpaf.py --merge -MR 1e3,1e3 -i sorted-restructured.paf > merged-sorted-restructured.paf
    Then call bubbles on each.

    Can also use awk for SAME operations:
    e.g. for -SR 0.6,1.66,0.75,0.75,0.5
    awk '$2/$7 >= 0.6 && $2/$7 <= 1.66 && $10/$2 >=0.75 && $10/$7 >= 0.75 && $10/$11 >= 0.5' aln.paf > same.paf
    Or if you want to merge first, do it w/ filterpaf.py then
    awk '$2/$7 >= 0.6 && $2/$7 <= 1.66 && $10/$2 >=0.75 && $10/$7 >= 0.75 && $10/$11 >= 0.5' merge.paf > same.paf 

    One application of --same is:
    In an initial asm, Canu or Falcon might give bubble tigs -- or maybe you perform a less-stringent filtering in the last step to yield more tigs.
    You can mark these tigs for removal with something in their name, but use them as part of downstream steps, such as evals, polishing, scaffolding etc.
    Ideally, you have kept track of these "rmtigs" throughout all these steps, many of which change the names or require simpler names.
    If not, you can map the original rmtigs back to the final assembly and use --same to determine which tigs in the final asm descend fromt the rmtigs.
    If scaffolding was done, rmtigs might have been combined together. To search for these as well as those above, one can also use awk in combination w/ BEDtools.
    (i) Merge with filterpaf (or an awk/bedtools version of it)
    (ii) awk and bedtools to collapse queries that align on the same target within some distance before "same" calling:
     awk 'OFS="\t" {print $6,$8,$9,$7,$1,$3,$4,$2, $5, $10,$11,$12}' mergedel | sortBed -i - | mergeBed -d 1000 -i - -c 4,5,6,7,8,9,10,11,12 -o distinct,collapse,min,sum,sum,collapse,sum,sum,mean | awk '{OFS="\t"; print $5,$8,$6,$7,$9,$1,$4,$2,$3,$10,$11,$12}' >targetmerged.paf
    (iii) same-calling
    awk '$2/$7 >= 0.6 && $2/$7 <= 1.66 && $10/$2 >=0.75 && $10/$7 >= 0.75 && $10/$11 >= 0.5' targetmerged.paf > same.paf
    (iv) collect target names for re-marking as rmtigs
    (v) You may also want to note where rmtigs were integrated inside large contigs during a scaffolding step (as this may be a mis-assembly)

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--paf', '-i',
                   type=str, required=True,
                   help='''Path to input PAF file.''')

parser.add_argument('--same', '-S',
                   action='store_true', default=False,
                   help='''Return lines where there was a "self-self" alignment as defined by some rules.
                        NOTE: this does not look at sequence names for the self-self definition.
                        It looks at sequence lengths and alignments.''')


parser.add_argument('--same_rules', '-SR',
                   type=str, default='1,1,1,1,1',
                   help='''Define "self-self" alignment with these rules.
                        Provide comma-sep list of: minLenRatio, maxLenRatio, minPctQueryMatch, minPctTargetMatch, minPctAlnMatch.
                        Default: 1,1,1,1,1.
                        Useful: 0.6,1.66,0.75,0.75,0.5.
                        It is treated as: (LenRatio within minLenRatio-maxLenRatio) AND (minPctQueryMatch AND minPctTargetMatch AND minPctAlnMatch)''')


parser.add_argument('--bubble', '-B',
                   action='store_true', default=False,
                   help='''Return lines where there the query was a "bubble"/redundant contig alignment as defined by some rules.
                        ''')


parser.add_argument('--bubble_rules', '-BR',
                   type=str, default='0.8,0.8,1',
                   help='''Define "bubble" alignment with these rules.
                        Provide comma-sep list of: minPctQueryInTarget,minPctAlnMatch,maxPctOfTarget.
                        Default: 0.8,0.8,1.
                        That is: at least 80 pct of query must be aligned in target with at least 80 pct matches in the alignment, and the query can be no longer than the target.
                        ''')

##parser.add_argument('--reverse_bubble', '-RB',
##                   action='store_true', default=False,
##                   help='''Treat query and target the opposite of default.
##                        ''')
##
##parser.add_argument('--any_bubble', '-AB',
##                   action='store_true', default=False,
##                   help='''Look for a bubble in both directions for each PAF record - i.e. default and reverse_bubble.
##                        ''')

parser.add_argument('--same_seq_allowed', '-SS',
                   action='store_true', default=False,
                   help='''By default, when searching for bubbles, alignments can not be self-self as determined by name.
                        In some rare cases, different sequences might have the same name or the user might want to turn off this behavior for their own ends.
                        This flag allows bubbles to be discovered from self-self alignments.
                        ''')

parser.add_argument('--merge', '-M',
                   action='store_true', default=False,
                   help='''Return lines where there the query was a "merge" alignment as defined by some rules.
                        ''')


parser.add_argument('--merge_rules', '-MR',
                   type=str, default='1e9,1e9',
                   help='''Define "merge" alignments with these rules.
                        Provide comma-sep list of: maxQueryGap, maxTargetGap.
                        Default: 1e9,1e9.
                        That is: merge anything along the same sequences (in most cases where seqs are <<< 1e9 bp).
                        
                        ''')

parser.add_argument('--hardmerge', '-H',
                   action='store_true', default=False,
                   help='''To be used with --merge.
                    This will discard lines where the new/old query is encompassed by the old/new query regardless of target sequence.
                    In other words, it will give preference to a longer query alignment that encompasses a shorter one.
                    This can help strongly collapse a PAF file into alignments you care about.

                    You may actually want to run a regular merge first, followed by a hard merge.
                    Otherwise, this may not give a chance for downstream merges to be longer than an upstream merge.
                        ''')

parser.add_argument('--presorted', '-P',
                   action='store_true', default=False,
                   help='''Assume PAF is pre-sorted by qname, qstart, qend, tname, tstart, tend.
                        This can be done at commandline via "sort -k1,1 -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n".
                        This can save substantial time.
                        ''')

parser.add_argument('--sort', '-s',
                   action='store_true', default=False,
                   help='''Sort by qname, qstart, qend, tname, tstart, tend.
                        This can also be done at commandline via "sort -k1,1 -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n".
                        If --merge AND --sort are invoked, --sort is actually "silent" b/c --merge will sort (unless --presorted used) far above.
                        Then the merging output will come BEFORE the sort-only output in a series of if/elif statements.
                        In other words, one never needs to use --sort with --merge, but may want to use --presorted with --merge.
                        ''')

parser.add_argument('--verbose', '-v',
                   type=int, default=False,
                   help='''Provide int 1 or 2. Default: False.''')

args = parser.parse_args()


###
same_rules = dict(zip(['minLenRatio', 'maxLenRatio', 'minPctQueryMatch', 'minPctTargetMatch', 'minPctAlnMatch'], [float(e) for e in args.same_rules.strip().split(',')]))
bubble_rules = dict(zip(['minPctQueryInTarget','minPctAlnMatch','maxPctOfTarget'], [float(e) for e in args.bubble_rules.strip().split(',')]))
merge_rules = dict(zip(['maxQueryGap','maxTargetGap'], [float(e) for e in args.merge_rules.strip().split(',')]))

###

class Paf(object):
    def __init__(self, paffile):
        f2i = lambda x: int(round(float(x)))
        #self.fmt = {0:str, 1:int, 2:int, 3:int, 4:str, 5:str, 6:int, 7:int, 8:int, 9:int, 10:int, 11:int}
        self.fmt = {0:str, 1:f2i, 2:f2i, 3:f2i, 4:str, 5:str, 6:f2i, 7:f2i, 8:f2i, 9:f2i, 10:f2i, 11:f2i}
        self.file = paffile
        self.key = {'query':0, 'qlen':1,'qstart':2,'qend':3,'strand':4,'target':5,'tlen':6,'tstart':7,'tend':8,'match':9,'alnlen':10,'mapq':11}
        with open(self.file) as fh:
            # only grabs columns 1-12
            self.paf = [[self.fmt[i](linelist[i]) for i in range(len(linelist))] for linelist in [line.strip().split()[:12] for line in fh.readlines()]]
        
        self.iterpaf = iter(self.paf)

    def __iter__(self):
        return self

    def next(self):
        return self.iterpaf.next()

    def line2txt(self,line):
        return '\t'.join([str(e) for e in line])

    def get_dataframe(self):
        self.pafdf = pandas.DataFrame(self.paf, columns=['query', 'qlen','qstart','qend','strand','target','tlen','tstart','tend','match','alnlen','mapq']).sort_values(by=['target','tstart','tend'])
        return self.pafdf

    def reset_iter(self):
        ''' To use after iter complete, or after set iter to something else'''
        self.iterpaf = iter(self.paf)

    def set_iter(self, l):
        ''' iter over customized version of paf'''
        self.iterpaf = iter(l)
        
    def get_paf(self):
        return self.paf
    
    def sort_paf(self, key=None):
        '''Example key:
            key=lambda x: x[5] to sort on target names
            or
            key=lambda x: (x[5],x[7],x[8]) to sort on target names and start/ends'''
        if key is None:
            self.paf.sort()
        else:
            self.paf.sort(key=key)

    def sort_paf_by_query(self):
        self.sort_paf(key=lambda x: (x[0],x[2],x[3]))

    def sort_paf_by_target(self):
        self.sort_paf(key=lambda x: (x[5],x[7],x[8]))

    def sort_paf_by_query_then_target(self):
        self.sort_paf(key=lambda x: (x[0],x[2],x[3],x[5],x[7],x[8]))

    def merge_adj_ident_queries(self, maxqgap=1e9, maxtgap=1e9, presorted=False, hardmerge=False, verbose=True):
        if not presorted:
            self.sort_paf_by_query_then_target()
        iterable = iter(self.paf)
        custpaf = []
        curr = iterable.next()
        try:
            nqueries = 1
            i = 0
            while iterable:
                i+=1
                if verbose >= 2:
                    sys.stderr.write("Record "+str(i)+"\n")
                old = curr
                curr = iterable.next()
                qgap = curr[2] - old[3]
                tgap = curr[7] - old[8]
                #tests
                samequery = old[0] == curr[0]
                sametarget = old[5] == curr[5]
                sane_qgap = qgap <= maxqgap
                sane_tgap = tgap <= maxtgap
                old_q_enveloped = (old[2] >= curr[2]) and (old[3] <= curr[3])
                curr_q_enveloped = (curr[2] >= old[2]) and (curr[3] <= old[3])
                #old_t_enveloped = (old[7] >= curr[7]) and (old[8] <= curr[8])
                #curr_t_enveloped = (curr[7] >= old[7]) and (curr[8] <= old[8])
                #
                if samequery and sametarget and old_q_enveloped:
                    ## new q completely encompasses old record, then keep curr and discard old
                    ## this can mean old t was encompassed or not, but we will take the new t coords either way
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Old query enveloped... \n")
                    curr = curr
                    nqueries += 1
                elif samequery and sametarget and curr_q_enveloped: ## old q completely encompasses new record, then keep old and discard curr
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Curr query enveloped... \n")
                        sys.stderr.write('  '.join([str(e) for e in old])+"\n")
                        sys.stderr.write('  '.join([str(e) for e in curr])+"\n")
                        sys.stderr.write(' '.join([str(e) for e in [curr_q_enveloped, curr[2], old[2], curr[3], old[3]]])+"\n")
                    curr = old
                    nqueries += 1
                elif samequery and hardmerge and old_q_enveloped:
                    ## new q completely encompasses old record, then keep curr and discard old
                    ## this can mean old t was encompassed or not, but we will take the new t coords either way
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Old query enveloped... \n")
                    curr = curr
                    nqueries += 1
                elif samequery and hardmerge and curr_q_enveloped: ## old q completely encompasses new record, then keep old and discard curr
                    if verbose:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write("Curr query enveloped... \n")
                        sys.stderr.write('  '.join([str(e) for e in old])+"\n")
                        sys.stderr.write('  '.join([str(e) for e in curr])+"\n")
                        sys.stderr.write(' '.join([str(e) for e in [curr_q_enveloped, curr[2], old[2], curr[3], old[3]]])+"\n")
                    curr = old
                    nqueries += 1
                ## THE NEXT ARE COMMENTED OUT B/C THIS WAS INTENDED TO BE QUERY-BASED WHEREAS THE NEXT TWO LINES CAN RESULT IN NON-MERGING OF OBVIOUS QUERY MERGES
##                elif samequery and sametarget and old_t_enveloped:
##                    if verbose:
##                        sys.stderr.write("Old target enveloped... \n")
##                    curr = curr
##                    nqueries += 1
##                elif samequery and sametarget and curr_t_enveloped:
##                    if verbose:
##                        sys.stderr.write("Curr target enveloped... \n")
##                        sys.stderr.write('  '.join([str(e) for e in old])+"\n")
##                        sys.stderr.write('  '.join([str(e) for e in curr])+"\n")
##                        
##                    curr = old
##                    nqueries += 1
                elif samequery and sametarget and sane_qgap and sane_tgap: ##COMES AFTER ENVELOPED IFs ON PURPOSE
                    nqueries += 1
                    # then merge and set merge -> curr (curr will be set to "old" in next iter)
                    q = old[0]
                    qlen = old[1]
                    qstart = min(old[2], curr[2]) #should be old[2]
                    qend = max(old[3], curr[3]) #should be curr[3]
                    strand = old[4]
                    t = old[5]
                    tlen = old[6]
                    tstart = min(old[7], curr[7]) #should be old[7]
                    tend = max(old[8], curr[8]) #should be curr[8]
                    summatch = old[9] + curr[9] 
                    sumaln = old[10] + curr[10]
                    matchrate = float(summatch)/sumaln
                    errorrate = 1-matchrate
                    gaplen = max(qgap,tgap)
                    summatch_w_gap = round(summatch + min(0, matchrate*gaplen)) #adds 0 for positive gaplens, subtracts for neg gap lens (alns overlapped and matches inflated)
                    sumaln_w_gap = round(sumaln + gaplen) ## Adds gaplen for pos, subtracts for neg
                    mapq = round((old[11]+curr[11])/2.0) ## NOTE this is not weighted for longer alignments, etc -- uniform weighting
                    merged = [q, qlen, qstart, qend, strand, t, tlen, tstart, tend, summatch_w_gap, sumaln_w_gap, mapq]
                    #
                    if verbose >= 2:
                        sys.stderr.write("Record "+str(i)+"\n")
                        sys.stderr.write(self.line2txt(old)+"\n")
                        sys.stderr.write(self.line2txt(curr)+"\n")
                        sys.stderr.write(" ".join([str(e) for e in [qgap, tgap, summatch, sumaln, matchrate, gaplen, summatch_w_gap, sumaln_w_gap,mapq]])+"\n")
                    curr = merged
                    if verbose:
                        sys.stderr.write(self.line2txt(curr)+"\n")
                    
                else:
                    #return old
                    custpaf.append(old+[nqueries])
                    nqueries = 1
                        
        except StopIteration:
            pass
        #Don't forget last one
        custpaf.append(old+[nqueries])
        return custpaf
        
def verbose(msg):
    if not msg.endswith('\n'):
        msg += '\n'
    if args.verbose:
        sys.stderr.write(msg)


verbose("FilterPaf:: Loading PAF....\n")
paf = Paf(args.paf)


## Both bubbles and same operations can be preceded by the merge operation
## If neither are selected, the merge is returned later
verbose("FilterPaf:: Processing....\n")


## PRE-MERGING (AND SORTING IF NOT PRESORTED) 
if args.merge:
    rules = merge_rules
    paf.set_iter( paf.merge_adj_ident_queries(maxqgap=merge_rules['maxQueryGap'], maxtgap=merge_rules['maxTargetGap'], presorted=args.presorted, hardmerge=args.hardmerge, verbose=args.verbose) )



## IF/ELIF/ELSE STATEMENTS
if args.same:
    rules = same_rules
    #print rules
    #save = []
    for line in paf:
        qlen = float(line[paf.key['qlen']])
        tlen = float(line[paf.key['tlen']])
        match = float(line[paf.key['match']])
        alnlen = float(line[paf.key['alnlen']])
        pass_len_ratio = (tlen/qlen >= rules['minLenRatio'] and tlen/qlen <= rules['maxLenRatio'])
        pass_query_aln = (match/qlen >= rules['minPctQueryMatch'])
        pass_target_aln = (match/tlen >= rules['minPctTargetMatch'])
        pass_aln = (match/alnlen >= rules['minPctAlnMatch'])
        #print qlen, tlen, match, alnlen, pass_len_ratio, tlen/qlen, pass_query_aln, match/qlen, pass_target_aln, match/tlen, pass_aln, match/alnlen
        if pass_len_ratio and (pass_query_aln and pass_target_aln and pass_aln):
            print paf.line2txt(line)

        
        
elif args.bubble:
    rules = bubble_rules
    #print rules
    if args.reverse_bubble:
        pass
    for line in paf:
        qlen = float(line[paf.key['qlen']])
        tlen = float(line[paf.key['tlen']])
        qalnlen = float(line[paf.key['qend']]) - float(line[paf.key['qstart']])
        match = float(line[paf.key['match']])
        alnlen = float(line[paf.key['alnlen']])
        if args.same_seq_allowed:
            not_same_seq = True
        else:
            not_same_seq = line[paf.key['query']] != line[paf.key['target']]
        min_query = (qalnlen/qlen >= rules['minPctQueryInTarget'])
        pass_aln = (match/alnlen >= rules['minPctAlnMatch'])
        max_target = (qlen/tlen <= rules['maxPctOfTarget'])
        #print qlen, tlen, match, alnlen, pass_len_ratio, tlen/qlen, pass_query_aln, match/qlen, pass_target_aln, match/tlen, pass_aln, match/alnlen
        if not_same_seq and min_query and pass_aln and max_target:
            print paf.line2txt(line)

elif args.merge: #Return merged iterable (merged iterable set above prior to --same and --bubble)
    for line in paf:
        print paf.line2txt(line)

elif args.sort:
    ## if --merge AND --sort are invoked, --sort is actually "silent" b/c --merge will sort (unless --presorted used) far above
    ## Then it will come BEFORE args.sort in the if/elif statements here.
    if not args.presorted:
        paf.sort_paf_by_query_then_target()
    for line in paf:
        print paf.line2txt(line)
else:
    for line in paf:
        print paf.line2txt(line)
