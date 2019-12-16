#!/usr/bin/env python2.7
import sys, argparse, re
from collections import defaultdict

parser = argparse.ArgumentParser(description="""
    This is experimental - and may not work appropriately with ALL BED files as it was only tested on ones I needed.
    
    Input: Sorted BED file with 4 columns: chr, start, end, name.
        To unix sort: sort -k1,1 -k2,2n -k3,3n
        Or: sortBed -i -

        OR better:
        Give a pre-computed slice file:
        awk 'OFS="\\t" {print $1,$2,"s",NR,$4,$5"\\n"$1,$3,"e",NR}' file.bed | sort -k1,1 -k2,2n > file.slice.txt

        Above AWK assumes 5 columns (chr, start, end, name, score).
        If no score col, do:
        awk 'OFS="\\t" {print $1,$2,"s",NR,$4,0"\\n"$1,$3,"e",NR}' file.bed | sort -k1,1 -k2,2n > file.slice.txt
    Output:
        Updated BED file s.t. overlapping intervals get comma-separated names for names that occurred there.
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--bedfile', '-f', 
                   type=str, default=False,
                   help='''BED file. Can be "stdin", "-", or "<()" as well. ''')

parser.add_argument('--slicefile', '-s',
                   type=str, default=False,
                   help='''Slice file as described above. ''')

parser.add_argument('--maxScoreOnly', '-M',
                   action='store_true', default=False,
                   help='''For each sliced and diced interval, only report the name of the element that had the highest score over those bases (where score in column 5)''')

parser.add_argument('--minScoreOnly', '-m',
                   action='store_true', default=False,
                   help='''For each sliced and diced interval, only report the name of the element that had the lowest score over those bases.''')

parser.add_argument('--sortedNames', '-S',
                   action='store_true', default=False,
                   help='''Does not work with -m/-M. By default, the collapsed names and scores in same order are given in 4th and 5th column.
                    This pushes those to the 5th and 6th columns, and provides a different names column for the 4th where
                    the names are sorted s.t. any time a set of names occurs in the column, the should specify the same permutation.
                    This makes it easier to computer sums of intersecting intervals with AWK, for example.''')

args = parser.parse_args()


patterns = defaultdict(int)





## FUNCTIONS
def get_connection(fh):
    stdin = fh in ['stdin', '-'] or fh[:1] == '<('
    if stdin:
        return sys.stdin
    else:
        return open(fh), stdin

def close_connection(f, stdin):
    if not stdin:
        f.close()


class Interval(object):
    def __init__(self, line):
        self.line = line.strip().split()
        self.linelen = len(self.line)
        self.chr = self.line[0]
        self.start = int(self.line[1])
        self.end = int(self.line[2])
        self.name = self.line[3]
        if self.linelen >=5:
            self.score = float(self.line[4])
        if self.linelen >= 6:
            self.strand = self.line[5]
            
    def __str__(self):
        return '\t'.join(self.line).strip()
        

class BED4(object):
    def __init__(self, fh, all_in_mem=True, learn_chr=False):
        self.all_in_mem = all_in_mem
        self.fn = fh
        self.resetiter()

    def __iter__(self):
        return self

    def next(self):
        try:
            return Interval(self.iterbed4.next())
        except Exception as e:
            raise StopIteration

    def process_line(line):
        return Interval(line)
    def close(self):
        if not self.closed:
            close_connection(self.bed4, self.stdin)
            self.closed = True
    def resetiter(self):
        self.bed4, self.stdin = get_connection(self.fn)
        self.closed = False
        self.iterbed4 = self.bed4
        if self.all_in_mem:
            self.iterbed4 = iter(self.bed4.readlines())
            self.close()

    
def slice_and_dice(bed, maxScoreOnly=False, minScoreOnly=False):
    #break_points = defaultdict(set)
    #curr_chr = None
    base_names = defaultdict(dict)
    base_scores = defaultdict(dict)
    for b in bed:
        #break_points[b.chr].add(b.start)
        #break_points[b.chr].add(b.end)
        for i in range(b.start,b.end):
            if i not in base_names[b.chr].keys():
                base_names[b.chr][i] = set([])
                if maxScoreOnly and b.linelen >= 5:
                    base_scores[b.chr][i] = float('-inf')
                elif minScoreOnly and b.linelen >= 5:
                    base_scores[b.chr][i] = float('inf')
                else:
                    base_scores[b.chr][i] = 0
                    
            if maxScoreOnly and b.linelen >= 5:
                if b.score > base_scores[b.chr][i]:
                    base_names[b.chr][i] = set([b.name])
                    base_scores[b.chr][i] = b.score
            elif minScoreOnly and b.linelen >= 5:
                if b.score < base_scores[b.chr][i]:
                    base_names[b.chr][i] = set([b.name])
                    base_scores[b.chr][i] = b.score
            else:
                base_names[b.chr][i].add(b.name)
                if b.linelen >= 5:
                    base_scores[b.chr][i] += b.score
                else:
                    base_scores[b.chr][i] = 0

    lines = []
    for chrm in base_names.keys():
        start = None
        name = None
        end = None
        score = None
        for base in base_names[chrm]:
            if start is None:
                start = base
            if name is None:
                name = base_names[chrm][base]
            if score is None:
                score = base_scores[chrm][base]
            if end is None:
                end = base+1 ## since the last base WAS labeled as such and BED interval doesn't include last
            if base_names[chrm][base] != name:
                lines.append(Interval('\t'.join([str(e) for e in [chrm, start, end, ','.join(sorted(list(name))), score]])))
                name = base_names[chrm][base]
                score = base_scores[chrm][base]
                start = base
                end = base+1
            else:
                end = base + 1
        lines.append(Interval('\t'.join([str(e) for e in [chrm, start, end, ','.join(sorted(list(name))), score] ])))
                            
    return lines
    
        
def slice2(bed):
    # can have it create a slice file on the fly in the future...
    # for now I prefer just giving the slice file
    chrom = None
    start = {}
    end = {}
    name = {}
    score = {}
    overlappingIntervalCount = 0
    for b in bed:
        if chrom is None:
            chrom = b.chr
        if b.chr != chrom:
            process_cluster(chrom, start, end, name, score)
        else: #same chr
            if overlappingIntervalCount > 0 and b.start > max(end.values()): # new upcoming cluster, process old cluster
                process_cluster(chrom, start, end, name, score)
            else:
                start[overlappingIntervalCount] = b.start
                end[overlappingIntervalCount] = b.end
                name[overlappingIntervalCount] = b.name
                score[overlappingIntervalCount] = b.score if b.linelen >= 5 else 0
    #Process last cluster
    process_cluster(chrom, start, end, name, score)
                


class SliceLine(object):
    def __init__(self,line):
        self.line = line.strip().split()
        self.linelen = len(line)
        self.chr = self.line[0]
        self.coord = int(self.line[1])
        self.instruct = self.line[2]
        self.idx = int(self.line[3])
        if self.line[2] == 's':
            self.name = self.line[4]
            self.score = self.line[5]
    def __str__(self):
        return '\t'.join(self.line)

def slice3(slice_fh, maxScoreOnly=False, minScoreOnly=False, sortedNames=False):
    chrom = None
    #start = {}
    #end = {}
    name = {}
    score = {}
    #idxs = []
    last_coord = None
    with open(slice_fh) as f:
        for line in f:
            line = SliceLine(line)
            if chrom is None:
                chrom = line.chr
            if line.chr != chrom:
                assert name == {}
                last_coord = None
                chrom = line.chr
            if last_coord is not None and line.coord > last_coord:
                if maxScoreOnly:
                    idx = max(score, key=score.get)
                    names = name[idx]
                    scores = score[idx]
                elif minScoreOnly:
                    idx = min(score, key=score.get)
                    names = name[idx]
                    scores = score[idx]
                else:
                    #names = ','.join([str(e) for e in name.values()])
                    #scores = ','.join([str(e) for e in score.values()])
                    names = ','.join([str(name[e]) for e in sorted(name.keys())])
                    scores = ','.join([str(score[e]) for e in sorted(score.keys())])
                if sortedNames:
                    snames = ','.join([str(e) for e in sorted(name.values())])
                    bed = [chrom, last_coord, line.coord, snames, names, scores]
                    ## names and scores should be in same order whereas snames makes it easier to compute sums outside of script
                else:
                    bed = [chrom, last_coord, line.coord, names, scores] ## names and scores should be in same order
                print '\t'.join([str(e) for e in bed])
            if line.instruct == 's':
                name[line.idx] = line.name
                score[line.idx] = line.score
            if line.instruct == 'e':
                name.pop(line.idx)
                score.pop(line.idx)
            if name == {}:
                last_coord = None
            else:
                last_coord = line.coord
            
                
            

def main(args):
    assert not (args.maxScoreOnly and args.minScoreOnly)
    assert not (args.bedfile and args.slicefile)
    if args.bedfile:
        # painfully slow for the moment
        # can have it create a slice file on the fly in the future...
        # for now I prefer just giving the slice file
        # That will be in the fxn slice2 above
        for line in slice_and_dice(BED4(args.bedfile), maxScoreOnly=args.maxScoreOnly, minScoreOnly=args.minScoreOnly):
            print line
    elif args.slicefile:
        slice3(args.slicefile, maxScoreOnly=args.maxScoreOnly, minScoreOnly=args.minScoreOnly, sortedNames=args.sortedNames)


## EXECUTE:






main(args)
