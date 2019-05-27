#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    G


    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('sets', nargs='*', 
                   type= str, 
                   help='''Path to set files (single column files).
                        ''')
outtype = parser.add_mutually_exclusive_group(required=True)

outtype.add_argument('-i', '--intersect', action='store_true', default=False,
                    help='''Return intersection.''')
outtype.add_argument('-u', '--union', action='store_true', default=False,
                    help='''Return union.''')
outtype.add_argument('-d', '--diff', type=int, default=False,
                    help='''Return those unique to given index. Need to provide 0-based integer.''')

parser.add_argument('-s', '--sort',
                   action='store_true', default=False,
                   help='''Sort output (as strings).
                        ''')

args = parser.parse_args()





     
class SetEntry(object):
    def __init__(self, entry):
        self.entry = entry
        with open(entry) as f:
            self.set = set([e.strip() for e in f.readlines()])
    def get(self):
        return self.set


class Sets(object):
    # Will hold dictionary of unique names
    # Each unique name can be associated with multiple entries
    # This will  have functions to extract summary info from each unique name
    def __init__(self):
        self.all = {}
        self.idx = -1
    def add(self, entry):
        self.idx += 1
        self.all[self.idx] = SetEntry(entry)
    def get_intersection(self):
        assert self.idx >= 0
        ans = self.all[0].get()
        if self.idx >= 1:
            for i in range(1,self.idx+1):
                ans = ans.intersection(self.all[i].get())
        return ans
    def get_union(self):
        assert self.idx >= 0
        ans = self.all[0].get()
        if self.idx >= 1:
            for i in range(1,self.idx+1):
                ans = ans.union(self.all[i].get())
        return ans
    def get_difference(self, idx):
        assert self.idx >= 0
        assert idx <= self.idx
        ans = self.all[idx].get()
        if self.idx >= 1:
            others = range(self.idx+1)
            others.pop(idx)
            for i in others:
                ans = ans.difference(self.all[i].get())
        return ans

        


def print_set(x, sort=False):
    if sort:
        print '\n'.join(sorted(list(x)))
    else:
        print '\n'.join(list(x))

# Read in sets
sets = Sets()
for setpath in args.sets:
    sets.add(setpath)

# Return ans
try:
    if args.intersect:
        print_set( sets.get_intersection(), args.sort)
    elif args.union:
        print_set( sets.get_union() , args.sort)
    elif args.diff is not False:
        print_set( sets.get_difference(args.diff) , args.sort)
        
    
except IOError:
    # Broken Pipe
    pass

