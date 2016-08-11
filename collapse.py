#!/usr/bin/env python
import sys
import argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Take in a table/txt file and collapse unique elements of 1 column
while doing operation on elements of a 2nd column.

No need for it to be sorted.
It will treat any occurence of an element in column 1 as part of same set.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--input', '-i',
                   type=str, required=True,
                   help='''Path to input text file''')

parser.add_argument('--delimiter', '-d',
                   type=str, default="\t",
                   help='''Delimiter. Default is tab.''')

parser.add_argument('--column','-c', type=int, default=1,
                    help=''' Column to collapse into unique instances..''')

parser.add_argument('--column2', '-c2', type=int, default=2,
                    help='''Column to perform collapsin operation on...''')

parser.add_argument('--operation', '-o', type=str, default='sum',
                    help=''' Operation to perform on column2. Default: sum.
Options: max, min, mean, median, sum, list.''')

parser.add_argument('--skip', '-s', type=str, default='#',
                    help='''Skip lines that start with given string. ''')

parser.add_argument('--header','-H', action='store_true', default=False,
                    help='''Print header line''')

args = parser.parse_args()


c1 = args.column-1
c2 = args.column2-1

fxns = args.operation.split(",")

ops = []
for e in fxns:
    if e == 'sum':
        fxn = sum
    elif e == 'max':
        fxn = max
    elif e == 'min':
        fxn = min
    elif e == 'mean':
        fxn = np.mean
    elif e == 'median':
        fxn = np.median
    elif e == 'list':
        fxn = lambda x: (",").join([str(e) for e in x])
    ops.append(fxn)

d = defaultdict(list)

if args.input == "-" or args.input == "stdin":
    f = sys.stdin
else:
    f = open(args.input)
    
for line in f:
    if args.skip:
        if line.startswith(args.skip):
            continue
    line = line.strip().split(args.delimiter)
    d[line[c1]].append(float(line[c2]))

if args.header:
    print (args.delimiter).join(["#c1_element", "number_found"] + fxns)
for e in sorted(d.keys()):
    ans = []
    for op in ops:
        ans.append( str(op(d[e])) )
    length = len(d[e])
    
    print (args.delimiter).join([ str(e), str(length) ] + ans)
