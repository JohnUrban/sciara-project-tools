#!/usr/bin/env python2.7
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
Options: max, min, mean, median, sum, list, distinct.''')

parser.add_argument('--addcol', '-c3', type=str, default=False,
                    help='''Legacy (see --othercols too): When using ONLY max or min, also report these columns -- provide comma-separated list.''')

parser.add_argument('--othercols', '-c4', type=str, default=False,
                    help='''Collapse these other columns (comma-sep list) into their own columns in output as well. To get column of line for min or max see --addcol.''')


parser.add_argument('--skip', '-s', type=str, default='#',
                    help='''Skip lines that start with given string. ''')

parser.add_argument('--header','-H', action='store_true', default=False,
                    help='''Print header line''')

parser.add_argument('--strings', '-str', action='store_true', default=False,
                    help=''' Default treats elements in columns as floats. Treat as strings instead --operation must be set to list....''')

args = parser.parse_args()

if args.strings:
    assert args.operation in ["list", "distinct"]

c1 = args.column-1
c2 = args.column2-1

fxns = args.operation.split(",")

addcol=False
othercols=False
if args.addcol and len(fxns) == 1 and fxns[0] in ('min','max','list'):
    addcol=True
    cols = [int(e)-1 for e in args.addcol.strip().split(",")]
    D = {}
    for col in cols:
        D[col] = defaultdict(list)
if args.othercols:
    othercols=True
    kols = [int(e)-1 for e in args.othercols.strip().split(",")]
    DD = {}
    for col in kols:
        DD[col] = defaultdict(list)
##    if fxns[0] == 'max':
##        addfxn = max
##    elif fxns[0] == 'min':
##        addfxn = min

ops = []
listfxn = lambda x: (",").join([str(e) for e in x])
distinctfxn = lambda x: (",").join([str(e) for e in list(set(x))])
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
        fxn = listfxn
    elif e == 'distinct':
        fxn = distinctfxn
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
    if args.strings:
        d[line[c1]].append(str(line[c2]))
    else:
        d[line[c1]].append(float(line[c2]))
    if addcol:
        for col in cols:
            D[col][line[c1]].append(line[col])
    if othercols:
        for col in kols:
            DD[col][line[c1]].append(line[col])

if args.header:
    print (args.delimiter).join(["#c1_element", "number_found"] + fxns)
for e in sorted(d.keys()):
    ans = []
    for op in ops:
        ans.append( str(op(d[e])) )
    length = len(d[e])
    if addcol:
##        print d[e]
##        print ans
        idx = d[e].index(float(ans[0]))
        addcols = []
        for col in cols:
            addcols.append( D[col][e][idx] )
        print (args.delimiter).join([ str(e), str(length) ] + ans + addcols )                   
    elif othercols:
        othercols = []
        for col in kols:
            othercols.append( listfxn(DD[col][e]) )
        print (args.delimiter).join([ str(e), str(length) ] + ans + othercols )                   

    else:    
        print (args.delimiter).join([ str(e), str(length) ] + ans)
