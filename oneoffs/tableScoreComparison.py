#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    Given 2 tables and specifying columns for name and score, report:
        number elements unique to each file
        number elements that had highest score for each file

    Use case:
    Two different BLASTP outputs on same gene set,
        which gave better hits?

    If something appears more than once, highest score is used.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('table', metavar='table', nargs=2,
                   type= str, 
                   help='''Paths to table 1 and table 2.''')

parser.add_argument('-n', '--name_column', type=int, default=1,
                    help='''Column where elements/names are found.''')

parser.add_argument('-s', '--score_column', type=int, default=2,
                    help='''Column where scores are found.''')


parser.add_argument('-d', '--delim', type=str, default='\t',
                    help='''Delimiter. Default = tab.''')

args = parser.parse_args()





lines = {}
scores = {}
score_col = args.score_column - 1
name_col = args.name_column - 1


t = -1
tables = {}
for i in range(len(args.table)):
    table = open(args.table[i])
    tables[i] = {}
    for line in table:
        line = line.strip().split(args.delim)
        score = float(line[score_col])
        name = line[name_col]
        try:
            if score > tables[i][name]:
                tables[i][name] = score
        except:
            tables[i][name] = score
    table.close()




## Number unique to each
##found = set([])
##for i in range(len(args.table)):
##    found = found.union(set(tables[i].keys))
##print "Found", len(found), "names."
t1 = set(tables[0].keys())
t2 = set(tables[1].keys())
union = t1.union(t2)
print "Found", len(union), "elements in Tables 1 and 2."
print "Table 1 had", len(t1.difference(t2)), "unique elements."
print "Table 2 had", len(t2.difference(t1)), "unique elements."
print "Tables 1 and 2 shared", len(t1.intersection(t2)), "elements."

counts = {0:0, 1:0, "same":0}
for name in list(union):
    if name in t1 and name in t2:
        if tables[0][name] > tables[1][name]:
            counts[0] += 1
        if tables[0][name] < tables[1][name]:
            counts[1] += 1
        else:
            counts["same"] += 1

print "Table 1 had", counts[0], "higher scoring shared elements."
print "Table 2 had", counts[1], "higher scoring shared elements."
print "Tables 1 and 2 had", counts["same"], "equally scoring shared elements."

## BELOW DOES NOT WORK WITH STDIN
##for table in args.table:
##    tablefile = opentable(table)
##    for line in tablefile:
##        line = line.strip().split(args.delim)
##        try:
##            print (args.delim).join(lines[line[name_col]])
##            lines.pop(line[name_col])
##        except:
##            ## The intention is to cause it to fail silently after seeing element once and printing
##            ## However, this can cause it to fail silently in unanticipated scenarios as well.
##            pass
##    closetable(tablefile)


