#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    Given table,
        keep only 1 line for each element in name column 
        that has the highest score in score column.

    Use case:
    BLASTP output on gene set, keep entry with highest bit score.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('table', metavar='table', nargs='*',
                   type= str, 
                   help='''Paths to as many tables as you want. They will be treated as one giant table.
                        This can also be left empty for standard in or specified as - or stdin.''')

parser.add_argument('-n', '--name_column', type=int, default=1,
                    help='''Column where elements/names are found.''')

parser.add_argument('-s', '--score_column', type=int, default=2,
                    help='''Column where scores are found.''')


parser.add_argument('-d', '--delim', type=str, default='\t',
                    help='''Delimiter. Default = tab.''')

args = parser.parse_args()




if len(args.table) == 0 or args.table[0] in ('stdin', '-') or args.table[0].startswith('<('):
    args.table = [sys.stdin]
    def opentable(x):
        return x
    def closetable(x):
        return x
else:
    def opentable(x):
        return open(x)
    def closetable(x):
        return x.close()

lines = {}
scores = {}
score_col = args.score_column - 1
name_col = args.name_column - 1
order = []

for table in args.table:
    tablefile = opentable(table)
    for line in tablefile:
        line = line.strip().split(args.delim)
        score = float(line[score_col])
        order.append(line[name_col])
        try:
            if score > scores[line[name_col]]:
                scores[line[name_col]] = score
                lines[line[name_col]] = line
        except:
            scores[line[name_col]] = score
            lines[line[name_col]] = line
    closetable(tablefile)




## Output in order of appearance in input
for name in order:
        try:
            print (args.delim).join(lines[name])
            lines.pop(name)
        except:
            ## The intention is to cause it to fail silently after seeing element once and printing
            ## However, this can cause it to fail silently in unanticipated scenarios as well.
            pass


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


