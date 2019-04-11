#!/usr/bin/env python2.7
import sys, argparse, re
from collections import defaultdict

parser = argparse.ArgumentParser(description="""
    This is just a utility like grep or awk.
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--file', '-f',
                   type=str, required=True,
                   help='''File to search. Can be "stdin", "-", or "<()" as well. ''')

parser.add_argument('--patterns', '-p',
                   type=str, required=True,
                   help='''File of patterns.''')

parser.add_argument('--column', '-c',
                   type=int, default=0,
                   help='''Column to search. 1-based. Default is to match the entire line.
Given a column number, it will require a match to the column.
Using --scan will look for it in entire line or column rather than requiring the line or column to be exact match.''')

parser.add_argument('--patterns_column', '-C',
                   type=int, default=0,
                   help='''1-based Column to extract patterns from in the file of patterns. Default is to look at the entire line as a string.''')

parser.add_argument('--v', '-v',
                   action='store_true', default=False,
                   help='''Same as grep -v: return lines NOT matched.''')

parser.add_argument('--scan', '-S',
                   action='store_true', default=False,
                   help='''Scan method: look for exact match in entire line or column rather than requiring the line or column to be exact match..''')
parser.add_argument('--delim', '-s',
                   type=str, default='\t',
                   help='''Delimiter of input file. Default: tab.''')
args = parser.parse_args()


patterns = defaultdict(int)



pcolumn = args.patterns_column - 1
entire_line = pcolumn < 0
with open(args.patterns) as f:
    for line in f:
        if entire_line:
            patterns[line.strip()] = 1
        else:
            line = line.strip().split()
            patterns[line[pcolumn]] = 1

patterns_regex = re.compile('|'.join(patterns.keys()))
# determine pattern in function
if args.scan:
    def pattern_in(x):
        if len(re.findall(patterns_regex, x)) > 0:
            return True
        return False
else:
    def pattern_in(x): ## returns 1 if part of dict construction above, 0 otherwise
        return patterns[searchString]

column = args.column - 1
stdin = args.file in ['stdin', '-'] or args.file[:1] == '<('
if stdin:
    f = sys.stdin
else:
    f = open(args.file)
for line in f:
    #line = line.strip().split(args.delim)
    line = line.strip()
    # Get string to search
    searchString = line.strip().split(args.delim)[column] if column >= 0 else line
    if not args.v and pattern_in(searchString):
        #print (args.delim).join( line )
        print line
    elif args.v and not pattern_in(searchString):
        #print (args.delim).join( line )
        print line
        
if not stdin:
    f.close()
