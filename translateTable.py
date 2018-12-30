#!/usr/bin/env python2.7
import sys, pandas
import argparse
import numpy as np
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

 Input:
     - table1 with N columns where 1 or more need translation from A->B
         - specify columns to translate
     - table2 with M columns where 2 columns are A and B key:val pairs
         - specify which coluns
 Output:
     - translated table1
     
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--table', '-i',
                   type=str, required=True,
                   help='''Path to input table1 file to be translated.''')

parser.add_argument('--columns', '-c',
                   type=str, required=True,
                   help='''Comma-separated list of columns (1-based) in table1 that need to be translated.''')

parser.add_argument('--dict', '-d',
                   type=str, required=True,
                   help='''Path to input table2 file that contains key and value columns.''')

parser.add_argument('--key', '-k',
                   type=int, required=True,
                   help='''Key/A column (1-based) in table2 - i.e. the values that need to be translated to something else.''')

parser.add_argument('--value', '-v',
                   type=int, required=True,
                   help='''Value/B column (1-based) in table2 - i.e. the replacement values.''')

parser.add_argument('--input_delimiter', '-ID',
                   type=str, default="\t",
                   help='''Input delimiter in table1-- default = tab.''')

parser.add_argument('--output_delimiter', '-OD',
                   type=str, default="\t",
                   help='''Output delimiter in translated table1-- default = tab.''')


parser.add_argument('--force', '-F',
                   action='store_true', default=False,
                   help='''By default, this script will fail if a key is encountered that was not anticipated. --force tells it to return the input key and keep moving.
                    In other words, there may be untranslated values in the output table that were not anticipated by your input dictionary.
                    Thie can be good, but is not default to avoid a "silent failure" in pipelines that don't want this.''')


args = parser.parse_args()


# Pythonize columns
kcol = args.key - 1
vcol = args.value - 1
tcols = [int(e)-1 for e in args.columns.strip().split(',')]


# Create key:val dict
transdict = {}
with open(args.dict) as d:
    for line in d:
        line = line.strip().split()
        transdict[line[kcol]] = line[vcol]

# Create key:val function 
if args.force:
    # that returns input value upon failure
    def translate(x):
        try:
            return transdict[x]
        except KeyError:
            return x
else:
    # that fails upon failure
    def translate(x):
        return transdict[x]

    
# Translate table
with open(args.table) as table:
    for line in table:
        line = line.strip().split(args.input_delimiter)
        print (args.output_delimiter).join([translate(line[i]) if i in tcols else line[i] for i in range(len(line)) ])
    
