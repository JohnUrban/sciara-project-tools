#!/usr/bin/env python2.7
import sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="""
    By John Urban. Apr 17, 2017.

    Input: a tab-delimited file, sorted on the column you'd like to split on.

    Output: Multiple tab-delimited files, split on adjacently-identical elements in your target column.

    Example:
    Input = BED file of BLAST hits for transcript sequences denoted in column 4 -- sorted on those transcript names in col 4.
    Output = One BED file for each transcript named in column 4 of input. 
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument("-i", "--input",
                   type= str, default=False, required=True,
                   help='''Path to tab-delim file. Can be '-' or 'stdin' to specify standard input.''')

parser.add_argument("-k", "--column",
                   type=int, default="4",
                   help='''1-based column number containing sorted elements to split on.''')

parser.add_argument("-s", "--splitelements",
                   type=str, default=False,
                   help='''By default, it assumes the element in the target column is in the final form for all operations.
This option allows you to provide a delimiter character to split the string on, and a 0-based index to tell it which sub-element to take as the 'final form'.
Provide those 2 parameters separated by a comma.
Example: '-s :,0' means to split the target elements at ':' and take the 0th element.''')

nonadj = parser.add_mutually_exclusive_group()

nonadj.add_argument("-n", "--notadjacent", action='store_true', default=False,
                   help='''By default, it is assumed the file is sorted and all lines to be split off together into a new file are adjacent.
This assumption allows 1 output file to be opened at a time.
This flag indicates that elements to be split off may not be adjacent.
A different strategy is needed.
The lines of the input file can be stored in memory and marked for each output file until the end of the input is reached.
Then one output file can be opened at a time.
This can be problematic if the file is too large to go into memory.
A different strategy is needed again. See --multipleopenfiles. ''')

nonadj.add_argument("-m", "--multipleopenfiles", action='store_true', default=False,
                   help='''By default, it is assumed the file is sorted and all lines to be split off together into a new file are adjacent.
That allows 1 output file to be opened at a time.
Another strategy would be to keep all output files open until the end of the input file is reached.
This would be useful if there are identical elements in the target column that are not adjacent, but should be included in the same output file.
--multipleopenfiles takes an approach where new output files are opened
each and every time a new element in the target column is found.
Non-adjacent lines with identical target elements can then be written to the same file.
This can cause problems of too many files being opened at the same time.
If that is the case and if you have enough memroy to load the entire input into memory,
then try --notadjacent instead.''')


parser.add_argument("-p", "--prependinputname", action='store_true', default=False,
                   help='''By default, the output filenames only have the element name.
This prepends the input file prefix to those output filenames, which can be useful to know where a collection
of output files came from.''')

parser.add_argument("-P", "--prefix", type=str, default=False,
                   help='''By default, the output filenames only have the element name.
This prepends a given prefix to those output filenames, which can be useful to know where a collection
of output files came from.
If used with --prependinputname, this will go in front of it.''')

parser.add_argument("-S", "--suffix", type=str, default=False,
                   help='''By default, the output filenames only have the element name.
This appends a given suffix to those output filenames (before file extension).
''')

parser.add_argument("-E", "--file_extension", type=str, default='.bed',
                   help='''File extension of output files. Default: .bed -- TODO allow automatically use input ext (and default when stdin).''')



parser.add_argument('--dir', '-d',
                   type= str, default=False,
                   help='''Directory to write output files to. Default: pwd. ''',)


args = parser.parse_args()

## put col number in pywanese
col = args.column - 1

ext = args.file_extension

outdir = ""
if args.dir:
    outdir = args.dir
    if not outdir.endswith("/"):
        outdir += "/"

if args.input in ("-", "stdin"):
    f = sys.stdin
else:
    f = open(args.input)

pre = ''
if args.prependinputname:
    pre += ('.').join(args.input.strip().split('.')[:-1]) + "-"
if args.prefix:
    pre += args.prefix + "-"

sfx = ""
if args.suffix:
    sfx = "-" + args.suffix

splitelements = False
if args.splitelements:
    splitelements = args.splitelements.split(",")
    spliton = splitelements[0]
    splitidx = int(splitelements[1])

def adjacent_singlefile(f):
    old = None
    for line in f:
        current = line.strip().split()[col]
        if splitelements:
            current = current.split(spliton)[splitidx]
        if old == None:
            out = open(outdir + pre + current + sfx + ext, 'w')
        elif current != old:
            out.close()
            out = open(outdir + pre + current + sfx + ext, 'w')
            old = current
        out.write(line)
        old = current
    out.close()

def nonadjacent_singlefile(f):
    outfiles = defaultdict(list)
    for line in f:
        current = line.strip().split()[col]
        if splitelements:
            current = current.split(spliton)[splitidx]
        outfiles[current].append( line )
    for current in outfiles.keys():
        out = open(outdir + pre + current + sfx + ext, 'w')
        for line in outfiles[current]:
            out.write( line )
        out.close()

def nonadjacent_multifile(f):
    outfiles = {}
    for line in f:
        current = line.strip().split()[col]
        if splitelements:
            current = current.split(spliton)[splitidx]
        if current not in outfiles.keys():
            outfiles[current] = open(outdir + pre + current + sfx + ext, 'w')
        outfiles[current].write( line )
    for current in outfiles.keys():
        outfiles[current].close()


if args.notadjacent:
    nonadjacent_singlefile(f)
elif args.multipleopenfiles:
    nonadjacent_multifile(f)
else:
    adjacent_singlefile(f)
