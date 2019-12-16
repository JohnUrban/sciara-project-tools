#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION - combined individual KineticsTools GFFs.

Although designed for one purpose, it can be used on any GFF.

Assumes the first has all the sequences in header.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('gff', metavar='gff', nargs='+',
                   type= str, 
                   help='''Path to as many GFFs as need be.''')

parser.add_argument('--header', '-H', action='store_true', default=False, help='''Return only combined header''')
parser.add_argument('--entries','-E', action='store_true', default=False, help='''Return only GFF entries, not header.''')
parser.add_argument('--files','-F', action='store_true', default=False, help='''Return only GFF file names.''')

args = parser.parse_args()


d = {}

gfflines = []
srclines = []
seqlines = []

##print args.gff

sys.stderr.write("Found " + str(len(args.gff)) + " GFF files.\n")
if args.files:
    for line in args.gff:
        print line
    quit()


if not args.entries:
    ## Learn and print header
    
    ## First pass through GFFs: learn the combined header
    with open(args.gff[0]) as f:
            for line in f:
                if line:
                    if not line.startswith('#'):
                        break
                    if line.startswith('##source'):
                        if line not in srclines:
                            srclines.append(line)
                    elif line.startswith('##gff'):
                        if line not in srclines:
                            gfflines.append(line)
                    else:
                        seqlines.append(line)


    if len(args.gff) > 1:
        for gff in args.gff[1:]:
            #print gff
            with open(gff) as f:
                for line in f:
                    if line:
                        if not (line.startswith('##source') or line.startswith('##gff')):
                            break
                        if line.startswith('##source'):
                            if line not in srclines:
                                srclines.append(line)

    ## Print header
    for line in gfflines+srclines+seqlines:
        print line.strip()

## grab GFF entry lines
if not args.header:
    for gff in args.gff:
        with open(gff) as f:
                for line in f:
                    if line:
                        if line.startswith('#'):
                            continue
                        print line.strip()    

                        

