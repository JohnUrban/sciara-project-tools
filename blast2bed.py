#!/usr/bin/env python2.7
import sys, argparse

parser = argparse.ArgumentParser(description="""

    Updates on Jan 4, 2019 may break your pipeline.
    It shouldn't. The output should be the same, and newer outputs are obtained by --force.
    Nonetheless, use --old to play it safe if you want.

    By John Urban.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument("-b", "--blastfile",
                   type= str, default=False, required=True,
                   help='''Tab-delim blast output file.''')
parser.add_argument("-k", "--columns",
                   type=str, default="1,2,3",
                   help='''Provide 3, 4, or 6 1-based, comma-separated column numbers for (in this order):
sequence name,
start,
end (for the subject),
strand,
name,
score.
Note: confusingly, this means if you set your blast file to look like BED6 (qseqid qstart qend stitle bitscore sstrand), you need to provide: 1,2,3,6,4,5.
Old_Default: 2,13,14,8.
New_Default: 1,2,3
The first 3 are required.
It is allowed to do the first 4 or all 6 as well.
If strand, name, or score is missing from BLAST file, then put -1 as the value.''')

parser.add_argument("-K", "--othercolumns",
                   type=str, default=False, 
                   help='''1-based Column number for other columns you want to trail.
Note: For strand information to be included in output, either use --force OR put that column in this list as well.''')

parser.add_argument("-F", "--force",
                   action='store_true', default=False, 
                   help='''By default, this script is mostly concerned with the first 3 BED columns (and strand).
This option tells it to try to make the first 6 BED-consistent columns: chr, start, end, name, score, strand, ..., .
In particular, this flag is useful if not using -B/--otherBEDcolumns.
''')

parser.add_argument("-D", "--delim",
                   type=str, default="\t", 
                   help='''Delimiter in BLAST file. Default is tab.
''')

parser.add_argument("-O", "--old",
                   action='store_true', default=False, 
                   help='''Restore older functionality to keep pipeline working.
''')
args = parser.parse_args()


###################################################################################
''' Process special delimiter args '''
###################################################################################
if args.delim == 'space':
    args.delim = ' '
elif args.delim == 'tab':
    args.delim = '\t'
elif args.delim == 'newline':
    args.delim = '\n'



###################################################################################
''' Process column information (make pythonish).'''
###################################################################################
cols = [int(e) - 1 for e in args.columns.split(",")]
name = cols[0]
start = cols[1]
end = cols[2]
ncols = len(cols)
strand = False
bedname = False
bedscore = False
if ncols > 3:
    try:
        strand = cols [3]
    except:
        strand = False
    try:
        bedname = cols[4]
    except:
        bedname = False
    try:
        bedscore = cols[5]
    except:
        bedscore = False

if strand < 0:
    strand = False
if bedname < 0:
    bedname = False
if bedscore < 0:
    bedscore = False

others = []
if args.othercolumns:
    others = [int(e) - 1 for e in args.othercolumns.split(",")]
if args.blastfile == "-" or args.blastfile =="stdin":
    blastfile = sys.stdin
else:
    blastfile = open(args.blastfile, 'r')

###################################################################################
''' Definitions/Functions '''
###################################################################################
def pos_strand(entry):
    return (strand != False and entry[strand] == "plus") or (int(entry[start]) < int(entry[end]))

def neg_strand(entry):
    return (strand != False and entry[strand] == "minus") or (int(entry[end]) < int(entry[start]))

def get_coords(entry):
    if pos_strand(entry):
        return [entry[name], str(int(entry[start])-1), entry[end]], '+'
    elif neg_strand(entry):
        return [entry[name], str(int(entry[end])-1), entry[start]], '-'


def outfmt(entry, outname=".", outscore=".", force=False):
    ##Note: it is using env variables
    if bedname:
        outname = ('_').join(entry[bedname].split())
    if bedscore:
        outscore = entry[bedscore]
    coords, sstrand = get_coords(entry)
    otherinfo = [entry[e] for e in others]
    otherbed = [outname, outscore, sstrand] if bedname or bedscore or force else []
    return ("\t").join(coords + otherbed + otherinfo)


###################################################################################
'''EXECUTE'''
###################################################################################

if not args.old:
    for entry in blastfile:
        print outfmt( entry.strip().split(args.delim) )

else:
    ## OLDER WAY
    for entry in blastfile:
        entry = entry.strip().split() ## args.delim not used here to keep with old way
        if (strand != False and entry[strand] == "plus") or (int(entry[start]) < int(entry[end])):
            print ("\t").join([entry[name], str(int(entry[start])-1), entry[end]] + [entry[e] for e in others])
        elif (strand != False and entry[strand] == "minus") or (int(entry[end]) < int(entry[start])):
            print ("\t").join([entry[name], str(int(entry[end])-1), entry[start]] + [entry[e] for e in others])
