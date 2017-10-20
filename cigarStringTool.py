#!/usr/bin/env python2.7
import sys
import argparse
import re
from collections import defaultdict
import pysam


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Count cigar string items.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--sam', '-s',
                   type= str, default=False,
                   help='''Input file in SAM format.''')
parser_input.add_argument('--bam', '-b',
                   type= str, default=False,
                   help='''Input file in BAM format.''')
parser_input.add_argument('--cmdline', '-c',
                   type= str, default=False,
                   help='''Input CIGAR strings (comma-separated if multiple). ''')

parser.add_argument('-t','--table',default=False,action='store_true',
                    help='''Put out tab-sep 9-column output of:
(1)queryname, (2)n-softclip, (3)n-match, (4)n-del, (5)n-ins, (6)query-len, (7)alignment length, (8)lenth of aligned segment in query, (9)length of aligned segment in reference''')

args = parser.parse_args()



def cigarCalculator(CIGAR, compiled_regex):
    parsedCIGAR = re.finditer(compiled_regex, CIGAR)
    CIGARcounts = defaultdict(int)
    for e in parsedCIGAR:
        CIGARcounts[str(e.group(0))[-1]] += int(str(e.group(0))[:-1])
        
    return ("\t").join([str(k)+":"+str(CIGARcounts[k]) for k in CIGARcounts.keys()]) + "\tlength:"+str(sum(CIGARcounts.values())) + "\talnlength:"+str(sum([CIGARcounts["M"], CIGARcounts["D"], CIGARcounts["I"]]))

## note using the cigar tuples is faster than cigarstrng approach for BAM/SAM
def samCigarCalculator(CIGAR, translate):
    CIGARcounts = defaultdict(int)
    for e in CIGAR:
        CIGARcounts[translate[e[0]]] += int(e[1])
    addons = "\ttotalSum:"+str(sum(CIGARcounts.values())) + "\tsumMDI_alnlen:"+str(sum([CIGARcounts["M"], CIGARcounts["D"], CIGARcounts["I"]])) + "\tsumMI_alnreadlen:"+str(sum([CIGARcounts["M"], CIGARcounts["I"]])) + "\tsumMD_alnreflen:"+str(sum([CIGARcounts["M"], CIGARcounts["D"]]))
    return ("\t").join([str(k)+":"+str(CIGARcounts[k]) for k in CIGARcounts.keys()]) + addons


def tableCigarCalculator(CIGAR, translate):
    CIGARcounts = defaultdict(int)
    for e in CIGAR:
        CIGARcounts[translate[e[0]]] += int(e[1])
    addons = "\t"+str(sum(CIGARcounts.values())) + "\t"+str(sum([CIGARcounts["M"], CIGARcounts["D"], CIGARcounts["I"]])) + "\t"+str(sum([CIGARcounts["M"], CIGARcounts["I"]])) + "\t"+str(sum([CIGARcounts["M"], CIGARcounts["D"]]))
    return ("\t").join([str(CIGARcounts[k]) for k in ["S","M","D","I"]]) + addons

if args.sam or args.bam:
    regexp = "[0-9]{1,}[A-Z]"
    compregexp = re.compile(regexp)
    if args.sam:
        samfile = pysam.Samfile(args.sam, "r")
    elif args.bam:
        samfile = pysam.Samfile(args.bam, "rb")
    translate = {0:"M", 1:"I",2:"D",3:"N",4:"S",5:"H",6:"P",7:"=",8:"X"}
    if args.table:
        for e in samfile:
            print e.qname + "\t" + tableCigarCalculator(e.cigar,translate) 
    else:
        for e in samfile:
            addons = "\ta/t/q-lens:" + str(e.alen) +","+ str(e.tlen) +","+ str(e.qlen)
            tags = ("\t").join([str(k[0])+":"+str(k[1]) for k in e.tags])
            print e.qname + "\t" + samCigarCalculator(e.cigar,translate) + addons +"\t" + tags
        
    
elif args.cmdline:
    regexp = "[0-9]{1,}[A-Z]"
    compregexp = re.compile(regexp)
    CIGARS = args.cmdline.split(",")
    for CIGAR in CIGARS:
        print cigarCalculator(CIGAR, compregexp)
