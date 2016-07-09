#!/usr/bin/env python
from Bio import SeqIO
import sys
import random
import datetime
import argparse

parser = argparse.ArgumentParser(description="""

Given a fasta file, return the fasta sequences with falcon-compatible headers.
New fasta names:
>m000_000/i/0_length OriginalHeader

In FALCON, DAZZ_DB/fasta2DB is very particular about fasta headers.
It's aimed at PacBio data.
It is looking for:
>name/wellnumber/number_number
It seems the "name" part can be any string.
However, it should be the same for all fasta entries in a file (can be different if newer versions of Falcon used-- See below).

Also, it seems that the num_num part can be a string of any 2 integers int1_int2.
There does not even seem to be a need for int1 < int2.
However, I suggest 0_readlen - why not give it some information?

The 'wellnumber' should be unique for each read.
As far as I can tell, well number can be any integer from 0 to zillions.
(Any read from same well and SMRT cell for PacBio is from same molecule).
DBsplit without -a will exclude all reads of the same wellnumber except the longest one.

Finally, it appears that you can add any additional information after a space (separating the pacbio-header from other info).
e.g.:
>name/wellnumber/number_number AnyOtherInfo

Older versions of fasta2DB required "name" to be exacly the same for every sequence in the file where each file was from a different SMRT cell.
For example, in a file named:
m131124_075049_42156_c100592722550000001823106305221480_s1_p0.fasta
Every subread entry name in the fasta file starts with:
>m131124_075049_42156_c100592722550000001823106305221480_s1_p0/well_number/start_end

In newer versions of fasta2DB, a file may contain the data from multiple SMRT cells provided the reads for each SMRT cell are consecutive
in the file.'
For example, in a file:
all_subreads.fasta
Made from combining two SMRT cell fasta files:
m131124_075049_42156_c100592722550000001823106305221480_s1_p0.fasta
m131124_110952_42156_c100592722550000001823106305221481_s1_p0.fasta
Should have all the reads from the first file appear before all the reads in the second file.
e.g. (without showing sequences after headers):
>m131124_075049_42156_c100592722550000001823106305221480_s1_p0/1/0_1000
>m131124_075049_42156_c100592722550000001823106305221480_s1_p0/2/0_1000
>m131124_075049_42156_c100592722550000001823106305221480_s1_p0/3/0_1000
>m131124_110952_42156_c100592722550000001823106305221481_s1_p0/1/0_1000
>m131124_110952_42156_c100592722550000001823106305221481_s1_p0/2/0_1000
>m131124_110952_42156_c100592722550000001823106305221481_s1_p03/0_1000

The movie name (the constant for each SMRT cell) seems like it can be ANY string, and does not have to start with "m".
For example - the following can be the constant part of every fasta name:
>m000_000
Or even:
>foo_bar



What is in a PacBio name?
Example PacBio header:
>m131124_075049_42156_c100592722550000001823106305221480_s1_p0/767/0_11853
>1222222_222222 33333 4444444444444444444444444444444444 55 66 777 8888888
1 = m = movie
2 = yymmdd_hhmmss
3 = instrument serial number
4 = SMRT cell barcode
5 = set number
6 = part number
7 = ZMW hole number
8 = subread coordinates on entire unrolled read (start_stop)


John Urban (2016)
    """, formatter_class= argparse.RawTextHelpFormatter)


inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('--fastx', '-f',
                   type= str,
                   help='''Path to fasta file.''')
inputtype.add_argument('--stdin', action='store_true',
                   help='''Fastx is coming from stdin stream.''')


filetype = parser.add_mutually_exclusive_group()
filetype.add_argument('--fa', action='store_true',
                   help='''Explicitly define as fasta input.''')
filetype.add_argument('--fq', action='store_true',
                   help='''Explicitly define as fastq input.''')



parser.add_argument('--moviename', '-n',
                   type=str, default='m000_000',
                   help='''This will be used as the movie name constant in all falconized fasta names.
As far as I can tell, it can be any string.
Default: m000_000''')

args = parser.parse_args()


################################################################################
''' Check args'''
################################################################################
# file or stdin?
if args.stdin:
    fastxFile = sys.stdin
else:
    fastxFile = open(args.fastx)

# FASTA or FASTQ?
if args.fa:
    fastx = "fasta"
elif args.fq:
    fastx = "fastq"
elif args.fastx: ## can figure out from file (right now not from stdin due to non-seekability)
    line1 = fastxFile.next()[0]
    if line1[0] == ">":
        fastx = "fasta"
    elif line1[0] == "@":
        fastx = "fastq"
    fastxFile.seek(0)
else:
    print "Expected fasta or fastq. File given needs to be reformatted if user thinks it is."
    quit()



################################################################################
'''Run'''
################################################################################


i=0
for fa in SeqIO.parse(args.fastx, fastx):
    i+=1
    print ">"+ ("/").join([args.moviename, str(i), str(0)+"_"+str(len(fa))]) + " " + fa.description
    print str(fa.seq)

