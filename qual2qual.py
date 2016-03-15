#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
import cStringIO as StringIO

parser = argparse.ArgumentParser(description = """

Given a fasta file and qual file that has encodings (like fastq) - generate qual file qith integers instead
--and vice versa

Version 2. Dec 16, 2015. John Urban

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('--fasta', '-f',
                   type= str, required=True,
                   help='''Path to fasta file.''')
parser.add_argument('--qual', '-q',
                   type= str, required=True,
                   help='''Path to qual file.''')
parser.add_argument('--spec', '-s',
                   type= str, required=False, default="",
                   help='''Since encoded quals can have ">" and this can be the first symbol on a line,
one may need more information to distinguish names. If all fasta entries have a shared feature in their name
(e.g. Contig), provide that to increase the specificity and ensure this program works.
If they do not have a shared feature, it is advised that one is given.
If there are several things where each alone does not capture all entries, but together they do -- then provide all
as a comma-separated list (e.g. Contig,Scaffold,Chr,Singleton)''')



args = parser.parse_args()

###################################################
################ FUNCTIONS ########################
###################################################
def make(f, specificity=""):
    if specificity != "":
        specificity = specificity.split(",")
    def checkspec(line, specificity=""):
        if specificity == "":
            return True
        for e in specificity:
            if e in line:
                return True
        return False
    F=open(f)
    d = {}
    seq = ""
    name = ">"
    for line in F:
        if line[0] == ">" and len(seq) == 0:
            #begin building first seq
            seqname = line.strip()[1:] 
            seq = ""
        elif line[0] == ">" and len(seq) >= 0 and checkspec(line, specificity):
            # if ">" is part of name, begin building subsequent seq
            # (else it is quality score and current sequence should continue building)
            d[seqname] = seq #enter last sequence into dict
            seqname = line.strip()[1:] # start new seq
            seq = ""
        else:
            #continue building current seq
            seq += line.strip()
    #last seq
    d[seqname] = seq
    F.close()
    return d





###################################################
################ EXECUTE ##########################
###################################################
fa = make(args.fasta, specificity=args.spec)
qual = make(args.qual, specificity=args.spec)

fq = ''
for name in fa.keys():
    fq += "@"+name+"\n"+fa[name]+"\n+\n"+qual[name]+"\n"

fqfile = StringIO.StringIO(fq)
SeqIO.convert(fqfile, "fastq", sys.stdout, "qual")
