#!/usr/bin/env python2.7
import sys, argparse, re
from collections import defaultdict
from itertools import product
from Bio import SeqIO


parser = argparse.ArgumentParser(description="""
    Like faCount from kentUtilities... only:
    - slower
    - allows specification of kmer sizes

    Only counts N for 1-mers.
    For k>1, kmers with N are not reported.
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('fastx', metavar='fastx', nargs='+',
                   type= str, 
                   help='''Path to fastx.''')

parser.add_argument('-k', '--k', type=str, default='1,2',
                    help = '''Comma-sep kmer-sizes. Default: 1,2.''')

parser.add_argument('--revcomp', default=False, action='store_true',
                    help = '''Also count kmers from reverse complement of each sequence.''')

parser.add_argument('--globalonly', '--summarize', '-g', default=False, action='store_true',
                    help = '''Only return output summarized across all sequences (global).''')

##parser.add_argument('--asKinToolsBgModel', default=False, action='store_true',
##                    help = '''Returns global output as a text file structured as what I've been using for kinToolsBG models -- likely not needed by anyone else.''')
##

parser.add_argument('--outfile', default=None, type=str,
                    help = '''Defaults to stdout.''')


parser.add_argument('--fastx_type', default='fasta', type=str,
                    help = '''Defaults to fasta. Can be fastq as well.''')
args = parser.parse_args()


def kmercount_in_string(string, kmerdict, klist):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k (a list of ksizes)'''
    if type(klist) == str:
        klist = list(klist)
    string = string.upper()
    stringlen = len(string)
    kmerdict['len'] = stringlen
    for i in range(stringlen-min(klist)+1):
        for k in klist:
            if stringlen-i >= k:
                kmerdict[string[i:i+k]] += 1
    return kmerdict

def write_out_entry(fhout, msg):
    ''' '''
    fhout.write(msg + '\n')
    

def getkmers(klist):
    kmers = {} # dict of k and kmers
    for k in klist:
        kmers[k] = [''.join(e) for e in product('ACGT',repeat=k)]
        if k == 1:
            kmers[k].append('N')
    return kmers


def get_header(kmers):
       return 'name\tlength\t' + '\t'.join([ '\t'.join(kmers[k]) for k in sorted(kmers.keys())])

def get_entry_str(name, kmers, kmerdict):
    return '\t'.join([name,str(kmerdict['len'])]+[str(kmerdict[kmer]) for k in sorted(kmers.keys()) for kmer in kmers[k] ])


def kmercount_in_fastx(fhlist, fastx='fasta', klist=[1,2], rev_comp=False, outfile=None, globalonly=False):
    ## ESTABLISH kmers that will be reported
    kmers = getkmers(klist)

    ## OPEN out connection
    fhout = open(args.outfile, 'w') if args.outfile is not None else sys.stdout

    ## HEADER
    write_out_entry(fhout, msg=get_header(kmers))

    ## Initialize global kmerdict for tally of all entries at end
    global_kmerdict = defaultdict(int)
    ## PROCESS
    for fh in fhlist:
        if fh in ('-', 'stdin') or fh.startswith('<('):
            fh = sys.stdin
        for fa in SeqIO.parse(fh, fastx):
            # renew kmerdict
            kmerdict = defaultdict(int)
            # collect kmer info
            if fa is not None:
                kmerdict = kmercount_in_string(str(fa.seq), kmerdict, klist)
                if rev_comp:
                    kmerdict = kmercount_in_string(str(fa.reverse_complement().seq), kmerdict, klist)
            # report kmer info for this entry
            if not globalonly:
                write_out_entry(fhout, msg = get_entry_str(str(fa.name), kmers, kmerdict) )

            # Add to global kmerdict
            for kmer in kmerdict.keys():
                global_kmerdict[kmer] += kmerdict[kmer]

    ## Report global
    write_out_entry(fhout, msg = get_entry_str('global', kmers, global_kmerdict) )
    
    ## CLOSE out connection
    if args.outfile is not None:
        fhout.close()

    ## Return global
    return global_kmerdict

def run(args):
    ## 1. Process args
    klist = [int(e) for e in args.k.strip().split(',')]
    

    

    ## 2. process
    global_kmerdict = kmercount_in_fastx(fhlist=args.fastx,
                       fastx=args.fastx_type,
                       klist=klist,
                       rev_comp = args.revcomp,
                       outfile=args.outfile,
                       globalonly=args.globalonly)





try:
    run(args)
except IOError:
    pass
