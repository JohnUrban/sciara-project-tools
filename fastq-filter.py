#!/usr/bin/env python
import gzip, sys, argparse
import numpy as np
from Bio import SeqIO
from fqEncoding import detectEncoding, rosettaStone, baseQC, baseQCmeans, baseQCquantile
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in fastq or fastq.sz file, performs filtering.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument("--orphans",
                   type= str, default=False,
                   help='''Input file of unpaired reads: fastq or fastq.gz''')
parser_input.add_argument('--paired',
                   type= str, default=False,
                   help='''Input files of paired reads: fastq or fastq.gz (comma-separated). ''')
parser_input.add_argument('--interleaved',
                   type= str, default=False,
                   help='''Input file of interleaved paired reads: fastq or fastq.gz (comma-separated). ''')

parser.add_argument('-N', "--filterN",
                   action="store_true", default=False,
                   help='''Remove reads with any N content''')

parser.add_argument('-Q', "--filterQ",
                   type=int, default=False,
                   help='''Remove read whose average quality goes below given cutoff. Default: False (no filtering based on quality). Range 0-40.''')

parser.add_argument('-B', "--filterBoth",
                   action="store_true", default=False,
                   help='''Always remove both reads in a pair if one is filtered.''')
parser_sep_or_append_orphans = parser.add_mutually_exclusive_group()
parser_sep_or_append_orphans.add_argument("--orphan-name", dest="orphan_name",
                   type=str, default=False,
                   help='''In paired mode, the default is to write orphans from both fwd and rev reads to same file.
                            The default filename is basically a concatenation of the two paired read file names (or interleaved file name twice).
                            This allows you to give it a better/shorter prefix -- so just provide a better/shorter prefix for combined orphans.
                            .fq.gz will be appended to whatever is given.''')
parser_sep_or_append_orphans.add_argument('-sep', "--separate-orphans", dest="separate_orphans",
                   action="store_true", default=False,
                   help='''In paired mode, the default is to write orphans from both fwd and rev reads to same file.
                            This flag will have them written to separate fwd and rev orphan files.
                            Cannot use with append_orphans''')
parser_sep_or_append_orphans.add_argument('-app', "--append-orphans", dest="append_orphans",
                   type=str, default=False,
                   help='''In paired mode, the default is to write orphans from both fwd and rev reads to same file.
                            This flag will have them appended to a pre-existing orpans file that you specify.
                            Cannot use with separate_orphans''')
parser.add_argument('-c', "--count",
                   action="store_true", default=False,
                   help='''Instead of filtering, just return counts of reads that pass/fail filters. Dry run.''')
parser.add_argument('-qh', "--qhist",
                   type=str, default=False,
                   help='''Provides a histogram of mean quality scores in reads. Can give an idea of how cutoffs might work. Provide filename.jpg or filename.pdf''')
parser.add_argument('-nc', "--Ncount",
                   action="store_true", default=False,
                   help='''It will just count how many reads have Ns. This is the same as specifying -N -c. This overrides use of -Q flag (sets it to False).''')
args = parser.parse_args()


def open_file_connection(fh):
    ## Open Connection to fq.gz -- this will also work with unzipped fq files
    connection = gzip.open(fh)
    try:
        connection.next()
        connection.seek(0)
    except IOError:
        connection = open(fh)
    return connection

def open_file_connection_appendmode(fh):
    ## Open Connection to fq.gz -- this will also work with unzipped fq files
    connection = gzip.open(fh)
    gz = True
    try:
        connection.next()
        connection.seek(0)
    except IOError:
        connection = open(fh)
        gz = False
    connection.close()
    if gz:
        connection = gzip.open(fh, 'ab')
    else:
        connection = open(fh, 'a')
    return connection

def get_next_read(connection):
    read = dict()
    read["name"] = connection.next().strip()
    read["seq"] = connection.next().strip()
    connection.next()
    read["qual"] = connection.next().strip()
    return read

def write_read(read, connection):
    read = ("\n").join([read["name"], read["seq"], "+", read["qual"]]) + "\n"
    connection.write(read)

def has_N(seq):
    if "N" in seq.upper():
        return True
    return False

def qual_mean(qual, qualTranslator):
    scores = []
    for char in qual:
        scores.append(qualTranslator[char])
    return np.mean(scores)

def below_Q(qual, Qcutoff, qualTranslator):
    if qual_mean(qual, qualTranslator) < Qcutoff:
        return True
    return False


def discard(args, read, qualTranslator):
    N = False
    Q = False
    if args.filterN:
        N = has_N(read["seq"])
    if args.filterQ:
        Q = below_Q(read["qual"], args.filterQ, qualTranslator)
    return (N or Q)

class EncodingError(Exception):
    def __init__(self):
        self.message = "Encoding between 2 files was detected as inconsistent"
    def __str__(self):
        return self.message

def mean_qual_plot(plot_file, means):
    plt.hist(means)
    if plot_file.endswith(".pdf") or plot_file.endswith(".jpg"):
            plt.savefig(plot_file)
    else:
            print "Unrecognized extension for %s! Try .pdf or .jpg" % (plot_file)
            sys.exit()
    plt.close()

def mean_qual_scatter(plot_file, means1, means2):
    plt.scatter(means1, means2, s=2)
    if plot_file.endswith(".pdf") or plot_file.endswith(".jpg"):
            plt.savefig(plot_file)
    else:
            print "Unrecognized extension for %s! Try .pdf or .jpg" % (plot_file)
            sys.exit()
    plt.close()

## input names
if args.orphans:
    f = args.orphans
    fh = open_file_connection(args.orphans)
elif args.interleaved:
    f = args.interleaved
    fh = open_file_connection(args.interleaved)
    fh2 = fh
elif args.paired:
    f, f2 = args.paired.split(",")
    fh = open_file_connection(f)
    fh2 = open_file_connection(f2)


## detect encoding
if args.filterQ or args.qhist:
    encoding1 = detectEncoding(fh)
    qualTranslator = rosettaStone(encoding1)
    if args.paired:
        encoding2 = detectEncoding(fh2)
        if encoding2 != encoding1:
            raise EncodingError
else: qualTranslator = None


### plot quality score means, then quit -- in future, just do at same time as filtering if want...
if args.qhist:
    means1 = []
    means2 = []
    try:
        while fh:
            read1 = get_next_read(fh)
            means1.append(qual_mean(read1["qual"], qualTranslator))
            if args.paired or args.interleaved:
                read2 = get_next_read(fh2)
                means2.append(qual_mean(read2["qual"], qualTranslator))
    except StopIteration:
        pass
    fh.close()
    mean_qual_plot("fwd_"+args.qhist, means1)
    meansout = open("means_"+(".").join(args.qhist.split(".")[:-1])+".txt", 'w')
    meansout.write((" ").join([str(e) for e in means1])+"\n")
    if args.paired:
        fh2.close()
        mean_qual_plot("rev_"+args.qhist, means2)
        mean_qual_scatter("scatter_"+args.qhist, means1, means2)
        meansout.write((" ").join([str(e) for e in means2])+"\n")
    meansout.close()
    quit()

## count reads with an N -- same as -N -c
if args.Ncount:
    args.filterN = True
    args.filterQ = False
    args.count = True

## output names
if not args.count:
    out = (".").join([e for e in f.split(".") if e != "fq" and e != "gz" and e != "fastq"])
    fo = gzip.open(out+".filtered.fq.gz", 'wb')
    if args.paired or args.interleaved:
        if args.paired:
            out2 = (".").join([e for e in f2.split(".") if e != "fq" and e != "gz" and e != "fastq"])
            fo2 = gzip.open(out2 +".filtered.fq.gz", 'wb')
        elif args.interleaved:
            out2 = out
            fo2 = fo
        if args.separate_orphans:
            fo_orphan1 = gzip.open(out+".filtered.orphans.fq.gz", 'wb')
            fo_orphan2 = gzip.open(out2+".filtered.orphans.fq.gz", 'wb')
        elif args.append_orphans:
            fo_orphan1 = open_file_connection_appendmode(args.append_orphans)
            fo_orphan2 = fo_orphan1
        else: #default combine
            if not args.orphan_name:
                args.orphan_name = out+"_"+out2+".filtered.combined_orphans"
            fo_orphan1 = gzip.open(args.orphan_name+".fq.gz", 'wb')
            fo_orphan2 = fo_orphan1


## main tasks
try:
    count = {'numread':0, 'rm1':0, 'rm2':0,'filterboth':0, 'rmboth':0}
    while fh:
        read1 = get_next_read(fh)
        count['numread'] += 1
        remove1 = discard(args, read1, qualTranslator)
        if args.paired or args.interleaved:
            read2 = get_next_read(fh2) ## nead to read it out no matter what for bookkeeping --e.g.staying in same place as read1 (paired) and/or moving up to the next read1 (interleaved)
            if args.filterBoth and remove1:
                count['filterboth'] += 1
                continue
            else:
                remove2 = discard(args, read2, qualTranslator)
            if args.filterBoth and remove2:
                count['filterboth'] += 1
                continue
        if remove1 and args.orphans:
            count['rm1'] += 1
            continue
        elif not remove1 and args.orphans: #keep read1, no pairs, and not dry run, then write
            if not args.count:
                write_read(read1, fo)
        elif remove1 and remove2:
            count['rmboth'] += 1
            continue
        # either read1 was F, r2 ws F, or both were F
        elif remove1 and (args.paired or args.interleaved):
            count['rm1'] += 1
            if not args.count:
                write_read(read2, fo_orphan2)
        elif remove2:
            count['rm2'] += 1
            if not args.count:
                write_read(read1, fo_orphan1)
        #both were Flse
        elif not remove1 and not remove2:
            if not args.count:
                write_read(read1, fo)
                write_read(read2, fo2)
        elif args.count:
            pass
        else:
            print "unanticipated outcome - debugging time"
                
        
except StopIteration:
    pass

## REPORT TO STDERR
if args.paired or args.interleaved:
    readtype = "paired end"
else:
    readtype = "single end"
msg = "Encountered %d %s reads." % (count['numread'], readtype)
sys.stderr.write(msg + "\n")
remain = count['numread']-count['rmboth']-count['filterboth']-count['rm1']-count['rm2']
msg = "%d %s reads remain (%f percent)" % (remain, readtype, 100.0*remain/count['numread'])
sys.stderr.write(msg + "\n")
if args.orphans:
    msg = "%d single end reads (%f percent) were discarded" %  (count['rm1'], 100.0*count['rm1']/count['numread'])
    sys.stderr.write(msg + "\n")
elif args.paired or args.interleaved:
    msg = "%d pairs (%f percent)were removed because at least one in the pair offended (filter both set)" % (count['filterboth'], 100.0*count['filterboth']/count['numread'])
    sys.stderr.write(msg + "\n")
    msg = "%d pairs (%f percent) were removed because both offended" %  (count['rmboth'], 100.0*count['rmboth']/count['numread'])
    sys.stderr.write(msg + "\n")
    msg = "%d fwd reads (%f percent) were discarded (making %d rev read orphans)" %  (count['rm1'], 100.0*count['rm1']/count['numread'], count['rm1'])
    sys.stderr.write(msg + "\n")
    msg = "%d rev reads (%f percent) were discarded (making %d fwd read orphans)" %  (count['rm2'], 100.0*count['rm2']/count['numread'], count['rm2'])
    sys.stderr.write(msg + "\n")
                                                                                                      


#closing up
if not args.count:
    fh.close()
    fo.close()
    if args.paired:
        fh2.close()
        fo2.close()
        fo_orphan1.close()
        fo_orphan2.close()
