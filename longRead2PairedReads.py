#!/usr/bin/env python2.7
import sys, gzip, argparse, random
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Given a fasta/fastq file of long reads (or any sequences), returns a pair of fastx files of paired-reads generated from the long reads.

TODO:
-- give each subread of a given long read, unique names. Now they all have the long read name.
-- randomize fragment start positions: identify # of fragments that would be sampled regularly, and sample that number randomly instead of in fixed steps
-- randomize insert size with given mean abd std dev: int(random.normalvariate(mu, sigma))
-- or just allow comma-separted multiple frag sizes to randomly select along each read....

    """, formatter_class= argparse.RawTextHelpFormatter)


inputtype = parser.add_mutually_exclusive_group(required=True)
inputtype.add_argument('--fastx', '-f',
                   type= str,
                   help='''Path to fasta/fastq file. Can be gzipped (file.fa.gz)''')
inputtype.add_argument('--stdin', action='store_true',
                   help='''Fastx is coming from stdin stream. Cannot be gzipped as stdin (for now).''')


filetype = parser.add_mutually_exclusive_group()
filetype.add_argument('--fa', action='store_true',
                   help='''Explicitly define as fasta input.''')
filetype.add_argument('--fq', action='store_true',
                   help='''Explicitly define as fastq input.''')


outtype = parser.add_mutually_exclusive_group()
outtype.add_argument('--ofa', action='store_true',
                   help='''Return output as fasta.''')
outtype.add_argument('--ofq', action='store_true',
                   help='''Return output as fastq.''')

parser.add_argument('-z', '--gzip', action="store_true", default=False,
                    help=''' Return paired reads in gzipped files.''')

parser.add_argument('-l','--fraglen',type=int,default=300,
                    help='''specify fragment length of synthetic paired-end reads.
This will be the length from beginning of read 1 to end of read2.
Can also think of it as the sliding window size.
Default: 300.''')

approach = parser.add_mutually_exclusive_group(required=True)
#default=50,
approach.add_argument('-s','--step',type=int, default=False,
                    help='''Specify the step size of sliding window used to generate pairs. default: 50.
Estimate number of read pairs by: numLongReads*(meanLongReadLen - fraglen + 1)/stepsize.
Example: 1.5 million long reads of mean length = 5000, fraglen for PE reads = 300, step size = 50,
will generate around 141 million PE reads.''')

approach.add_argument('-N','--num',type=int, default=False,
                    help='''Specify how many PE reads to sample from each eligible long read.''')

##parser.add_argument('-U','--unique',action="store_true",default=False,
##                    help='''Only has effect when used with -N/--num. This says to only sample unique PE reads from each long
##read. If ''')

parser.add_argument('-r','--readlen',type=int,default=100,
                    help='''specify read length of synthetic paired-end reads. This will be the length each read.
Must be <= fragment length.
Default: 300.''')
parser.add_argument('-S','--start',type=int,default=0,
                    help='''Use this so it starts at position S in each read instead of the first position. This is 0-based.
So, e.g., subtract 1 from 500 to get 499, if you want it to start at the 500th position.''')

parser.add_argument('-R','--random',action='store_true', default=False,
                    help='''Use this in conjunction with --start so that it does not start at the same position for every read.
Instead it will start at a randomly sampled position to start at up to the given distance.
For example, start anywhere in the first 500 bp.
Default: 0 --> Meaning it always starts at the first position.''')

parser.add_argument('--seed',type=int,default=-1,
                    help='''Specify a seed for random sampling (for start site selection above).
This will enable you to do the same sampling time and time again if needed.
Default: -1 --> means a seed will not be used at all.''')

parser.add_argument('-o','--outprefix',type=str,default="out-long2paired",
                    help='Specify output prefix for paired-read files. Default: long2paired-out')

parser.add_argument('-n','--nostats',action="store_true",default=False,
                    help='''By default - # long reads used, # not used, # PE reads generated,
and other stats are put into text files with outprefix name.
Can turn this off.''')

parser.add_argument('-u','--unused',action="store_true",default=False,
                    help='Report names and lengths of unused long reads.')

args = parser.parse_args()


#########
if args.stdin and not (args.fa or args.fq):
    print "When using --stdin, filetype should be explicitly defined --fa or --fq."
    quit()

if args.stdin:
    fastxFile = sys.stdin
    input_gzipped = False
else:
    fastxFile = gzip.open(args.fastx)
    input_gzipped = True
    try:
        fastxFile.next()
        fastxFile.seek(0)
    except IOError:
        input_gzipped = False
        fastxFile = open(args.fastx, 'r')



if args.fa:
    in_fastx = "fasta"
elif args.fq:
    in_fastx = "fastq"
elif args.fastx: ## can figure out from file (right now not from stdin due to non-seekability)
    line1 = fastxFile.next()[0]
    if line1[0] == ">":
        in_fastx = "fasta"
    elif line1[0] == "@":
        in_fastx = "fastq"
    fastxFile.seek(0)
else:
    print "Expected fasta or fastq. File given needs to be reformatted if user thinks it is."
    quit()
if args.ofa:
    out_fastx = "fasta"
elif args.ofq:
    out_fastx = "fastq"
else: #default to in_fastx
    out_fastx = in_fastx


### functions
assert args.readlen <= args.fraglen

def fq_out(record, name, seq, start, end):
    return name+"\n"+seq[start:end]+"\n+\n"+str

## execute

if args.gzip:
    r1 = gzip.open(args.outprefix+"-1."+out_fastx+".gz", 'w')
    r2 = gzip.open(args.outprefix+"-2."+out_fastx+".gz", 'w')
else: ##do no gzip output
    r1 = open(args.outprefix+"-1."+out_fastx, 'w')
    r2 = open(args.outprefix+"-2."+out_fastx, 'w')

unused_readnames = {}
num_used = 0
num_unused = 0
total_long = 0
total_pairs = 0

if args.seed > -1:
    random.seed(args.seed)

for record in SeqIO.parse(fastxFile, in_fastx):
    total_long += 1
##    print record.id
##    print len(record.seq)
    orig_name = record.name
    if args.step:
        if args.random:
            start = random.randint(0, args.start)
        else:
            start = args.start
        if len(record.seq) >= args.fraglen+start:
            num_used += 1
            for i in range(start, len(record.seq)-args.fraglen+1, args.step):
                total_pairs += 1
                record.id = orig_name+":"+str(i)+"-"+str(i+args.fraglen)+"\n"
                record.description = str(args.readlen) ## otherwise, it was giving id and description as name, which was same info twice
                SeqIO.write(record[i:i+args.readlen], r1, out_fastx)
                SeqIO.write(record[i+args.fraglen-args.readlen:i+args.fraglen].reverse_complement(id=True, name=True, description=True), r2, out_fastx)
    ##            SeqIO.write(record[i+args.fraglen-args.readlen:i+args.fraglen], r2, out_fastx)
                ## still need to make it take right end of fragment (right now it just revcomps beginnin)
                ## also -- need to name the subreads of a given long read slightly different names....
        else:
            num_unused += 1
            unused_readnames[record.id] = len(record.seq)
    elif args.num:
        for j in range(args.num):
            total_pairs += 1
            ## Pick random spot
            seqlen = len(record.seq)
            if seqlen >= args.fraglen:
                i = random.randint(0, seqlen-args.fraglen)
                record.id = orig_name+":"+str(i)+"-"+str(i+args.fraglen)+"\n"
                record.description = str(args.readlen)
                SeqIO.write(record[i:i+args.readlen], r1, out_fastx)
                SeqIO.write(record[i+args.fraglen-args.readlen:i+args.fraglen].reverse_complement(id=True, name=True, description=True), r2, out_fastx)
            else:
                num_unused += 1
                unused_readnames[record.id] = len(record.seq)               
                


if not args.stdin:
    fastxFile.close()

r1.close()
r2.close()

### stats/info file(s)
stats = open(args.outprefix+".stats.txt", 'w')
stats.write("total number long reads:\t"+str(total_long)+"\n")
stats.write("number long reads used:\t"+str(num_used)+"\n")
stats.write("number long reads unused:\t"+str(num_unused)+"\n")
stats.write("number pairs made:\t"+str(total_pairs)+"\n")
stats.write("frag lengths:\t"+str(args.fraglen)+"\n")
stats.write("read lengths:\t"+str(args.readlen)+"\n")
stats.close()
if args.unused:
    unused = open(args.outprefix+".unused.txt",'w')
    for r in unused_readnames.keys():
        unused.write(r+"\t"+str(unused_readnames[r])+"\n")
    unused.close()
            



