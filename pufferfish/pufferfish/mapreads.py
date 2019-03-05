from FileListClass import *
import gzip, bz2
import logging, os
logger = logging.getLogger('pufferfish')

##bowtie2 --very-sensitive -N 1 -x sciara -U $reads 2>>$prefix.txt | samtools view -bSh - 2>>$prefix.txt | samtools sort - $prefix 2>>$prefix.txt 
##
##samtools index $prefix.bam 2>>$prefix.txt 


def detect_fastx_type(fastx):
    try:
        assert os.path.exists(fastx)
    except AssertionError as e:
        e.args += (fastx, 'File does not exist')
        raise
    try:
        connection = gzip.open(fastx)
        line = connection.next()
    except IOError:
        connection.close()
        connection = bz2.BZ2File(fastx)
        try:
            line = connection.readline()
        except IOError:
            connection.close()
            connection = open(fastx,'r')
            line = connection.next()
    ans = fastx_check(line)
    connection.close()
    return ans

def fastx_check(line):
    if line[0] == ">":
        return "fa"
    elif line[0] == "@":
        return "fq"


        
def make_bt2(fasta, prefix):
    cmd = 'bowtie2-build ' + fasta + ' ' + prefix
    return cmd, prefix

def bowtie2(bt2, fastx, threads, errfile):
    bowtie2fq = 'bowtie2 --threads %s -q --very-sensitive -N 1 -x %s -U %s 2>%s'
    bowtie2fa = 'bowtie2 --threads %s -f --very-sensitive -N 1 -x %s -U %s 2>%s'
    fx = detect_fastx_type(fastx)
    if fx == 'fa':
        return bowtie2fa % (threads, bt2, fastx, errfile)
    elif fx == 'fq':
        return bowtie2fq % (threads, bt2, fastx, errfile)


## IT APPEARS SAMTOOLS NO LONGER NEEDS THE "-"
## FOR SORT, IT MAY BE BAD TO PUT THERE...
def samtools_view(errfile): # removes unmapped reads
    return 'samtools view -bSh -F 4 - 2>>%s' % (errfile)

##def samtools_sort(bamprefix, errfile):
##    return 'samtools sort - -o %s.bam 2>>%s' % (bamprefix, errfile)
def samtools_sort(bamprefix, threads, errfile):
    return 'samtools sort --threads %s -o %s.bam 2>>%s' % (threads, bamprefix, errfile)


def samtools_index(bamfile, errfile):
    return 'samtools index %s 2>>%s' % (bamfile, errfile)




def map_reads(fastx, bt2, threads, dry=False):
    prefix = os.path.basename(fastx).split(".")[0]
    errfile = prefix + ".mapread.err"
    mapreads = (" | ").join([bowtie2(bt2,fastx,threads,errfile), samtools_view(errfile), samtools_sort(prefix, threads, errfile)])
    bam = prefix + ".bam"
    run_cmd( mapreads, dry )
    cmd = samtools_index(bam,errfile)
    run_cmd( cmd, dry )


def run_cmd(cmd, dry=True):
    if dry:
        print cmd 
    else:
        os.system( cmd )


##def report_commands(parser, args):
##    f = open("pufferfish.mapreads.commands.used.txt",'w')
##    if args.ref_fasta:
##        bt2build,bt2prefix = make_bt2(args.ref_fasta,dry_run=True)
####        f.write(bt2build + "\n")
##        run_cmd(bt2build, f, dry=True)
##    else:
##        bt2prefix = args.bt2
##    for fastx in args.fastxfiles:
##        bt2, sort = map_reads(fastx, bt2prefix, dry_run=True)   
##        f.write(bt2 + "\n" + sort + "\n")
##    f.close()
        
def run(parser, args):
##    report_commands(parser, args)
    
##    if not args.dry:
    if args.ref_fasta:
        prefix = args.ref_fasta.split(".fa")[0]
        cmd, args.bt2 = make_bt2(args.ref_fasta, prefix)
        run_cmd( cmd, args.dry )
    for fastx in args.fastxfiles:
        map_reads(fastx, args.bt2, str(args.threads), args.dry)
    
    
