#!/usr/bin/env python
import sys, os, argparse
import multiprocessing

parser = argparse.ArgumentParser(description = '''
This program is a wrapper over several programs in a pipeline to use Pilon
for genome polishing.
One can start from almost anywhere in the pipeline by specifying the right
flag.

Requires:
bowtie2
samtools/htslib 1.3
picard tools
pilon

NOTE:
This was originally designed to be used with Picard 1.95.
It has been updated to work with 2.1.1.
It can work with either --picardversion.
The semantics changed somewhere between 1.95 and 2.1.1, from:
java -jar MarkDuplicates.jar
to:
java -jar picard.jar MarkDuplicates
So when using --picardversion, specify "--picardversion 1" if your version uses the former and "--picardversion 2" if it uses the latter.
Note:
As of Picard version 2.0.1 (Nov. 2015), Picard requires Java 1.8 (jdk8u66).

Run whole pipeline:
pilon-pipeline.py --pilonjar pathtopilonjar --reference pathtoreferencetobepolished --R1 pathtoR1.fq --R2 pathtoR2.fq  --mkdupjar MarkDuplicates.jar  --picardversion 1
python pilon-pipeline.py --dry --pilonjar pathtopilonjar --reference pathtoreferencetobepolished --R1 pathtoR1.fq --R2 pathtoR2.fq --mkdupjar picard.jar --picardversion 2

''', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--reference', type=str, default=False, required=True,
                    help='''Reference is used in Pilon (to be polished). Also, if a bowtie2 index needs to be built, it will use this.''')

parser.add_argument('--bt2', type=str, default=False,
                    help=''' If a bowtie2 index is already built, specify path to bt2 shared prefix.''')

parser.add_argument('--R1', type=str, default=False,
                    help='''If need to bowtie2 map PE reads to bt2 index,
specify path to file with #1 in pair, "R1"''')

parser.add_argument('--R2', type=str, default=False,
                    help='''If need to bowtie2 map PE reads to bt2 index,
specify path to file with #2 in pair, "R2"''')

parser.add_argument('--bam', type=str, default=False,
                    help='''If already mapped reads, but did not mark duplicates,
use this flag to specify path to the BAM file.''')

parser.add_argument('--mkdup_bam', type=str, default=False,
                    help='''If already mapped reads and marked duplicates,
use this flag and specify path to BAM file with marked duplicates.''')

parser.add_argument('--flagstat', action='store_true', default=False,
                    help='''Use this flag and specify path to desied BAM file, if you want to
to include SAMtools flagstat analysis on MarkDup BAM.''')

parser.add_argument('--dry', action='store_true', default=False,
                    help='''Do not run commands, only print them to screen''')

parser.add_argument('--minins', type=int, default=0,
                    help='''Minimum insert size to expect when mapping paired reads. Default: 0.''')
parser.add_argument('--maxins', type=int, default=1000,
                    help='''Maximum insert size to expect when mapping paired reads. Default: 1000.''')
parser.add_argument('--threads', type=int, default=0,
                    help='''Number of threads to use when mapping paired reads -- also used for samtools sort. Default: auto-detect.''')
parser.add_argument('--stm', type=str, default='768M',
                    help='''Max memory per thread to use when using samtools sort. Default (from samtools): 768M.''')
parser.add_argument('--mkdupjar', type=str, default='', required=True,
                    help='''Path to MarkDuplicates.jar from Picard-Tools.
e.g. /path/to/picard-tools-1.95/MarkDuplicates.jar
This is required.
If skipping this step, just enter anything - it will be ignored.''')
parser.add_argument('--picardversion', type=int, default=2,
                    help = '''Picard version. Provide integer. Options are [1,2].
This was originally designed to be used with Picard 1.95.
It has been updated to work with 2.1.1.
It can work with either --picardversion.
The semantics changed somewhere between 1.95 and 2.1.1, from:
java -jar MarkDuplicates.jar
-->to-->
java -jar picard.jar MarkDuplicates ... 
So when using --picardversion, specify "--picardversion 1" if your version uses the former and "--picardversion 2" if it uses the latter.
Note:
As of Picard version 2.0.1 (Nov. 2015), Picard requires Java 1.8 (jdk8u66).
''')
parser.add_argument('--xmxmem', type=str, default='64g',
                    help='''Provide a string with adjacent integer_letter such as 32g or 64g. Default: 64g.''')
parser.add_argument('--rmdup', action='store_true', default=False,
                    help='''Specify this flag, if you'd like to remove the duplicates rather than just marking them.''')
parser.add_argument('--not_pre_sorted', action='store_true', default=False,
                    help='''Specify this flag if you have given a BAM file that is not sorted.
Usually you will have sorted your BAM file (e.g. as part of mapping the reads).
If you are not supplying a BAM file, and therefore running the pipeline from an earlier step,
then it is automatically sorted.
This flag will mostly not be used as usually BAMs will be pre-sorted.
Using this flag will cause Picard to assume it is not sorted, and that step may take a lot longer.''')

parser.add_argument('--pilonmem', type=str, default='254g',
                    help='''Provide a string with adjacent integer_letter such as 128g or 254g. Default: 254g.''')
parser.add_argument('--pilonjar', type=str, default='', required=True,
                    help='''Path to pilon.jar.
e.g. /users/jurban/software/pilon/1.13.diploid-test/pilon-1.13-10-g51626f4.jar''')
parser.add_argument('--fix', type=str, default='bases',
                    help='''What to fix with Pilon. Default: bases.
              A comma-separated list of categories of issues to try to fix:
                "bases": try to fix individual bases and small indels;
                "gaps": try to fill gaps;
                "local": try to detect and fix local misassemblies;
                "all": all of the above;
                "none": none of the above; new fasta file will not be written.
              The following are experimental fix types:
                "amb": fix ambiguous bases in fasta output (to most likely alternative).
                "breaks": allow local reassembly to open new gaps (with "local").
                "novel": assemble novel sequence from unaligned non-jump reads.''')

parser.add_argument('--clean', action='store_true', default=False,
                    help='''Clean up intermediate files, such as the BAM without duplicates marked. Default: False.''')


parser.add_argument('--slurm', action='store_true', default=False,
                    help='''If using SLURM grid system, then using this flag
will submit each job to the grid consecutively. NOT IMPLEMENTED YET.''')

args = parser.parse_args()

## auto-detect number of cpus available
if args.threads == 0:
    args.threads = multiprocessing.cpu_count()

def run_cmd(cmd, dry=True):
    if dry:
        print cmd
    else:
        os.system( cmd )

def bowtie2_build(fasta, prefix):
    return "bowtie2-build " + fasta + " " + prefix


def bowtie2(R1,R2,index, minins=0, maxins=1000, bt2p=1, stp=1, outbam='reads.bam', stm='768M'):
    return (' ').join(['bowtie2 -p', str(bt2p), '--very-sensitive -N 1 --minins', str(minins), '--maxins', str(maxins), '-x', index, '-1', str(R1), '-2', str(R2), '| samtools sort --threads', str(stp), '-m', str(stm), '>', outbam])

def samtools_index(bam):
    return 'samtools index ' + bam

def samtools_flagstat(bam):
    pre = bam.split(".bam")[0]
    return (" ").join(["samtools flagstat", bam, ">", pre+".flagstats.txt"])

def picard_mark_duplicates1(jar, inbam, outbam, metricsfilename, mem='64g', rmdup=False, presorted=True):
    if rmdup:
        RMDUP="true"
    else:
        RMDUP="false"
    if presorted:
        PRESORTED="true"
    else:
        PRESORTED="false"
    return (" ").join(["java -Xmx" + str(mem), "-jar", jar, "INPUT="+inbam, "OUTPUT="+outbam, "METRICS_FILE="+metricsfilename, "REMOVE_DUPLICATES="+RMDUP, "ASSUME_SORTED="+PRESORTED])

def picard_mark_duplicates2(jar, inbam, outbam, metricsfilename, mem='64g', rmdup=False, presorted=True):
    if rmdup:
        RMDUP="true"
    else:
        RMDUP="false"
    if presorted:
        PRESORTED="true"
    else:
        PRESORTED="false"
    return (" ").join(["java -Xmx" + str(mem), "-jar", jar, "MarkDuplicates INPUT="+inbam, "OUTPUT="+outbam, "METRICS_FILE="+metricsfilename, "REMOVE_DUPLICATES="+RMDUP, "ASSUME_SORTED="+PRESORTED])

def get_picard_fxn(version):
    if version == 1:
        return picard_mark_duplicates1
    elif version == 2:
        return picard_mark_duplicates2

def pilon(jar, genome, infrags, out, mem='254g',fixes='bases'):
    return (" ").join(["java -Xmx"+mem, '-jar', jar, '--genome', genome, '--output', out,  '--changes', '--frags', infrags, '--diploid', '--fix', fixes]) 


def clean(f):
    return "rm " + f


if not args.bt2 and not (args.bam or args.mkdup_bam):
    bt2prefix = os.path.basename(args.reference).split(".")[0]
    cmd = bowtie2_build(args.reference, bt2prefix)
    run_cmd( cmd , args.dry )
elif args.bt2:
    bt2prefix = args.bt2


##if args.bt2 and args.R1 and args.R2 and not args.bam:
if args.R1 and args.R2 and not args.bam:
    ## map reads
    bam = os.path.basename(args.R1).split(".")[0] + ".bam"
    cmd = bowtie2(R1 = args.R1, R2 = args.R2, index = bt2prefix, minins = args.minins, maxins = args.maxins, bt2p = args.threads, stp = args.threads, outbam=bam, stm=args.stm)
    run_cmd( cmd , args.dry )
elif args.bam:
    bam = args.bam
else:
    bam = False

if args.flagstat and bam:
    #flagstat allreads bam file
    cmd = samtools_flagstat(bam)
    run_cmd( cmd , args.dry )
    
    
##if args.bam and not args.mkdup_bam:
if not args.mkdup_bam:
    ## mkdup
    if args.not_pre_sorted:
        presorted=False
    else:
        presorted=True
    picard_mark_duplicates = get_picard_fxn(args.picardversion)
    mkdupbam = bam.split(".bam")[0] + "markdup.bam"
    metricsfilename = bam.split(".bam")[0] + "markdup.metrics.txt"
    cmd = picard_mark_duplicates(jar=args.mkdupjar, inbam=bam, outbam=mkdupbam, metricsfilename=metricsfilename, mem=args.xmxmem, rmdup=args.rmdup, presorted=presorted)
    run_cmd( cmd, args.dry )
    
else:
    #note
    mkdupbam = args.mkdup_bam

if not os.path.exists(mkdupbam + '.bai'):
    samtools_index(mkdupbam)

if args.clean and bam:
    cmd = clean(bam)
    run_cmd( cmd, args.dry )
    if os.path.exists(bam+'.bai'):
        cmd = clean(bam+'.bai')
        run_cmd( cmd, args.dry )

if args.flagstat:
    #flagstat mkdup file
    cmd = samtools_flagstat(mkdupbam)
    run_cmd( cmd , args.dry )

## RUN PILON
##out = args.reference.split(".fa")[0] + ".pilonpolished"
out = os.path.basename(args.reference).split(".fa")[0] + ".pilon"
cmd = pilon(jar = args.pilonjar, genome = args.reference, infrags = mkdupbam, out=out, mem=args.pilonmem,fixes=args.fix)
run_cmd( cmd, args.dry )
    
