from FileListClass import *
import logging, os
logger = logging.getLogger('pufferfish')



def bedtools_makewindows(genome, window, step, outbed):
    return "bedtools makewindows -g " + genome + " -w " + window + " -s " + step + " > " + outbed

## from BEDtools 2.21.0
##def coverageBed(bam, windowbed):
##    #optionally give path to windowbed made with bedtools_makewindows
##    return "coverageBed -abam " + bam + " -b " + windowbed

##updated for BEDtools 2.25.0
def coverageBed(bam, windowbed, genome):
    #optionally give path to windowbed made with bedtools_makewindows
    #assumes BAM and windowBed are sorted
    return "coverageBed -a " + windowbed + " -b " + bam + "-sorted -G " + genome + "-counts"

def sortBed():
    return "sortBed -i - "

def cut(cols=[1,2,3,4]):
    return "cut -f " + (",").join([str(e) for e in cols])

def samtools_view(bam, mapq):
    return 'samtools view -b -h -F 4 -q ' + mapq + " " + bam

##def cov_pipeline(bam, mapq, windowbed, bedGraph):
##    return (" | ").join([samtools_view(bam,mapq), coverageBed("-", windowbed), sortBed(), cut()]) + " > "  + bedGraph

##sortBed not strictly necessary anymore but keeping anyway
def cov_pipeline(bam, mapq, windowbed, genome, bedGraph):
    return (" | ").join([samtools_view(bam,mapq), coverageBed("-", windowbed, genome), sortBed()]) + " > "  + bedGraph


def windowbed_name(window, step):
    wininfo = "w" + window + ".s" + step
    windowbed = wininfo + ".bed"
    return wininfo, windowbed
    
def bedgraph_name(bam, mapq, wininfo):
    prefix = os.path.basename(bam).split(".")[0]
    bedGraph = prefix + ".q" + mapq + "." + wininfo + ".bedGraph"
    return bedGraph

def picard_mark_duplicates(jar, inbam, outbam, metricsfilename, mem='4g', rm_opt_dup=True, rmdup=False, presorted=True):
    if rm_opt_dup:
        RMOPTDUP="true"
    else:
        RMOPTDUP="false"
    if rmdup:
        RMDUP="true"
    else:
        RMDUP="false"
    if presorted:
        PRESORTED="true"
    else:
        PRESORTED="false"
    return (" ").join(["java -Xmx" + str(mem), "-jar", jar, "MarkDuplicates INPUT="+inbam, "OUTPUT="+outbam, "METRICS_FILE="+metricsfilename, "REMOVE_SEQUENCING_DUPLICATE="+RMOPTDUP, "REMOVE_DUPLICATES="+RMDUP, "ASSUME_SORTED="+PRESORTED])


def report_commands(parser,args):
    f = open("pufferfish.getcov.commands.used.txt",'w')
    wininfo, windowbed = windowbed_name(args.window, args.step)
    makewindows_cmd = bedtools_makewindows(args.genome, args.window, args.step, windowbed)
    f.write(makewindows_cmd + "\n")
    for bam in FileList(args.bams, extension='.bam'):
        bedGraph = bedgraph_name(bam, args.mapq, wininfo)
        cov_cmd = cov_pipeline(bam, args.mapq, windowbed, args.genome, bedGraph)
        f.write(cov_cmd + "\n")
    rm_cmd = "rm " + windowbed
    f.write(rm_cmd)

def run_commands(parser, args):
    wininfo, windowbed = windowbed_name(args.window, args.step)
    makewindows_cmd = bedtools_makewindows(args.genome, args.window, args.step, windowbed)
    os.system(makewindows_cmd)
    for bam in FileList(args.bams, extension='.bam'):
        bedGraph = bedgraph_name(bam, args.mapq, wininfo)
        cov_cmd = cov_pipeline(bam, args.mapq, windowbed, args.genome, bedGraph)
        os.system(cov_cmd)
    rm_cmd = "rm " + windowbed
    os.system(rm_cmd)

def run_cmd(cmd, dry=True):
    if dry:
        print cmd 
    else:
        os.system( cmd )
       
def run(parser, args):
    if not args.force:
        assert os.path.exists(args.genome)
##    report_commands(parser, args)
##    if not args.dry:
##        run_commands(parser, args)
##        wininfo, windowbed = windowbed_name(args.window, args.step)
    ## Make windows
    cmd = bedtools_makewindows(args.genome, args.window, args.step, windowbed)
    run_cmd(cmd, args.dry)
    wininfo, windowbed = windowbed_name(args.window, args.step)
    for bam in FileList(args.bams, extension='.bam'):
        
        #filter based on mapq
        if args.mapq > 0:
            tmpbam = "pufferfish.getcov.tmp." + os.path.basename(bam).split(".bam")
            cmd = samtools_view(bam, args.mapq) + " > " + tmpbam
            run_cmd(cmd, args.dry)
        
        #filter out optical and PCR duplicates
        if args.mapq > 0:
            inbam = tmpbam
        else:
            inbam = bam
        if args.filterdup:
            outbam = os.path.basename(tmpbam).split(".bam")[0] + ".mkdup.bam"
            metricsfile = os.path.basename(tmpbam).split(".bam")[0] + ".mkdup.metrics.txt"
            rm_opt_dup = not args.keepopt
            cmd = picard_mark_duplicates(jar=args.filterdup, inbam=tmpbam, outbam=outbam, metricsfilename=metricsfile, mem=args.picardmem, rm_opt_dup=rm_opt_dup, rmdup=args.rmdup, presorted=True)
            run_cmd(cmd, args.dry)
        else:
            outbam = inbam
        
        #get coverage counts in each bin
        bedGraph = bedgraph_name(bam, args.mapq, wininfo)
        cmd = coverageBed(bam=outbam, windowbed=windowbed, genome=args.genome) + " > " + bedGraph
        run_cmd(cmd, args.dry)
##        cov_cmd = cov_pipeline(bam, args.mapq, windowbed, args.genome, bedGraph)
##        os.system(cov_cmd)
    if args.clean:
        #remove bin/window file
        cmd = "rm " + windowbed
        run_cmd(cmd, args.dry)
        
        
