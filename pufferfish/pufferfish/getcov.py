from FileListClass import *
import logging, os
logger = logging.getLogger('pufferfish')



def bedtools_makewindows(genome, window, step, outbed):
    return "bedtools makewindows -g " + genome + " -w " + window + " -s " + step + " > " + outbed

def coverageBed(bam, windowbed):
    #optionally give path to windowbed made with bedtools_makewindows
    return "coverageBed -abam " + bam + " -b " + windowbed

def sortBed():
    return "sortBed -i - "

def cut(cols=[1,2,3,4]):
    return "cut -f " + (",").join([str(e) for e in cols])

def samtools_view(bam, mapq):
    return 'samtools view -b -h -F 4 -q ' + mapq + " " + bam

def cov_pipeline(bam, mapq, windowbed, bedGraph):
    return (" | ").join([samtools_view(bam,mapq), coverageBed("-", windowbed), sortBed(), cut()]) + " > "  + bedGraph

def windowbed_name(window, step):
    wininfo = "w" + window + ".s" + step
    windowbed = wininfo + ".bed"
    return wininfo, windowbed
    
def bedgraph_name(bam, mapq, wininfo):
    prefix = os.path.basename(bam).split(".")[0]
    bedGraph = prefix + ".q" + mapq + "." + wininfo + ".bedGraph"
    return bedGraph


def report_commands(parser,args):
    f = open("pufferfish.getcov.commands.used.txt",'w')
    wininfo, windowbed = windowbed_name(args.window, args.step)
    makewindows_cmd = bedtools_makewindows(args.genome, args.window, args.step, windowbed)
    f.write(makewindows_cmd + "\n")
    for bam in FileList(args.bams, extension='.bam'):
        bedGraph = bedgraph_name(bam, args.mapq, wininfo)
        cov_cmd = cov_pipeline(bam, args.mapq, windowbed, bedGraph)
        f.write(cov_cmd + "\n")
    rm_cmd = "rm " + windowbed
    f.write(rm_cmd)

def run_commands(parser, args):
    wininfo, windowbed = windowbed_name(args.window, args.step)
    makewindows_cmd = bedtools_makewindows(args.genome, args.window, args.step, windowbed)
    os.system(makewindows_cmd)
    for bam in FileList(args.bams, extension='.bam'):
        bedGraph = bedgraph_name(bam, args.mapq, wininfo)
        cov_cmd = cov_pipeline(bam, args.mapq, windowbed, bedGraph)
        os.system(cov_cmd)
    rm_cmd = "rm " + windowbed
    os.system(rm_cmd)

def run(parser, args):
    assert os.path.exists(args.genome)
    report_commands(parser, args)
    if not args.dry:
        run_commands(parser, args)
        
        
