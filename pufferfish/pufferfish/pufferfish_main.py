#!/usr/bin/env python
import os.path
import sys
import argparse

#logger
import logging
logger = logging.getLogger('pufferfish')

# pufferfish imports
#import pufferfish.version
pfv = '0.0.0'


## FUNCTIONS
def run_subtool(parser, args):
    ## the function to be used in subparser.set_defaults(func=func)
    ## it selects the module that the sub-parser is made for
    ## then uses the "run" function in that module
    if args.command == 'mapreads':
        import mapreads as submodule
    elif args.command == 'getcov':
        import getcov as submodule
    elif args.command == 'findpuffs':
        import findpuffs as submodule
    elif args.command == 'dump':
        import pk2txt as submodule
    elif args.command == 'puffcn':
        import cn as submodule
    elif args.command == 'puffcnpy':
        import cnpy as submodule
    elif args.command == 'summits':
        import findsummits as submodule
    elif args.command == 'normalize':
        import normalize as submodule
    elif args.command == 'hmm':
        import generalhmm as submodule
    elif args.command == 'generate':
        import generate as submodule
    elif args.command == 'filter':
        import filterfish as submodule
    elif args.command == 'help':
        import helper as submodule
    # run the chosen submodule
    submodule.run(parser, args)


## This subclass is used as parser_class for submodule sub-parsers
class ArgumentParserWithDefaults(argparse.ArgumentParser):
    ## Child/sub of argparse.ArgumentParser class (the super)
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        ## Add arguments that can be used by all sub-commands
        self.add_argument("-q", "--quiet", help='''Do not output warnings to stderr.''',
                          action="store_true",
                          dest="quiet")



    
    
def main():
    logging.basicConfig()

    ## Create top-level parser
    parser = argparse.ArgumentParser(prog="pufferfish", description=''' PufferFish - HMM-based approach(es) to finding and analyzing developmentally regulated amplicons (genomic sites that are programmed to increase in copy number over time).''',
                                     formatter_class=argparse.RawTextHelpFormatter)#ArgumentDefaultsHelpFormatter
    parser.add_argument('-v', '--version', help='''Installed pufferfish version.''',
                        action='version',
                        version='%(prog)s ' + pfv)#str(pufferfish.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command',
                                       parser_class=ArgumentParserWithDefaults)


    ## Create a sub-command parser for mapreads
    parser_mapreads = subparsers.add_parser('mapreads',
                                            help='''Depends on Bowtie2 and SAMtools.
Maps fasta/fastq files to genome (can provide bt2 index or fasta reference (which will first be converted to bt2 index).
Maps reads to genome and filters out unmapped reads before sorting and indexing.''')
    parser_mapreads.add_argument('fastxfiles', metavar='fastxfiles', nargs='+', type=str,
                                 help='''Paths to as many fasta/fastq files as youd like to map. Can be gz or bz2.''')
    parser_mapreads_reftype = parser_mapreads.add_mutually_exclusive_group(required=True)
    parser_mapreads_reftype.add_argument('-b', '--bt2', type=str,
                               help='''Path to bt2 index prefix.''')
    parser_mapreads_reftype.add_argument('-r', '--ref_fasta', type=str,
                               help='''Path to reference fasta. A bt2 index will be generated in same dir.''')
    parser_mapreads.add_argument('--threads', '-p', type=int, default=1, help='''Number of parallel threads for bowtie2 to use. Default: 1.''')

    parser_mapreads.add_argument('--dry', action='store_true', default=False, help='''Only writes out the commands that will be used if set.''')

    parser_mapreads.set_defaults(func=run_subtool)


##    ## Create a sub-command parser for filterdup
##    parser_mapreads = subparsers.add_parser('filterdup',
##                                            help='''Depends on Picard Tools 2.1.1, BEDtools, pybedtools, pysam.
##Remove optical duplicates and marks PCR duplicates.
##All PCR duplicates except K at a given site are removed.
##K is determined by a binomial calculuation using a bin size and number of reads in a given bin.
##Then any duplicates in that bin are subject to filtering down to K.
##1. Remove optical duplicates and mark PCR duplicates.
##2. Make bins
##3. Get read count in those bins
##4. For each bin, check if there are marked reads. If so, calculate K and filter.
##5. write remaining reads as you go...
##''')
##    parser_filterdup.add_argument('bams', metavar='bams', nargs='+',
##                               type=str,
##                               help=''' Paths to BAM files that need duplicate filtering.''')
##    parser_filterdup.add_argument('-g', '--genome', type=str,
##                               help='''Path to BEDtools genome file describing reference reads were mapped to.''')
##    parser_filterdup.add_argument()
##    parser_filterdup.add_argument('--dry', action='store_true', default=False, help='''Only writes out the commands that will be used if set.''')
##
##    parser_filterdup.set_defaults(func=run_subtool)
##    

    ## Create sub-command parser for getcov
    ## TODO add filterdup possibility from macs2... rm pcr dups
    parser_getcov = subparsers.add_parser('getcov',
                                          help=''' Depends on SAMtools and BEDtools.''')
    parser_getcov.add_argument('bams', metavar='bams', nargs='+',
                               type=str,
                               help=''' Paths to as many bam files as you need to get coverage for.
Can include a file-of-filenames (FOFN) and tarballs as well.''')
    
    parser_getcov.add_argument('-f','--filterdup', type=str,
                               help='''Provide /path/to/picard.jar  (picard v2.1.1 or higher)''')
    parser_getcov.add_argument('-g', '--genome', type=str, required=True,
                               help='''Path to file.genome as needed and defined by BEDTools. See "bedtools makewindows" or "bedtools coverage"''')
    parser_getcov.add_argument('-w', '--window', type=str, default='500',
                               help='''Integer window size - will be counting in windows of this size. Default: 500.''')
    parser_getcov.add_argument('-s', '--step', type=str, default='500',
                               help='''Integer step size - will slide window over this much. Default: 500.''')
    parser_getcov.add_argument('-Q', '--mapq', type=str, default='0',
                               help='''Integer mapq cut-off - only include reads with mapping quality >= INT. Default: 0.''')
    parser_getcov.add_argument('-m', '--picardmem', type=str, default='4g',
                               help='''Provide memory needed/available to Picard MarkDuplicates as integer_letter string, such as 500m, 1g, 2g, 64g, etc. Default: 4g.''')
    parser_getcov.add_argument('--keepopt',action='store_true', default=False, help='''Optical duplicates are removed by default. This flag says to mark them instead.''')
    parser_getcov.add_argument('--rmdup',action='store_true', default=False, help='''PCR duplicates are marked by default. This flag will result in removing them (as well as optical duplicates).''')
    parser_getcov.add_argument('--dry',action='store_true', default=False, help='''Only writes out the commands that will be used if set.''')
    parser_getcov.add_argument('--clean',action='store_true',default=False,help='''Remove intermediate files... Default: False.''')
    parser_getcov.add_argument('--force',action='store_true',default=False,help='''Ignore assertions. Good for made-up filenames when debugging in dry-runs. Do not use this for real run. Default: False.''')
    parser_getcov.set_defaults(func=run_subtool)



    ## Create sub-command parser for findpuffs
    parser_findpuffs = subparsers.add_parser('findpuffs',
                                             help='''Take in getcov bedGraphs, do stuff.''')
    parser_findpuffs_input = parser_findpuffs.add_mutually_exclusive_group(required=True)
    parser_findpuffs_input.add_argument('-i','--input', type=str,
                                  help='''Input file -- a tab-sep file with 2 columns (stage number, filename) with a single line for all getcov bedGraph files you wish to include.
Example:
1\tstage1.bedGraph''')
    parser_findpuffs_input.add_argument('-ip','--inpickle', type=str,
                                        help='''Pickle file (e.g. data.pk) containing already processed getcov bedGraphs as MultiCovBed object.''')

    parser_findpuffs.add_argument('-op','--outpickle', type=str, default='data.fp.pk',
                                        help='''Name for output pickle file (e.g. data.fp.pk.gz) that will contain the MultiCovBed object made when this is run.
Pickled data is automatically gzipped. If .gz not at end of given filename, it will be added.
If the filename exists it will be erased and written over.
Default: data.fp.pk.gz''')
    
    parser_findpuffs.add_argument('-s1', '--smoothbeforemedian', default=False, type=int,
                                  help='''Smooth counts in bins before finding the median bin counts (and before median normalization).
Must provide integer window size to smooth in (should be longer than bin size)
-- e.g. if bin size = 500, smoothing bandwidth could be 10000.
One probably should not do both --smoothbeforemedian and --smoothbeforecor. Pick one or none.''')
    parser_findpuffs.add_argument('-s2', '--smoothbeforecor', default=False, type=int,
                                  help='''Smooth counts in bins after median normalization, but before finding correlations in each bin.
Must provide integer window size to smooth in (should be longer than bin size)
-- e.g. if bin size = 500, smoothing bandwidth could be 10000.
One probably should not do both --smoothbeforemedian and --smoothbeforecor. Pick one or none.''')

    parser_findpuffs.add_argument('-bw', '--corsmoothbandwidth', type=int, default=15000,
                                  help='''For smoothing correlation scores. Provide integer window size to smooth in (should be longer than bin size)
-- e.g. if bin size = 500, smoothing bandwidth could be 10000. Default: 15000.''')

    parser_findpuffs.add_argument('-mep', '--max_empty_bin_pct', type=float, default=0.4,
                                  help='''For filtering contigs out, contig is allowed to have up to X (proportion between 0 and 1) bins with 0 coverage. Default: 0.4.''')

    parser_findpuffs.add_argument('-mop', '--max_offending_samples_pct', type=float, default=0.4,
                                  help='''For filtering contigs out, contig is allowed to exceed max_empty_bin_pct in Y (proportion between 0 and 1) of the samples. Default: 0.4.''')


    parser_findpuffs.set_defaults(func=run_subtool)
    

    ## Create sub-command for dump
    parser_dump = subparsers.add_parser('dump',
                                                help=''' Take in pickled object containing MultiCovBed object where statepath has been found.
Output bedGraph of statepath and/or BED file containing coordinates of states.''')
    parser_dump.add_argument('-ip', '--inpickle', type=str, required=True,
                                    help='''Path to input pickle.''')
    parser_dump.add_argument('-p','--prefix', type=str, required=True,
                                     help='''Output files with provided --prefix''')
    parser_dump.add_argument('-c', '--counts', default=False, action='store_true',
                             help='''Output counts as expanded single-step, bedGraph-like file. Counts will be normalized.
Columns 4+ will have counts from files in order files were given.
Can use this with awk to create bedGraphs for each -- e.g. for i in {4..8}; do awk -v "i=$i" 'OFS="\t" {print $1,$2,$3,$i}' q30.w500.s500.default.counts.bedGraph > q30.w500.s500.default.counts.$i.bedGraph; done ''')
    parser_dump.add_argument('-fc', '--filtered_counts', default=False, action='store_true',
                             help='''Output counts from filtered contigs as expanded single-step, bedGraph-like file. Counts will likely NOT be normalized (as filtereing is done prior to median normalization).
Columns 4+ will have counts from files in order files were given.
Can use this with awk to create bedGraphs for each -- e.g. awk 'OFS="\t" {print $1,$2,$3,$4}' filtered_counts.bedGraph > file.filtered_counts.bedGraph ''')
    parser_dump.add_argument('-vc','--viterbi_collapsed', default=False, action='store_true',
                                     help='''Output viterbi statepath expanded single-step bedGraph with provided --prefix''')
    parser_dump.add_argument('-ve','--viterbi_expanded', default=False, action='store_true',
                                     help='''Output viterbi statepath collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-cc','--correlations_collapsed', default=False, action='store_true',
                                     help='''Output bin correlation scores as collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-ce','--correlations_expanded', default=False, action='store_true',
                                     help='''Output bin correlation scores as expanded single-step bedGraph with provided --prefix''')
    parser_dump.add_argument('-scc','--smoothed_correlations_collapsed', default=False, action='store_true',
                                     help='''Output smoothed bin correlation scores as collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-sce','--smoothed_correlations_expanded', default=False, action='store_true',
                                     help='''Output smoothed bin correlation scores as expanded single-step bedGraph with provided --prefix''')
    parser_dump.add_argument('-sc','--slopes_collapsed', default=False, action='store_true',
                                     help='''Output bin slopes as collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-se','--slopes_expanded', default=False, action='store_true',
                                     help='''Output bin slopes as expanded single-step bedGraph with provided --prefix''')

    parser_dump.set_defaults(func=run_subtool)


    ## create sub-sommand for puffcn (puff copy number)
    parser_puffcn = subparsers.add_parser('puffcn',
                                          help = '''Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample
(for additional Fold-enrichment normalization), define whether a region is best explained by cn=1,2,4,8,16,32,64.
Ideally, this can give an idea of where replication forks approximately reach from each firing.''')
    parser_puffcn.add_argument('-l','--latestage', type=str, required=True,
                                help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_puffcn.add_argument('-e','--earlystage', type=str, required=False, default=False,
                               help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')
    parser_puffcn.add_argument('--replace', action='store_true', default=False,
                                        help='''Turn on "replace" functionality. By default this will replace '.' in the count column of bedGraphs with '0'.
Use --replace_with and --replace_this to change.''')
    parser_puffcn.add_argument('--replace_this', type=str, default='.',
                                        help='''Used with --replace. Specify the character in count column to replace. Default = '.' ''')
    parser_puffcn.add_argument('--replace_with', type=str, default='0',
                                        help='''Used with --replace. Specify the character to replace the --replace_this character with.
Must be a string that can be converted to a float. Default = '0' ''')
    parser_puffcn_protocol = parser_puffcn.add_mutually_exclusive_group(required=True)
    parser_puffcn_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_puffcn_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then they are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_puffcn_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then they are smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available.
Then the HMM is run.
Note: if early is not present, this is same as protocol 4.''')
    parser_puffcn_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth.
Then the HMM is run.
Note: if early is not present, this is same as protocol 3.''')

    parser_puffcn_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth)
Then the HMM is run.
Note: if early is not present, this is same as protocol 6.''')

    parser_puffcn_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available. (i.e. smooth -> L/E)
Then the HMM is run.
Note: if early is not present, this is same as protocol 5.''')
    parser_puffcn_protocol.add_argument('-7', '--protocol7', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then the HMM is run.
Note: No median normalization or smoothing is performed. If only late is given, then this is just an identity/pass-through function.''')
    
    parser_puffcn.add_argument('-m', '--emodel', type=str, default='normal',
                               help='''Specify emissions model to assume for HMM. Options: normal, exponential, poisson, geometric, gamma. Default: normal.
Note that all you ever need to do is given the expected means and standard deviations for each state from sampled data.
The normal model will use those directly. Poisson will use the means as lambda. Exponential will use the means as B w/ rate 1/B. Geometric will also use the means as 1/mu.
Gamma will estimate alpha and beta (shape and scale) parameters from the means and standard deviations - if you have A/B in mind, convert to mu and sigma by A*B and (A*B^2)^0.5.
Note that the exponential is same as the gamma when shape is set to 1, making mu and sigma equal B. Thus, if you do gamma w/ muscale=1, then you should get same as exponential.''')
    parser_puffcn.add_argument('-p', '--path', type=str, default='viterbi',
                               help='''Specify whether to take state path defined by viterbi or posterior decoding. Options: viterbi, posterior. Default: viterbi.''')
    parser_puffcn.add_argument('-s', '--scale', action='store_true', default=False,
                               help='''Before going back into the HMM, re-scale counts back towards original magnitude by multiplying by median.
This will also scale the HMM state means by the median to cooperate.''')
    parser_puffcn.add_argument('-c', '--collapsed', action='store_true', default=False,
                               help='''Return collapsed variable-step bedGraph instead of expanded single-step bedGraph.
This is often a much smaller file.''')
    parser_puffcn.add_argument('-ps', '--pseudo', type=float, default=0.1,
                               help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
Should be between 0 and 1.
Should be small enough to not change other values much,
but big enough such that numbers divided by 0+pseudo do not become massive.
Default: 0.1.''')
    parser_puffcn.add_argument('-bw', '--bandwidth', type=int, default=2500,
                               help=''' If kernel smoothing, specify bandwidth (int).
Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
Default: 2500.''')
    parser_puffcn.add_argument('--impute', type=int, default=False,
                               help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
This bandwidth is generally longer than the one you would provide for regular smoothing.
Only bins with a count of 0 will take on smoothed (imputed) values.
Try: 10000.
NOTE: In practice, this lead to its own set of problems and I do not recommend using it in its current form.''')
    parser_puffcn.add_argument('--counts', type=str, default=False,
                               help=''' Use this flag and specify an output prefix for the final normalized late stage bin counts bedGraph.''')
    parser_puffcn.set_defaults(func=run_subtool)
    parser_puffcn.add_argument('--levels', action='store_true', default=False,
                               help=''' Use this flag to output levels bdg instead of state bdg.
Levels are the means 1,2,4,8,16,32,64.
These are obtained by 2**(state-1).''')


    parser_puffcn.add_argument('--mu', type=str, default='1,2,4,8,16,32,64',
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: 1,2,4,8,16,32,64
To change: Provide comma-seprated list of state means.
The number of states will be calculated from this list.
If changing state sigmas (used in normal model), it must have same number of states represented.

NOTE: If using exponential or geometric distribution, provide the expected mean RCN values of the states
    as you would for normal or poisson models. This script will automatically take their inverses to work
    in the exponential and geometric models.''')

    parser_puffcn.add_argument('--sigma', type=str, default=None,
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state sigmas (stdevs) were previously hard-coded.
This option allows some flexibility from the command-line to change the state sigmas.
Default: if not changed, defaults to square root of state means (Poisson-like).
To change: Provide comma-seprated list of state sigmas.
Alternatively: Use --mu_scale (default False) with a scaling factor multiplied against the MUs.
The number of states is calculated from this state mean list, which defaults to 7.
If changing state sigmas (used in normal model), it must have same number of states represented as state means.''')

    parser_puffcn.add_argument('--mu_scale', type=float, default=None,
                               help=''' See --sigma for more details on sigmas.
Use this to scale means (--mu) to use as stdevs (sigma) instead of taking square roots of means.
For example, --mu_scale 0.5 will use mu*0.5 as the stdev.''')

##    parser_puffcn.add_argument('--changestate', type=float, default=0.001,
##                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
##The default transition probabilities were previously hard-coded.
##For now, there are two parameters for transition probabilities: changing states or staying in a state.
##They are shard by all states.
##This option allows some flexibility from the command-line to change the probability of changing states.
##Default: 0.001.
##NOTE: To ensure transition probabilities from state i to j sum to 1,
##    the final transition probabilities will actually be X/sum(all trans i to j),
##    where X is 0.001 by default or what user has given.
##''')
##
##    parser_puffcn.add_argument('--samestate', type=float, default=0.999,
##                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
##The default transition probabilities were previously hard-coded.
##For now, there are two parameters for transition probabilities: changing states or staying in a state.
##They are shard by all states.
##This option allows some flexibility from the command-line to change the probability of staying in the current state (not changing).
##Default: 0.999.
##NOTE: To ensure transition probabilities from state i to j sum to 1,
##    the final transition probabilities will actually be X/sum(all trans i to j),
##    where X is 0.999 by default or what user has given.
##''')


    parser_puffcn.add_argument('--special_idx', type=int, default=0,
                               help='''Only for use if you're very familiar with the program (and change defaults).
The default state means is 1,2,4,8,16,32,64.
The default index for the mean that represents copy number 1 is 0.
In this lingo - CN=1 is the special state, and the 0-based index of the special state in that list is 0.
If you were to change parameters that affect where the special state is in the list, make sure to change this index.
This index is only used to help construct initial probabilities and transition probabilies.
If understood, it can be used to designate any single special state (not necessarily the one that corresponds to CN=1).
The other parameters to use with this are:
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--prob_leave_special
--prob_stay_special
--prob_other_to_special
--prob_other_to_other
--prob_other_to_self

Alternative to, an initial probability vector can be given with --initialprobs 
''')

    parser_puffcn.add_argument('--init_special', type=float, default=0.997,
                               help='''Probability of starting in the 'special state' (usually copy number = 1). Default: 0.997.
The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
''')



    parser_puffcn.add_argument('--leave_special_state', type=float, default=0.001,
                               help='''Probability of leaving the 'special state' (usually copy number = 1).
Default: 0.001.
If number is betwen 0 and 1, it will be assumed a probability.
If number given is > 1, then it will be treated as the average length (number of bins) of the special state.
For example, if 1000 is given, it will be 1/1000 = 0.001.
In terms of bp lengths, one would need to multiply n_bins * bin_length OR divide bp length by bin_length
Thus, if you want to see a change every 500 kb w/ 500 bp bins, then 500kb/500 = 1 kb = 1000 -- which will be interpreted as 0.001.
Or as another example, if you expect to see a change every 2 Mb with 100 bp bins, then 2e6/1e2 = 1e4 = 10 kb = 10000, interpreted as 0.0001.

The probability of staying in this state is the complement: 1-p
''')

    parser_puffcn.add_argument('--leave_other', type=str, default=None,
                               help='''Probability of leaving one of the other states.
This defaults to --leave_special_state making all transition probabilities out of states the same (0.001 by default).

To change, provide a probability of leaving (p).

If the  first number is betwen 0 and 1, it will be assumed a probability.
If the first number given is > 1, then it will be treated as the average length (number of bins).
For example, if 1000 is given, it will be 1/1000 = 0.001.

If only 1 number is given, then that is assumed to be the probability of transitioning to all the other states.

You can also give a comma-separated pair of 2 probabilities:
    prob of leaving to special state
    prob of leaving to another 'non-special' state.
Make sure the probabilities sum to what you expect the overall probability of leaving the state is...
    which should be p_to_special + p_to_nonspecial * (num_non_special-1) = p_to_special + p_to_nonspecial (nstates-2)

For example, in a 7-state model:
    0.001,0.0002 --> 0.001 + 0.0002 * 5 = 0.001 + 0.001 = 0.002
    OR
    0.001,0.001 --> 0.001 + 0.001 * 5 = 0.006

If the second number is > 1, the same rules apply as to the first number.

For other analyses, I've used:
0.00001,0.000000000001
OR
0.001,0.0000000001

The probability of staying in these states is the complement: 1-p1-p2

NOTE: the program forces the transition probabilities of a given state to sum to 1.
''')

    parser_puffcn.add_argument('--initialprobs', type=str, default=None,
                               help='''PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: [0.997, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]
The default will change with more or less states described w/ --mu and --sigma.
By default, the first state will start out as 0.997 as above, all other states will be (1-0.997)/n_other_states.
That behavior also changes with following parameters:
--special_idx -- determines which state (not necessarily first) will be given default 0.997 (OR other with --initcn1)
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--leave_special_state

--prob_other_to_special
--prob_other_to_other
--prob_other_to_self
To change the initial probs manually: Provide comma-separated list of initial probs -e.g.: '0.997,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005'
This must have same number of states represented as state means (--mu; default 7).

''')
##'0.997,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005'  [0.997, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]

    parser_puffcn.add_argument('--kmeans', type=int, default=None,
                               help='''PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default HMM parameters were previously hard-coded.
There are now other options available to tune the paramters.
This option is the first way to learn parameters.
It uses kmeans clustering of the data to estimate initial, transition, and emission probs.
For now, you need to make an assumption about k by providing an integer.
This model probably works best when you expect 2-3 states...
This option over-rides all other parameter options (which will be ignored).
''')

    parser_puffcn.set_defaults(func=run_subtool)


    ## create sub-command for summits
    parser_summits = subparsers.add_parser('summits',
                                           help=''' Find summits...''')
    parser_summits.add_argument('-l','--latestage', type=str, required=True,
                                help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_summits.add_argument('-e','--earlystage', type=str, required=False, default=False,
                               help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')
    parser_summits_protocol = parser_summits.add_mutually_exclusive_group(required=True)
    parser_summits_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_summits_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then they are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_summits_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then they are smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available.
Then the HMM is run.
Note: if early is not present, this is same as protocol 4.''')
    parser_summits_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth.
Then the HMM is run.
Note: if early is not present, this is same as protocol 3.''')

    parser_summits_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth)
Then the HMM is run.
Note: if early is not present, this is same as protocol 6.''')

    parser_summits_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available. (i.e. smooth -> L/E)
Then the HMM is run.
Note: if early is not present, this is same as protocol 5.''')
    
    parser_summits.add_argument('-ps', '--pseudo', type=float, default=0.1,
                               help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
Should be between 0 and 1.
Should be small enough to not change other values much,
but big enough such that numbers divided by 0+pseudo do not become massive.
Default: 0.1.''')
    parser_summits.add_argument('-bw', '--bandwidth', type=int, default=2500,
                               help=''' If kernel smoothing, specify bandwidth (int).
Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
Default: 2500.''')
    parser_summits.add_argument('--impute', type=int, default=False,
                               help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
This bandwidth is generally longer than the one you would provide for regular smoothing.
Only bins with a count of 0 will take on smoothed (imputed) values.
Try: 10000.''')
    
    parser_summits_regions = parser_summits.add_mutually_exclusive_group(required=True)

    parser_summits_regions.add_argument('--regions', type=str, default=False,
                                help = ''' Find summits in these regions - provide BED file.''')
    parser_summits_regions.add_argument('--states', type=str, default=False,
                                        help=''' Provide statepath bedGraph output by "cn" sub-command. Peak regions will be found automatically.''')
    parser_summits.add_argument('--thresh_state', type=int, default=1,
                                help=''' Used with --states. Only consider regions with states higher than state given. Default: 1.''')
    parser_summits.add_argument('--merge1', type=float, default=10e3,
                                help = '''Used with --states. After extracting only bins with higher state value than --thresh_state, merge bins if they are with in --merge1 bp from each other. Default: 10e3.''')
    parser_summits.add_argument('--minwidth', type=float, default=50e3,
                                help = '''After extracting bins with states > --thresh_state and merging remaining bins that are within --merge1 bp of each other,
only keep merged regions > --minwidth.''')
    parser_summits.add_argument('--merge2', type=float, default=40e3,
                                help = '''After (i) extracting bins with states > --thresh_state, (ii) merging remaining bins that are within --merge1 bp of each other,
(iii) retaining only merged regions > --minwidth, merge regions that are within --merge2 bp of each other.''')
    parser_summits.add_argument('--max_state_thresh', type=int, default=2,
                                help = ''' After (i) extracting bins with states > --thresh_state, (ii) merging remaining bins that are within --merge1 bp of each other,
(iii) retaining only merged regions > --minwidth, (iv) merging filtered regions that are within --merge2 bp of each other,
only retain the remaining merged regions if their maximum state is > --max_State_thresh''')

    
    parser_summits.set_defaults(func=run_subtool)


    ## create sub-command for normalize
    parser_normalize = subparsers.add_parser('normalize',
                                          help = '''Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample
(for additional Fold-enrichment normalization), just return the late-stage sample with normalized values as specified by protocol options below.''')

    parser_normalize.add_argument('-l','--latestage', type=str, required=True,
                                help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_normalize.add_argument('-e','--earlystage', type=str, required=False, default=False,
                               help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')
    parser_normalize.add_argument('--replace', action='store_true', default=False,
                                        help='''Turn on "replace" functionality. By default this will replace '.' in the count column of bedGraphs with '0'.
Use --replace_with and --replace_this to change.''')
    parser_normalize.add_argument('--replace_this', type=str, default='.',
                                        help='''Used with --replace. Specify the character in count column to replace. Default = '.' ''')
    parser_normalize.add_argument('--replace_with', type=str, default='0',
                                        help='''Used with --replace. Specify the character to replace the --replace_this character with.
Must be a string that can be converted to a float. Default = '0' ''')

    parser_normalize_protocol = parser_normalize.add_mutually_exclusive_group(required=True)
    parser_normalize_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_normalize_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then they are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_normalize_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then they are smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available.
Then the HMM is run.
Note: if early is not present, this is same as protocol 4.''')
    parser_normalize_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth.
Then the HMM is run.
Note: if early is not present, this is same as protocol 3.''')

    parser_normalize_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth)
Then the HMM is run.
Note: if early is not present, this is same as protocol 6.''')

    parser_normalize_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available. (i.e. smooth -> L/E)
Then the HMM is run.
Note: if early is not present, this is same as protocol 5.''')
    
    parser_normalize.add_argument('-c', '--collapsed', action='store_true', default=False,
                               help='''Return collapsed variable-step bedGraph instead of expanded single-step bedGraph.
This is often a much smaller file.''')
    
    parser_normalize.add_argument('-ps', '--pseudo', type=float, default=0.1,
                               help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
Should be between 0 and 1.
Should be small enough to not change other values much,
but big enough such that numbers divided by 0+pseudo do not become massive.
Default: 0.1.''')
    
    parser_normalize.add_argument('-bw', '--bandwidth', type=int, default=2500,
                               help=''' If kernel smoothing, specify bandwidth (int).
Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
Default: 2500.''')
    parser_normalize.add_argument('--impute', type=int, default=False,
                               help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
This bandwidth is generally longer than the one you would provide for regular smoothing.
Only bins with a count of 0 will take on smoothed (imputed) values.
Try: 10000.''')
        
##    parser_normalize.add_argument('--counts', type=str, default=False,
##                           help=''' Use this flag and specify an output prefix for the final normalized late stage bin counts bedGraph.''')


    parser_normalize.set_defaults(func=run_subtool)








###########GENERATE
        ## create sub-sommand for puffcn (puff copy number)
    parser_generate = subparsers.add_parser('generate',
                                          help = '''Generate emitted_data and statepath bedGraphs.''')

    parser_generate.add_argument('-f','-b', '-i', '--bedgraph', type=str, required=True,
                                help='''Provide path to bedGraph that contains the intervals in first 3 columns to return with generated data.''')


    parser_generate.add_argument('-m', '--emodel', type=str, default='normal',
                               help='''Specify emissions model to assume for HMM. Options: normal, exponential. Default: normal.''')

    parser_generate.add_argument('--mu', type=str, default='1,2,4,8,16,32,64',
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: 1,2,4,8,16,32,64
To change: Provide comma-seprated list of state means.
The number of states will be calculated from this list.
If changing state sigmas (used in normal model), it must have same number of states represented.

NOTE: If using exponential or geometric distribution, provide the expected mean RCN values of the states
    as you would for normal or poisson models. This script will automatically take their inverses to work
    in the exponential and geometric models.''')

    parser_generate.add_argument('--sigma', type=str, default=None,
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state sigmas (stdevs) were previously hard-coded.
This option allows some flexibility from the command-line to change the state sigmas.
Default: if not changed, defaults to square root of state means (Poisson-like).
To change: Provide comma-seprated list of state sigmas.
Alternatively: Use --mu_scale (default False) with a scaling factor multiplied against the MUs.
The number of states is calculated from this state mean list, which defaults to 7.
If changing state sigmas (used in normal model), it must have same number of states represented as state means.''')

    parser_generate.add_argument('--mu_scale', type=float, default=None,
                               help=''' See --sigma for more details on sigmas.
Use this to scale means (--mu) to use as stdevs (sigma) instead of taking square roots of means.
For example, --mu_scale 0.5 will use mu*0.5 as the stdev.''')


    parser_generate.add_argument('--special_idx', type=int, default=0,
                               help='''Only for use if you're very familiar with the program (and change defaults).
The default state means is 1,2,4,8,16,32,64.
The default index for the mean that represents copy number 1 is 0.
In this lingo - CN=1 is the special state, and the 0-based index of the special state in that list is 0.
If you were to change parameters that affect where the special state is in the list, make sure to change this index.
This index is only used to help construct initial probabilities and transition probabilies.
If understood, it can be used to designate any single special state (not necessarily the one that corresponds to CN=1).
The other parameters to use with this are:
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--prob_leave_special
--prob_stay_special
--prob_other_to_special
--prob_other_to_other
--prob_other_to_self

Alternative to, an initial probability vector can be given with --initialprobs 
''')

    parser_generate.add_argument('--init_special', type=float, default=0.997,
                               help='''Probability of starting in the 'special state' (usually copy number = 1). Default: 0.997.
The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
''')



    parser_generate.add_argument('--leave_special_state', type=float, default=0.001,
                               help='''Probability of leaving the 'special state' (usually copy number = 1).
Default: 0.001.
If number is betwen 0 and 1, it will be assumed a probability.
If number given is > 1, then it will be treated as the average length (number of bins) of the special state.
For example, if 1000 is given, it will be 1/1000 = 0.001.
In terms of bp lengths, one would need to multiply n_bins * bin_length OR divide bp length by bin_length
Thus, if you want to see a change every 500 kb w/ 500 bp bins, then 500kb/500 = 1 kb = 1000 -- which will be interpreted as 0.001.
Or as another example, if you expect to see a change every 2 Mb with 100 bp bins, then 2e6/1e2 = 1e4 = 10 kb = 10000, interpreted as 0.0001.

The probability of staying in this state is the complement: 1-p
''')

    parser_generate.add_argument('--leave_other', type=str, default=None,
                               help='''Probability of leaving one of the other states.
This defaults to --leave_special_state making all transition probabilities out of states the same (0.001 by default).

To change, provide a probability of leaving (p).

If the  first number is betwen 0 and 1, it will be assumed a probability.
If the first number given is > 1, then it will be treated as the average length (number of bins).
For example, if 1000 is given, it will be 1/1000 = 0.001.

If only 1 number is given, then that is assumed to be the probability of transitioning to all the other states.

You can also give a comma-separated pair of 2 probabilities:
    prob of leaving to special state
    prob of leaving to another 'non-special' state.
Make sure the probabilities sum to what you expect the overall probability of leaving the state is...
    which should be p_to_special + p_to_nonspecial * (num_non_special-1) = p_to_special + p_to_nonspecial (nstates-2)

For example, in a 7-state model:
    0.001,0.0002 --> 0.001 + 0.0002 * 5 = 0.001 + 0.001 = 0.002
    OR
    0.001,0.001 --> 0.001 + 0.001 * 5 = 0.006

If the second number is > 1, the same rules apply as to the first number.

For other analyses, I've used:
0.00001,0.000000000001
OR
0.001,0.0000000001

The probability of staying in these states is the complement: 1-p1-p2

NOTE: the program forces the transition probabilities of a given state to sum to 1.
''')

    parser_generate.add_argument('--initialprobs', type=str, default=None,
                               help='''PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: [0.997, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]
The default will change with more or less states described w/ --mu and --sigma.
By default, the first state will start out as 0.997 as above, all other states will be (1-0.997)/n_other_states.
That behavior also changes with following parameters:
--special_idx -- determines which state (not necessarily first) will be given default 0.997 (OR other with --initcn1)
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--leave_special_state

--prob_other_to_special
--prob_other_to_other
--prob_other_to_self
To change the initial probs manually: Provide comma-separated list of initial probs -e.g.: '0.997,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005'
This must have same number of states represented as state means (--mu; default 7).

''')

    
    parser_generate.set_defaults(func=run_subtool)





    ## create sub-command for filter/filterfish
    parser_filter = subparsers.add_parser('filter',
                                       help = '''Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample
    (for additional Fold-enrichment normalization), just return the late-stage sample with normalized values as specified by protocol options below.''')

    parser_filter.add_argument('--counts', type=str, default=False,
                            help=''' Use this flag and specify an output prefix for the final normalized late stage bin counts bedGraph.''')

    parser_filter_unit = parser_filter.add_mutually_exclusive_group()
    parser_filter_unit.add_argument('-sd1','--stdev_above', action='store_true', default=False, 
                             help='''Use value given as multiple of standard deviations above the mean.''')
    parser_filter_unit.add_argument('-sd2','--stdev_below', action='store_true', default=False, 
                             help='''Use value given as multiple of standard deviations BELOW the mean.''')
    parser_filter_unit.add_argument('-mu','--mean', action='store_true', default=False, 
                             help='''Use value given as multiple of the mean.''')
    parser_filter.add_argument('-V','--value', type=float, required=True,
                             help='''Value to filter on -- a float. Required.''')
    parser_filter.add_argument('-R','--relation', type=str, default=">",
                             help='''Relationship to value to filter on -- i.e. greater than, less than, etc. Accepted values are:
    gt, ge, lt, le, eq, ne -- respectively representing the relations >, >=, <, <=, ==, !=''')

    parser_filter.add_argument('-l','--latestage', type=str, required=True,
                             help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_filter.add_argument('-e','--earlystage', type=str, required=False, default=False,
                            help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')

    parser_filter_protocol = parser_filter.add_mutually_exclusive_group(required=True)

    parser_filter_protocol.add_argument('-s','--skipnorm', action='store_true', default=False,
                             help='''Use provided bedGraph (late option) directly -- skip any normalization procedure.''')
    parser_filter_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are median normalized.
    Then late stage is normalized to early stage if available..''')
    parser_filter_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
    Then they are median normalized.
    Then late stage is normalized to early stage if available.''')
    parser_filter_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first median normalized.
    Then they are smoothed with bandwidth given by --bandwidth.
    Then late stage is normalized to early stage if available.
    Note: if early is not present, this is same as protocol 4.''')
    parser_filter_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first median normalized.
    Then late stage is normalized to early stage if available.
    Then late/early is smoothed with bandwidth given by --bandwidth.
    Note: if early is not present, this is same as protocol 3.''')

    parser_filter_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                     help='''Late stage is normalized to early stage if available.
    Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth).
    Note: if early is not present, this is same as protocol 6.''')

    parser_filter_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
    Then late stage is normalized to early stage if available. (i.e. smooth -> L/E).
    Note: if early is not present, this is same as protocol 5.''')

    parser_filter.add_argument('-c', '--collapsed', action='store_true', default=False,
                            help='''Return collapsed variable-step bedGraph instead of expanded single-step bedGraph.
    This is often a much smaller file.''')

    parser_filter.add_argument('-ps', '--pseudo', type=float, default=0.1,
                            help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
    Should be between 0 and 1.
    Should be small enough to not change other values much,
    but big enough such that numbers divided by 0+pseudo do not become massive.
    Default: 0.1.''')

    parser_filter.add_argument('-bw', '--bandwidth', type=int, default=2500,
                        help=''' If kernel smoothing, specify bandwidth (int).
    Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
    Default: 2500.''')
    parser_filter.add_argument('--impute', type=int, default=False,
                            help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
    This bandwidth is generally longer than the one you would provide for regular smoothing.
    Only bins with a count of 0 will take on smoothed (imputed) values.
    Try: 10000.''')
     

    parser_filter.set_defaults(func=run_subtool)
 
 
 






    ## create sub-command for help
    parser_help = subparsers.add_parser('help', help=''' Gives more extensive guidance on using pufferfish.''')
    parser_help.set_defaults(func=run_subtool)


    
    ## parse the args and call the selected function
    args = parser.parse_args()


    


    ## check if args.quiet set (a default to all sub-modules from ArgumentParserWithDefaults class)
    if args.quiet:
        logger.setLevel(logging.ERROR)

    ## attempt to run args.func (which calls run_subtool() for all), catch errors
    try:
        args.func(parser, args)
    except IOError, e:
        ## often pipe will break for various reasons and will raise sigpipe error
        ## can import errno and do "if e.errno != errno.EPIPE:"
        ##   but errno.EPIPE = 32 -- so just use 32 here
        if e.errno != 32: ## ignore SIGPIPE
            raise


## Run main when this script used from command-line
if __name__ == "__main__":
    main()


## If needed parallelization
## Would have to parellize from here....
##    from joblib import Parallel, delayed
##    import time
##    from glob import glob
##    folder = "del"
##    files = glob('{}/*.txt'.format(folder))
##    def sleep(f):
##        print f
##        time.sleep(0.001)
##    ##    for f in files:
##    ##        sleep(f) #args.parallel
##    Parallel(n_jobs=2)(delayed(sleep)(f) for f in files)
