#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict

# -l 4  -m 33  -h 75

parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--assignments', '-a',
                   type= str, default=False, required=True,
                   help='''Path to assignments table output by: "convert-step5-allasm-to-table.py -f all.sizes --sort_by_size > table-sort.txt" .''')

parser.add_argument('--purge_flags', '-p',
                   type= str, default=False, 
                   help='''Path to purge_haplotigs flagging CSV.''')

parser.add_argument('--gaps', '-g',
                   type= str, default=False, 
                   help='''Path to BED file of N gaps in scaffolds. This is used to re-adjust the number of bases with 0 counts from BEDtools output (and/or purge).''')

parser.add_argument('--coverage', '-c',
                   type= str, default=False, 
                   help='''Path to "bedtools genomecov" coverage histogram output.''')

parser.add_argument('--lowlow', '-:',
                   type= int, default=1,
                   help='''Low cutoff used/desired.''')

parser.add_argument('--low', '-l',
                   type= int, default=4,
                   help='''Low cutoff used/desired.''')

parser.add_argument('--med', '-m',
                   type= int, default=33,
                   help='''Middle cutoff used/desired.''')

parser.add_argument('--high', '-H',
                   type= int, default=75,
                   help='''High cutoff used/desired.''')

parser.add_argument('--rewrite_purge', '-RP', type=str, default=False,
                    help=''' Rewrite purge CSV with regular contig names found in assignments file. Give outputfilename.''')

parser.add_argument('--rewrite_coverage', '-RC', type=str, default=False,
                    help=''' Rewrite coverage file. This will have gap adjustments as well if gapfile provided. Else "gap adjustments" will equal the unadjusted.''')

parser.add_argument('--callX_from_purge', '-CXP', type=float, default=False, help='''Output likely X names. Provide cutoff for what percent of contig bases should be considered haploid by purge_haplotigs - e.g. 51 or 80.
Also see --include_rmtigs, --include_redundant, --include_all_purge_class.
By default rmtigs, redundant tigs, and anything but KEEP, is ignored.''')

parser.add_argument('--include_rmtigs', action='store_true', default=False)
parser.add_argument('--include_redundant', action='store_true', default=False)
parser.add_argument('--include_all_purge_class', action='store_true', default=False)

parser.add_argument('--classify', type=str, default=False, help=''' Returns updated assignment table. Provide filename.''')
args = parser.parse_args()





class Genome(object):
    def __init__(self, assignments, flags=False, gaps=False, coverage=False):
        self.afh = assignments
        self.ffh = flags
        self.knownheader = ["contig","length","tax_label","is_remove_tig","contains_rmtig", "is_redundant","purge_haplotigs_class","info_name"]
        self.knownpurge = ["contig","contig_reassign","bases_hap_dip","bases_low_high","bases_all","perc_low_coverage","perc_hap_coverage","perc_dip_coverage","perc_high_coverage"]
        self.header = {}
        for i in range(len(self.knownheader)):
            self.header[self.knownheader[i]] = i-1 ## b/c not using element 0 later
        self.order = []
        self.genome = {}
        self.purge = {}
        self.gapsums = defaultdict(int)
        self.coverage = defaultdict(dict)
        
        self.add_assignments(assignments)
        if flags:
            self.add_purge_flags(flags)
        if gaps:
            self.add_gaps(gaps)
        if coverage:
            self.add_coverage(coverage)
            
        
    def add_assignments(self, fh):
        with open(fh) as f:
            header = f.next().strip().lstrip('#').split()
            assert header == self.knownheader
            for line in f:
                line = line.strip().split()
                self.genome[line[0]] = line[1:]
                self.order.append( line[0] )
                
        self.genomelength = sum( [int(self.genome[e][0]) for e in self.genome.keys()] )

    def _gettigname_from_purge(self, purgename):
        splitname = purgename.split('_')
        if len(splitname) < 3:
            ## In case this is used on non purge names
            return '_'.join(splitname)
        else:
            return '_'.join(splitname[1:3])

    def add_purge_flags(self, fh):
        with open(fh) as f:
            header = f.next().strip().lstrip('#').split(',')
            assert header == self.knownpurge
            pkeys=[]
            for line in f:
                line = line.strip().split(',')
                tigname = self._gettigname_from_purge(line[0])
                self.purge[tigname] = line[1:]
                pkeys.append(tigname)
            for key in self.genome.keys():
                if key not in pkeys:
                    self.purge[key] = ['S'] + ['0']*8

    def get_haploid_contigs_from_purge(self, perc_hap_cov_cut=51, reject_rmtigs=True, reject_redundant=True, include_only_keep=True):
        #print perc_hap_cov_cut, reject_rmtigs, reject_redundant, include_only_keep
        for contig in self.purge.keys():
            if float(self.purge[contig][5]) >= perc_hap_cov_cut:
                if reject_rmtigs and self.isRmTig(contig):
                    continue
                if reject_redundant and self.isRedundant(contig):
                    continue
                if include_only_keep and not self.purge_class_is_keep(contig):
                    continue
                print contig

    def add_gaps(self, fh):
        with open(fh) as f:
            keys=set([])
            for line in f:
                line = line.strip().split()
                self.gapsums[line[0]] += int(line[2]) - int(line[1])
                keys.add(line[0])
                
            for key in self.genome.keys():
                if key not in keys:
                    self.gapsums[key] = 0

    def _add_cov_line(self, line):
        line = line.strip().split()
        contig = self._gettigname_from_purge(line[0])
        depth = int(line[1])
        nbases = int(line[2])
        if len(line) > 3:
            length = int(line[3])
        else:
            if contig == 'genome':
                 length = self.genomelength 
            else:
                length = int(self.genome[contig][0])

        if contig == 'genome':
            assert length == self.genomelength 
        else:
            assert length == int(self.genome[contig][0])
        pct = nbases / float(length)
        # adjust length if gap info available
        seqlength = length - self.gapsums[contig]
        if depth == 0:
            nbases2 = nbases - self.gapsums[contig]
        else:
            nbases2 = nbases
        seqpct = nbases2 / float(seqlength)
        try:
            self.coverage[contig][depth] = [nbases, length, pct, self.gapsums[contig], seqlength, seqpct]
        except:
            sys.stderr.write("WARN: See _add_cov_line\n")
            quit()
        return contig

    def _get_cov_lines(self, contig):
        lines = ''
        for depth in sorted(self.coverage[contig]):
            line = '\t'.join([str(e) for e in [contig, depth]+ self.coverage[contig][depth]])
            lines += line + '\n'
        return lines

    def add_coverage(self, fh):
        ## computing gaps first will allow this to adjust for gaps
        with open(fh) as f:
            keys=set([])
            for line in f:
                key = self._add_cov_line(line)
                keys.add(key)
                
            for key in self.genome.keys():
                if key not in keys:
                    line = '\t'.join([key, '0', self.genome[key][0], self.genome[key][0], '100'])
                    self._add_cov_line(line)

    def isRmTig(self, contig):
        return self.genome[contig][self.header['is_remove_tig']] == 'yes'

    def isRedundant(self, contig):
        return self.genome[contig][self.header['is_redundant']] == 'yes'

    def isBacterial(self, contig):
        return self.genome[contig][self.header['tax_label']] == 'bacterial'



    def purge_class_is_keep(self, contig):
        return self.genome[contig][self.header['purge_haplotigs_class']] == 'KEEP'

    def rewrite_purge(self, fh):
        out = open(fh, 'w')
        out.write('#'+','.join(self.knownpurge) + '\n')
        for contig in self.purge.keys():
            line = ','.join([contig]+self.purge[contig])
            out.write(line + '\n')


    def rewrite_coverage(self, fh):
        out = open(fh, 'w')
        for contig in sorted(self.coverage.keys()):
            if contig != 'genome':
                out.write(self._get_cov_lines(contig))
        out.write(self._get_cov_lines('genome'))

    def _get_classify_line(self, contig, lowlow, low, mid, high):
        perc = {'zero':0, 'junklow':0, 'low':0,'haplo':0,'diplo':0,'high':0}
        for depth in sorted(self.coverage[contig].keys()):
            if depth == 0:
                perc["zero"] += self.coverage[contig][depth][5]
            elif depth > 0 and depth < lowlow:
                perc["junklow"] += self.coverage[contig][depth][5]
            elif depth >= lowlow and depth < low:
                perc['low'] += self.coverage[contig][depth][5]
            elif depth >= low and depth < mid:
                perc['haplo'] += self.coverage[contig][depth][5]
            elif depth >= mid and depth <= high:
                perc['diplo'] += self.coverage[contig][depth][5]
            elif depth > high:
                perc['high'] += self.coverage[contig][depth][5]
            else:
                sys.stderr.write('WARN: See _get_classify_lines....\n')
                quit()
        cov_class = 'None'
        cov_class_perc = 0
        covclasses=['zero', 'junklow', 'low','haplo','diplo','high']
        #print perc
        for covclass in covclasses:
            if perc[covclass] > cov_class_perc:
                #print perc[covclass] , cov_class_perc
                cov_class = covclass
                cov_class_perc = perc[covclass]
        my_purge_class = 'None'
        chr_class = 'None'
        if self.isRmTig(contig):
            my_purge_class = 'JUNK'
            chr_class = 'RMTIG'
        elif self.isRedundant(contig) and cov_class == 'haplo':
            my_purge_class = 'HAPLO'
            chr_class = 'REDHAP'
        elif self.isRedundant(contig) and cov_class == 'high':
            my_purge_class = 'REPEAT'
            chr_class = 'R'
        elif self.isRedundant(contig) and cov_class == 'low':
            my_purge_class = 'JUNK'
            chr_class = 'REDLOW'
        elif self.isRedundant(contig) and cov_class == 'lowlow':
            my_purge_class = 'JUNK_REDLOWLOW'
            chr_class = 'REDLOWLOW'
        elif self.isRedundant(contig) and cov_class == 'zero':
            my_purge_class = 'JUNK_REDZERO'
            chr_class = 'REDZERO'
        elif self.isRedundant(contig):
            my_purge_class = 'JUNK_RED' ## should not be able to get here (on purpose)
            chr_class = 'RED'
        elif cov_class == 'low':
            my_purge_class = 'JUNK_LOW'
            chr_class = 'L'
        elif cov_class == 'junklow':
            my_purge_class = 'JUNK_LOWLOW'
            chr_class = 'L'
        elif cov_class == 'zero':
            my_purge_class = 'JUNK_ZERO'
            chr_class = 'L'
        elif cov_class == 'high':
            my_purge_class = 'JUNK_HIGH'
            chr_class = 'R'
        elif cov_class == 'diplo':
            my_purge_class = 'KEEP_DIPLO'
            chr_class = 'A'
        elif cov_class == 'haplo':
            my_purge_class = 'KEEP_HAPLO'
            chr_class = 'X'
        else:
            my_purge_class = 'KEEP_CHECK'
            chr_class = 'None'

        if self.isBacterial(contig):
            chr_class = 'B'

        #print contig, cov_class, chr_class, my_purge_class
        out = [contig] + self.genome[contig] + [perc[e] for e in covclasses] + [cov_class, chr_class, my_purge_class]
        return '\t'.join( [str(e) for e in out] )

        

    def write_classification(self, fh, lowlow=1, low=4, mid=33, high=75):
        self.headerext = ['perc_zero','perc_junklow','perc_low', 'perc_haplo', 'perc_diplo', 'perc_high', 'cov_class', 'chr_class', 'my_purge_class']
        out = open(fh, 'w')
        out.write('#'+','.join(self.knownheader + self.headerext) + '\n')
        for contig in self.order:
            out.write(self._get_classify_line(contig, lowlow, low, mid, high)+'\n')

    


### EXECUTE
genome = Genome(assignments=args.assignments, flags=args.purge_flags, gaps=args.gaps, coverage=args.coverage)

if args.rewrite_purge:
    assert args.purge_flags
    genome.rewrite_purge(args.rewrite_purge)

if args.callX_from_purge:
    assert args.purge_flags
    genome.get_haploid_contigs_from_purge(perc_hap_cov_cut=args.callX_from_purge, reject_rmtigs= not args.include_rmtigs, reject_redundant=not args.include_redundant, include_only_keep=not args.include_all_purge_class)

if args.rewrite_coverage:
    assert args.coverage
    genome.rewrite_coverage(args.rewrite_coverage)

if args.classify:
    assert args.classify
    assert args.coverage
    genome.write_classification(args.classify)

### TESTS
##for key in genome.genome.keys():
##    print genome.genome[key]
##    print "rmtig", "redundant", "keep"
##    print genome.isRmTig(key), genome.isRedundant(key), genome.purge_class_is_keep(key)
##    print
