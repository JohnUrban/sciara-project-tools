#!/usr/bin/env python2.7
import sys, argparse, pickle, os, shutil
from collections import defaultdict
from itertools import islice, izip_longest
from Bio import SeqIO
import numpy.random as ran


parser = argparse.ArgumentParser(description="""
    Inputs
    1. modifications.csv from kinetics tools
    2. reference FASTA
    3. k = size
    4. cov position (can be 1-k).
    --> e.g. in a dimer use coverage from base at pos 1 (e.g. C in CG).
    --> e.g. in a trimer use coverage from base at pos 2 (e.g. C in GCG).

    modification.csv example lines:
    refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage,frac,fracLow,fracUp
    "contig_4",1983,1,T,10,3.450,1.572,0.677,5.100,3,,,
    "contig_4",1985,1,C,8,2.780,1.407,0.788,3.527,3,,,

    Output:
    Coverage levels over each kmer in format:
    contig, kmer, coverage, number of instances

    e.g.
    contig_1        AA      3       316

    Two possible ways:
    I. a. Read in Fasta and index every dimer in a dictionary {contig:{strand:{pos:kmer}}}
       b. Read in CSV and use the contig_pos to grab kmer there to assign cov in new dict: {contig:{kmer:{cov:num}}}

    II. a. Read in CSV making dictionary of {contig:{strand:{pos:cov}}}
        b. Read in Fasta, using dict1 to assign cov to kmers {contig:{kmer:{cov:num}}}
    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('-f', '--fasta', type=str, required=True,
                    help = '''Path to reference fasta.''')

parser.add_argument('-m', '--modcsv', type=str, required=True,
                    help = '''Path to modification.csv from kineticsTools.''')

parser.add_argument('-k', '--kmer', type=int, required=True,
                    help = '''Size of kmers. Priovide int. E.g. 2 for dimers, 3 for trimers.''')

parser.add_argument('-p', '--covpos', type=int, required=True,
                    help = '''0-based position in kmer to apply coverage to -- e.g. 0 for dimers (first base of 0,1) or 1 for trimers (mid base of 0,1,2).''')

parser.add_argument('-N', '--nlines', type=int, default=100000,
                    help = '''Read in this many lines from mod CSV to memory at a time... default 100000. Provide integer.''')

parser.add_argument('-D', '--debug_table', action='store_true', default=False,
                    help = '''Show a debug table to ensure you are counting kmers correctly.... This option is just for dev.''')

parser.add_argument('-Q', '--quiet', action='store_true', default=False,
                    help = '''Silence verbosity.''')


parser.add_argument('-L', '--lite', action='store_true', default=False,
                    help = '''Make this as memory-efficient as possible.''')
                    ## Can do by breaking up and pickling CSV data by sequence name
                    ## instead of pickling... it could also just process the KmerCov for each contig separately
                    ## except all seqs would need to be read in at once...
                    ## or an initial step would be to put indiv seq files in the tmpdir
args = parser.parse_args()



class ModCSV(object):
    def __init__(self, fh, nlines, convert_to_zero_base=True, debug=False, lite=False, tmpdir='tmpdir', runid=None, fasta=None, k=None, p=None, quiet=False, cleanupwhendone=False):
        self.fh = fh
        self.nlines = nlines
        self.dict = {}
        self.debug = debug
        self.lite = lite
        self.runid = str(runid) if runid is not None else str(ran.randint(1000000,9999999,1)[0])
        self.tmpdir = 'tmpdir-'+self.runid
        self.fasta = fasta
        self.quiet = quiet
        self.cleanupwhendone = cleanupwhendone
        self.k = k
        self.p = p
        self.reducePosBy = 1
        self.processed = []
        if not convert_to_zero_base:
            self.reducePosBy = 0
        self.currseq = None
        if self.lite:
            os.mkdir(self.tmpdir)
        if self.fasta is not None:
            self._parse_fasta()
        self._parse_file()
            

    def get(self, seqname, strand, pos):
        return self.dict[seqname][strand][pos]

    def clean(self):
        shutil.rmtree(self.tmpdir)

    def _reset_dict(self):
        self.dict = {}
        
    def _parse_fasta(self):
        for fa in SeqIO.parse(self.fasta, 'fasta'):
            with open("{}/{}.fasta".format(self.tmpdir,str(fa.name)), 'w') as fh:
                fh.write(">{}\n{}\n".format(str(fa.name), str(fa.seq)))

    def _parse_file(self):
        with open(self.fh) as f:
            self.header = f.next() ## burn this line off
            self._iterlines(f)
            # Ensuring that next was run on the last sequence in lite mode at end of file (EOF)
            self.line = dict(seqname="EOF")
            self._attemptKmerCov()
        if self.cleanupwhendone:
            shutil.rmtree(self.tmpdir)


    def _parse_line(self):
        self.line = self.line.strip().split(',')
        self.line = dict(seqname=self.line[0][1:-1],
                         pos=int(self.line[1]) - self.reducePosBy,
                         strand=int(self.line[2]),
                         cov=int(self.line[9]),
                         base=self.line[3])

    def _add_to_dict(self):
        try:
            self.dict[self.line['seqname']] 
        except KeyError:
            self.dict[self.line['seqname']] = {}
        try:
            self.dict[self.line['seqname']][self.line['strand']] 
        except KeyError:
            self.dict[self.line['seqname']][self.line['strand']] = {}
        # This should now be set up
        if self.debug:
            self.dict[self.line['seqname']][self.line['strand']][self.line['pos']] = (self.line['cov'], self.line['base'])
        else:
            self.dict[self.line['seqname']][self.line['strand']][self.line['pos']] = self.line['cov']
            
    def _iterlines(self, f):
        try:
            while f:
                #for line in islice(f, self.nlines):
                for next_n_lines in izip_longest(*[f] * self.nlines):
                    for self.line in next_n_lines:
                        self._parse_line()
                        self._attemptKmerCov()
                        self._add_to_dict()
                        self.currseq = self.line['seqname']
                                                  
                                                                 
        except:
            pass

    def _attemptKmerCov(self):
        ## Before adding to dict, see if current sequence can be processed
        if self.lite:
            if self.currseq is not None and self.line['seqname'] != self.currseq:
                ##print '{} {} {}'.format(self.line['seqname'], self.currseq, self.line['seqname']==self.currseq)
                #seqname has changed, process what was self.currseq
                assert self.currseq not in self.processed
                kmercov = KmerCov(fh='{}/{}.fasta'.format(self.tmpdir, self.currseq),
                                  covbypos=self,
                                  k=self.k,
                                  p=self.p,
                                  fastx='fasta')
                
                kmercov.print_table()
                sys.stderr.write('RunID {}: Finished {}....\n'.format(self.runid, self.currseq))
                self._reset_dict()
                self.processed.append(self.currseq)
        
                    


class KmerCov(object):
    def __init__(self, fh, covbypos, k=2, p=1, fastx='fasta'):
        '''
        covbypos:  ModCSV() object 
        '''
        #{contig:{kmer:{cov:num}}}
        self.fh = fh
        self.covbypos = covbypos #{contig:{strand:{pos:cov}}}
        self.debug = covbypos.debug
        self.k = k
        self.p = p
        self.fastx = fastx
        self.dict = {}
        self.ddict = {}
        self._parse_file()

    def print_table(self):
        if self.debug:
            self._debug_table()
        else:
            self._print_table()

    def _print_table(self):
        for seqname in sorted(self.dict.keys()):
            for kmer in sorted(self.dict[seqname].keys()):
                for cov in sorted(self.dict[seqname][kmer].keys()):
                    count = self.dict[seqname][kmer][cov]
                    print '\t'.join([str(e) for e in [seqname, kmer, cov, count]])

    def _debug_table(self):
        for seqname in sorted(self.ddict.keys()):
            for kmer in sorted(self.ddict[seqname].keys()):
                for cov in sorted(self.ddict[seqname][kmer].keys()):
                    count = self.ddict[seqname][kmer][cov]
                    print '\t'.join([str(e) for e in [seqname, kmer, cov, count]])


    def _parse_seq(self, seq, rev=False):
        correction = 0
        if rev:
            correction = self.falen - 1
        for self.pos in range(self.falen - self.k + 1):
            self.kmer = seq[self.pos : self.pos + self.k]
            # This needs to be self.pos + self.p
            try:
                pos = abs(correction - (self.pos + self.p))
                if self.debug:
                    self.cov, self.base = self.covbypos.get(self.seqname, self.strand, pos)
                else:
                    self.cov =  self.covbypos.get(self.seqname, self.strand, pos)
                    self.base = '.'
                #print self.kmer, self.base, self.p, self.pos, pos
            except: ## This was not a position found in ModCSV()
                self.cov = 0
                self.base = "."
                
        
            self._add_to_dict()
            self._debug_dict()

    def _parse_fwd(self):
        self._parse_seq(str(self.fa.seq))

    def _parse_rev(self):
        self._parse_seq(str(self.fa.reverse_complement().seq), rev=True)

    def _add_to_dict(self):
        #{contig:{kmer:{cov:num}}}
        try:
            self.dict[self.seqname] 
        except KeyError:
            self.dict[self.seqname] = {}
        try:
            self.dict[self.seqname][self.kmer] 
        except KeyError:
            self.dict[self.seqname][self.kmer] = {}
        try:
            self.dict[self.seqname][self.kmer][self.cov] += 1
        except KeyError:
            self.dict[self.seqname][self.kmer][self.cov] = 1
            
    def _debug_dict(self):
        #{contig:{kmer:{cov:num}}}
        try:
            self.ddict[self.seqname] 
        except KeyError:
            self.ddict[self.seqname] = {}
        try:
            self.ddict[self.seqname][self.kmer] 
        except KeyError:
            self.ddict[self.seqname][self.kmer] = {}
        try:
            self.ddict[self.seqname][self.kmer][self.base] += 1
        except KeyError:
            self.ddict[self.seqname][self.kmer][self.base] = 1
            
    def _parse_file(self):
        for self.fa in SeqIO.parse(self.fh, self.fastx):
            self.falen = len(self.fa)
            self.seqname = str(self.fa.name)
            self.strand = 0
            self._parse_fwd()
            self.strand = 1
            self._parse_rev()


def run(args):
    runid = str(ran.randint(1000000,9999999,1)[0])

    if not args.lite:
        if not args.quiet:
            sys.stderr.write('RunID {}: ModCSV step....\n'.format(runid))
        modcsv = ModCSV(fh=args.modcsv, nlines=args.nlines, debug=args.debug_table, lite=False, runid=runid)

        if not args.quiet:
            sys.stderr.write('RunID {}: KmerCov step....\n'.format(runid))
        kmercov = KmerCov(fh=args.fasta, covbypos=modcsv, k=args.kmer, p=args.covpos, fastx='fasta')

        if not args.quiet:
            sys.stderr.write('RunID {}: Print step....\n'.format(runid))
        kmercov.print_table()
    else:
        sys.stderr.write('RunID {}: Combined ModCSV+KmerCov steps....\n'.format(runid))
        ModCSV(args.modcsv, args.nlines, debug=args.debug_table, lite=True,
               runid=runid, fasta=args.fasta,
               k=args.kmer, p=args.covpos,cleanupwhendone=True)






try:
    run(args)
except IOError:
    pass

