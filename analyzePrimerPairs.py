#!/usr/bin/env python
import sys, os, subprocess
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Copy/paste Primer Blast results into text file.
    This goes through text file and separates primers into fasta file pairs.

    TODO:
    - go from fasta file of primers (file pairs and/or interleaved) -- also option to add pair number to name
    - go from tab-sep file (with col specified for name,rev,fwd) -- also option to add pair number to name
    - rank option (by least sum)
    - be able to use multiple bt2 and multiple blastdbs (i.e. from more than 1 assrmbly)
    - be able to do just primer or just product analyses
    - be able to specify entire pipeline --
    - some verbosity levels

    -TODO: product only still writes primer files... put a stop to that...
    
    """, formatter_class= argparse.RawTextHelpFormatter)

##parser_input = parser.add_mutually_exclusive_group()
parser.add_argument('-i', "--input",
                   type=str, required=False,
                   help='''Path to primer blast file.''')

parser.add_argument('-I', "--interleaved",
                   type=str, required=False,
                   help='''Path to interleaved fwd/rev primers fasta file.''')

parser.add_argument('-f', "--fwd",
                   type=str, required=False,
                   help='''Path to fwd primer fasta file.''')

parser.add_argument('-r', "--rev",
                   type=str, required=False,
                   help='''Path to reverse primer fasta file.''')

parser.add_argument('-t', "--tsv",
                   type=str, required=False, default=False,
                   help='''Path to reverse primer fasta file.''')

parser.add_argument('-fc', "--fwdcol",
                   type=int, required=False, default=False,
                   help='''1-based column number for forward primer when using --tsv''')

parser.add_argument('-rc', "--revcol",
                   type=int, required=False, default=False,
                   help='''1-based column number for reverse primer when using --tsv''')

parser.add_argument('-nc', "--namecol",
                   type=int, required=False, default=False,
                   help='''1-based column number for name column when using --tsv. Not using this will just result in primer pair named by number -- where relevant.''')

parser.add_argument('-o', "--outdir",
                   type=str, required=True,
                   help='''Path to output dir - where files will be deposited.''')

parser.add_argument('-bt2', "--bowtie2index",
                   type=str, required=False,
                   help='''Path to bowtie2 index (up to prefix of index files).''')

parser.add_argument('-bdb', "--blastdb",
                   type=str, default=None,
                   help='''Path to blast db (up to prefix of db files).''')

parser.add_argument('-ref', "--reference",
                   type=str, default=None,
                   help='''Path to reference.fasta file.''')

parser.add_argument('-p', "--product_only",
                   action="store_true", default=False,
                   help='''Specifies that the only thing you want is to obtain the PCR target product(s)
of a primer pair or set of primer pairs (in any of the the input file types).
Example use:
analyzePrimerPairs.py -I interleaved.fa -o ./ -bt2 bt2index-prefix -ref ref.fa''')


parser.add_argument('-P', "--product_spec_only",
                   type=str, default=False,
                   help='''TODO
Specifies that you are providing products and the only thing you want
is to check its specificity.
Example use:
analyzePrimerPairs.py -P products.fa -o ./ -bdb blast-index-prefix -ref ref.fa''')


parser.add_argument('-v', "--verbose",
                   action="store_true", default=False,
                   help='''Say stuff while doing stuff.''')

args = parser.parse_args()



if bool(args.tsv):
    if not (bool(args.fwdcol) and bool(args.revcol)):
        sys.stderr.write('\nInput argument error: --tsv, --fwdcol, --revcol must be given together\n\n')
        quit()
if bool(args.fwd):
    if not bool(args.rev):
        sys.stderr.write('\nInput argument error: --fwd, --rev must be given together\n\n')
        quit()


class Result(object):
    def __init__(self,result):
        self.result = result
        self.i = 0
    def next(self):
        self.i += 1
        return self.result[self.i-1]
    def reset(self):
        self.i = 0
    def __str__(self):
        return str(self.result)
        

class Pairs(object):
    def __init__(self,f,outdir,bt2,bdb=None,reference=None):
        '''f is Path to primer blast file.'''
        self.f = f
        if outdir[-1] != "/":
            outdir += "/"
        self.outdir = outdir
        self.pairs = {}
        self.i = 1
        self.info_paths = {}
        self.primer_paths = {}
        self.bt2 = bt2
        self.blastdb = bdb
        self.reference = reference
        self.primerSpecOut_paths = {}
        self.primerSpecResults = {}
        self.productSpecOut_paths = {}
        self.productSpecResults = {}
        self.results = None
        self.primer_test = False
        self.product_test = False
        
##    def readPrimerPair(self):#primer_blast
##        pair = {}
##        pair["name"] = self.fh.next()
##        pair["labels"] = self.fh.next()
##        pair["fwd"] = self.fh.next()
##        pair["rev"] = self.fh.next()
##        pair["length"] = self.fh.next()
##        return pair

    def readPrimerPairs(self):#could take any function "readPrimerPair"
        self.i = 1
        self.fh = open(self.f, 'r')
        while self.fh:
            try:
                pair = self.readPrimerPair()
                self.pairs[self.i] = pair
                self.i += 1
            except:
                break
        self.fh.close()

    def get_pair_nums(self):
        return sorted(self.pairs.keys())

##    def get_fwd(self,pair):#primer_blast
##        return self.pairs[pair]["fwd"].split()[2]
##
##    def get_rev(self,pair):#primer_blast
##        return self.pairs[pair]["rev"].split()[2]

    def get_primer_path(self,pair):
        if not self.primer_paths:
            self.write_primer_files()
        return self.primer_paths[pair]

##    def get_pair_info_path(self, pair):#primer_blast
##        if not self.info_paths:
##            self.write_pair_stats()
##        return self.info_paths[pair]

    def get_primerSpecResult(self,pair):
        return self.primerSpecResults[pair]
    
##    def write_pair_stats(self):#primer_blast
##        for pair in self.get_pair_nums():
##            path = self.outdir+"primer-pair-"+str(pair)+"-info.txt"
##            self.info_paths[pair] = path
##            o = open(path, 'w')
##            o.write(self.pairs[pair]["name"]+"\n")
##            o.write(self.pairs[pair]["labels"]+"\n")
##            o.write(self.pairs[pair]["fwd"]+"\n")
##            o.write(self.pairs[pair]["rev"]+"\n")
##            o.write(self.pairs[pair]["length"]+"\n")
##            o.write("\n")
##            o.close()

    def write_primer_files(self):
        for pair in self.get_pair_nums():
            path_fwd = self.outdir+"primer-1-"+str(pair)+".fa"
            path_rev = self.outdir+"primer-2-"+str(pair)+".fa"
            self.primer_paths[pair] = [path_fwd,path_rev]
            o = open(path_fwd, 'w')
            o.write(">p1\n")
            o.write(self.get_fwd(pair)+"\n")
            o.close()
            o = open(path_rev, 'w')
            o.write(">p2\n")
            o.write(self.get_rev(pair)+"\n")
            o.close()

    def primerSpecificity(self, verbose=False):
        for pair in pairs.get_pair_nums():
            if verbose:
                sys.stderr.write("\t Pair " + str(pair) + "\n")
            cmd = self._ps_cmd(pair)
            result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            self.primerSpecResults[pair] = result
        self.primer_test = True


    def productSpecificity(self, verbose=False):
        if self.blastdb is not None:
            for pair in pairs.get_pair_nums():
                if verbose:
                    sys.stderr.write("\t Pair " + str(pair) + "\n")
                cmd = self._target_cmd(pair)
                result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                self.productSpecResults[pair] = result
            self.product_test = True
            
    def writeProduct(self, verbose=False):
        if self.reference is not None:
            for pair in pairs.get_pair_nums():
                if verbose:
                    sys.stderr.write("\t Pair " + str(pair) + "\n")
                cmd = self._only_get_target_sequence_cmd(pair)
                result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                if verbose:
                    sys.stderr.write(result+"\n")
                    
                

    def write_primerSpecResults(self):
        if self.primer_test:
            for pair in self.get_pair_nums():
                path = self.outdir+"primer-pair-"+str(pair)+".specificity.txt"
                self.primerSpecOut_paths[pair] = path
                o = open(path, 'w')
                o.write(self.primerSpecResults[pair]+"\n")
                o.close()

    def write_productSpecResults(self):
        if self.product_test:
            for pair in self.get_pair_nums():
                path = self.outdir+"primer-pair-"+str(pair)+".product.txt"
                self.productSpecOut_paths[pair] = path
                o = open(path, 'w')
                o.write(self.productSpecResults[pair]+"\n")
                o.close()

    def parse_primerSpecResult(self, pair):
        if self.primer_test:
            result = self.primerSpecResults[pair].split("\n")
            i = 0
            while result[i] != "FINAL REPORT:":
                i+=1
            parsed_result = []
            for j in range(i,len(result)):
                parsed_result.append(result[j])
            parsed_result = Result(parsed_result)
            return parsed_result

    def parse_productSpecResult(self, pair):
        if self.product_test:
            result = self.productSpecResults[pair].split("\n")
            i = 0
            while result[i] != "PCR TARGET BLAST RESULTS":
                i+=1
            parsed_result = []
            for j in range(i,len(result)):
                parsed_result.append(result[j])
            parsed_result = Result(parsed_result)
            return parsed_result


    def analyze_primerSpecResult(self,pair):
        if self.primer_test:
            f = self.parse_primerSpecResult(pair)
            counts = []
            f.next() # burn off "FINAL REPORT:"
            f.next() #
            counts.append( f.next().split()[4] )
            counts.append( f.next().split()[4] )
            counts.append( f.next().split()[4] )
            counts.append( f.next().split()[4] )
            f.next()
            f.next()
            counts.append( f.next().split()[2] )
            counts.append( f.next().split()[2] )
            f.next()
            f.next()
            counts.append( f.next().split()[2] )
            counts.append( f.next().split()[2] )
            return counts

    def analyze_productSpecResult(self,pair):
        if self.product_test:
            f = self.parse_productSpecResult(pair)
            counts = []
            f.next() # burn off "PCR TARGET BLAST RESULTS"
            counts.append( f.next().split()[5] )
            return counts

    def analyze_primerSpecResults(self):
        if self.results is None:
            self._initialize_results()
        if self.primer_test:
            for pair in self.get_pair_nums():
                counts = self.analyze_primerSpecResult(pair)
                self.results[pair]['primers'] = counts

    def analyze_productSpecResults(self):
        if self.results is None:
            self._initialize_results()
        if self.product_test:
            for pair in self.get_pair_nums():
                counts = self.analyze_productSpecResult(pair)
                self.results[pair]['target'] = counts

    def analyze_results(self):
        if not self.primer_test:
            self.primerSpecificity()
        if not self.product_test:
            self.productSpecificity
        self.analyze_primerSpecResults()
        self.analyze_productSpecResults()

    def print_results(self):
        for pair in self.get_pair_nums():
            results = [str(pair), self.get_fwd(pair), self.get_rev(pair)] + [e for e in self.results[pair]['primers'] + self.results[pair]['target']]
            print ("\t").join(results)

    def _initialize_results(self):
        self.results = {}
        for pair in self.get_pair_nums():
            self.results[pair] = {'primers':[], 'target':[]}
            
    def get_fwd(self,pair):#all but primer blast (needs own)
        return self.pairs[pair]["fwd"]

    def get_rev(self,pair):#all but primer blast (needs own)
        return self.pairs[pair]["rev"]        

    def _ps_cmd(self, pair):
        return ["primerSpecificity.sh",  self.bt2, pairs.get_primer_path(pair)[0], pairs.get_primer_path(pair)[1]]

    def _target_cmd(self, pair):
        return ["productSpecificity.sh",  self.bt2, pairs.get_primer_path(pair)[0], pairs.get_primer_path(pair)[1], self.blastdb, self.reference]

    def _only_get_target_sequence_cmd(self, pair):
        return ["get_PCR_product_from_primers.sh",  self.bt2, pairs.get_primer_path(pair)[0], pairs.get_primer_path(pair)[1], self.reference]



class PrimerBlastFile(Pairs):
    def __init__(self,f,outdir,bt2,bdb=None,reference=None):
        '''f is Path to primer blast file.'''
        self.f = f
        if outdir[-1] != "/":
            outdir += "/"
        self.outdir = outdir
        self.pairs = {}
        self.i = 1
        self.info_paths = {}
        self.primer_paths = {}
        self.bt2 = bt2
        self.blastdb = bdb
        self.reference = reference
        self.primerSpecOut_paths = {}
        self.primerSpecResults = {}
        self.productSpecOut_paths = {}
        self.productSpecResults = {}
        self.results = None
        self.primer_test = False
        self.product_test = False
        
    def readPrimerPair(self):#primer_blast
        pair = {}
        pair["name"] = self.fh.next()
        pair["labels"] = self.fh.next()
        pair["fwd"] = self.fh.next()
        pair["rev"] = self.fh.next()
        pair["length"] = self.fh.next()
        return pair

    def get_fwd(self,pair):#primer_blast
        return self.pairs[pair]["fwd"].split()[2]

    def get_rev(self,pair):#primer_blast
        return self.pairs[pair]["rev"].split()[2]


    def get_pair_info_path(self, pair):#primer_blast
        if not self.info_paths:
            self.write_pair_stats()
        return self.info_paths[pair]

    
    def write_pair_stats(self):#primer_blast
        for pair in self.get_pair_nums():
            path = self.outdir+"primer-pair-"+str(pair)+"-info.txt"
            self.info_paths[pair] = path
            o = open(path, 'w')
            o.write(self.pairs[pair]["name"]+"\n")
            o.write(self.pairs[pair]["labels"]+"\n")
            o.write(self.pairs[pair]["fwd"]+"\n")
            o.write(self.pairs[pair]["rev"]+"\n")
            o.write(self.pairs[pair]["length"]+"\n")
            o.write("\n")
            o.close()

    

class PairedFastaFiles(Pairs):
    def __init__(self,f1,f2,outdir,bt2,bdb=None,reference=None):
        '''f is Path to primer blast file.'''
        self.f1 = f1
        self.f2 = f2
        if outdir[-1] != "/":
            outdir += "/"
        self.outdir = outdir
        self.pairs = {}
        self.i = 1
        self.info_paths = {}
        self.primer_paths = {}
        self.bt2 = bt2
        self.blastdb = bdb
        self.reference = reference
        self.primerSpecOut_paths = {}
        self.primerSpecResults = {}
        self.productSpecOut_paths = {}
        self.productSpecResults = {}
        self.results = None
        self.primer_test = False
        self.product_test = False
        
    def readPrimerPair(self, j):#paired_fasta
        pair = {}
        pair["name"] = self.fwd[j].name
##        pair["labels"] = self.fh.next()
        pair["fwd"] = str(self.fwd[j].seq)
        pair["rev"] = str(self.rev[j].seq)
##        pair["length"] = self.fh.next()
        return pair

    def readPrimerPairs(self):#could take any function "readPrimerPair"
        self.i = 1
        self.fwd = [e for e in SeqIO.parse(self.f1, "fasta")]
        self.rev = [e for e in SeqIO.parse(self.f2, "fasta")]
        for j in range(len(self.fwd)):
            try:
                pair = self.readPrimerPair(j)
                self.pairs[self.i] = pair
                self.i += 1
            except:
                break


class InterleavedFastaFile(Pairs):
    def __init__(self,f,outdir,bt2,bdb=None,reference=None):
        '''f is Path to primer blast file.'''
        self.f = f
        if outdir[-1] != "/":
            outdir += "/"
        self.outdir = outdir
        self.pairs = {}
        self.i = 1
        self.info_paths = {}
        self.primer_paths = {}
        self.bt2 = bt2
        self.blastdb = bdb
        self.reference = reference
        self.primerSpecOut_paths = {}
        self.primerSpecResults = {}
        self.productSpecOut_paths = {}
        self.productSpecResults = {}
        self.results = None
        self.primer_test = False
        self.product_test = False
        
    def readPrimerPair(self, j):#paired_fasta
        pair = {}
        pair["name"] = self.fasta[j].name
##        pair["labels"] = self.fh.next()
        pair["fwd"] = str(self.fasta[j].seq)
        pair["rev"] = str(self.fasta[j+1].seq)
##        pair["length"] = self.fh.next()
        return pair

    def readPrimerPairs(self):#could take any function "readPrimerPair"
        self.i = 1
        self.fasta = [e for e in SeqIO.parse(self.f, "fasta")]
        for j in range(0, len(self.fasta), 2):
            try:
                pair = self.readPrimerPair(j)
                self.pairs[self.i] = pair
                self.i += 1
            except:
                break
          

class TSV(Pairs):
    def __init__(self,f,fwd_col, rev_col, outdir,bt2,bdb=None,reference=None,name_col=False):
        '''f is Path to primer blast file.'''
        self.f = f
        if outdir[-1] != "/":
            outdir += "/"
        self.outdir = outdir
        self.pairs = {}
        self.i = 1
        self.info_paths = {}
        self.primer_paths = {}
        self.bt2 = bt2
        self.blastdb = bdb
        self.reference = reference
        self.primerSpecOut_paths = {}
        self.primerSpecResults = {}
        self.productSpecOut_paths = {}
        self.productSpecResults = {}
        self.results = None
        self.primer_test = False
        self.product_test = False
        if name_col:
            self.name_col = name_col-1 #pythonese
        else:
            self.name_col = name_col
        self.fwd_col = fwd_col-1 #pythonese
        self.rev_col = rev_col-1 #pythonese
        
    def readPrimerPair(self, j):#paired_fasta
        pair = {}
        if self.name_col:
            pair["name"] = self.tsv[j].strip().split()[self.name_col]
        else:
            pair['name'] = "primer-pair-"+str(j)
##        pair["labels"] = self.fh.next()
        pair["fwd"] = self.tsv[j].strip().split()[self.fwd_col]
        pair["rev"] = self.tsv[j].strip().split()[self.rev_col]
##        pair["length"] = self.fh.next()
        return pair

    def readPrimerPairs(self):#could take any function "readPrimerPair"
        self.i = 1
        self.fh = open(self.f, 'r')
        self.tsv = self.fh.readlines()
        for j in range(len(self.tsv)):
            try:
                pair = self.readPrimerPair(j)
                self.pairs[self.i] = pair
                self.i += 1
            except:
                break



    
##pairs = Pairs(args.input, args.outdir, args.bowtie2index, args.blastdb, args.reference)

if args.input:
    pairs = PrimerBlastFile(args.input, args.outdir, args.bowtie2index, args.blastdb, args.reference)
elif args.interleaved:
    pairs = InterleavedFastaFile(args.interleaved, args.outdir, args.bowtie2index, args.blastdb, args.reference)
elif args.fwd and args.rev:
    pairs = PairedFastaFiles(args.fwd, args.rev, args.outdir, args.bowtie2index, args.blastdb, args.reference)
elif args.tsv and args.fwdcol and args.revcol:
    pairs = TSV(args.tsv, args.fwdcol, args.revcol, args.outdir, args.bowtie2index, args.blastdb, args.reference, args.namecol)
elif args.product_spec_only:
    pass


if args.verbose:
    sys.stderr.write("Reading primer pairs in...\n")
pairs.readPrimerPairs()

if args.product_only:
    if args.verbose:
        sys.stderr.write("Writing PCR target prduct file ...\n")
    pairs.writeProduct(verbose=args.verbose)
    quit()

if args.input:
    if args.verbose:
        sys.stderr.write("Writing out primer pair information files in...\n")
    pairs.write_pair_stats()

if args.verbose:
    sys.stderr.write("Writing primer fastas...\n")
pairs.write_primer_files()

if args.verbose:
    sys.stderr.write("Checking primer specificity...\n")
pairs.primerSpecificity(verbose=args.verbose)

if args.verbose:
    sys.stderr.write("Writing primer specificity results...\n")
pairs.write_primerSpecResults()

if args.verbose:
    sys.stderr.write("Checking product specificity...\n")
pairs.productSpecificity(verbose=args.verbose)

if args.verbose:
    sys.stderr.write("Writing product specificity results...\n")
pairs.write_productSpecResults()

if args.verbose:
    sys.stderr.write("Analyzing results...\n")
pairs.analyze_results()


pairs.print_results()





################################################################################################
##def runfail(cmd):
##    sys.stderr.write( "Running : %s \n" % cmd )
##    if not 0 == os.system(cmd):
##        sys.exit("Failed : %s " % cmd)           
##
##def bowtie2_from_files(bt2, fwd, rev):
##    pass
##
##def blastn_from_files(bdb, product):
##    pass
##
##def run_primerSpecificity_script(bt2, fwd, rev, bdb, ref):
##    pass

## QC on arguments ##
##start_path = os.getcwd()
##if not os.path.exists(args.outdir):
##    sys.exit("Outdir does not exist.")
##if not os.path.isabs(query_file):
##    query_file = os.path.join(start_path, query_file)
##if not no_ref and not os.path.isabs(ref_file):
##    ref_file = os.path.join(start_path, ref_file)
##for pair in pairs.get_pair_nums():
##    pairs.parse_SpecResult(pair)
##    print
    
##class Pairs(object):
##    def __init__(self,f,outdir,bt2,bdb=None,reference=None):
##        '''f is Path to primer blast file.'''
##        self.f = f
##        if outdir[-1] != "/":
##            outdir += "/"
##        self.outdir = outdir
##        self.pairs = {}
##        self.i = 1
##        self.info_paths = {}
##        self.primer_paths = {}
##        self.bt2 = bt2
##        self.blastdb = bdb
##        self.reference = reference
##        self.primerSpecOut_paths = {}
##        self.primerSpecResults = {}
##        self.productSpecOut_paths = {}
##        self.productSpecResults = {}
##        self.results = None
##        self.primer_test = False
##        self.product_test = False
##        
##    def readPrimerPair(self):#primer_blast
##        pair = {}
##        pair["name"] = self.fh.next()
##        pair["labels"] = self.fh.next()
##        pair["fwd"] = self.fh.next()
##        pair["rev"] = self.fh.next()
##        pair["length"] = self.fh.next()
##        return pair
##
##    def readPrimerPairs(self):#could take any function "readPrimerPair"
##        self.i = 1
##        self.fh = open(self.f, 'r')
##        while self.fh:
##            try:
##                pair = self.readPrimerPair()
##                self.pairs[self.i] = pair
##                self.i += 1
##            except:
##                break
##        self.fh.close()
##
##    def get_pair_nums(self):
##        return sorted(self.pairs.keys())
##
##    def get_fwd(self,pair):#primer_blast
##        return self.pairs[pair]["fwd"].split()[2]
##
##    def get_rev(self,pair):#primer_blast
##        return self.pairs[pair]["rev"].split()[2]
##
##    def get_primer_path(self,pair):
##        if not self.primer_paths:
##            self.write_primer_files()
##        return self.primer_paths[pair]
##
##    def get_pair_info_path(self, pair):#primer_blast
##        if not self.info_paths:
##            self.write_pair_stats()
##        return self.info_paths[pair]
##
##    def get_primerSpecResult(self,pair):
##        return self.primerSpecResults[pair]
##    
##    def write_pair_stats(self):#primer_blast
##        for pair in self.get_pair_nums():
##            path = self.outdir+"primer-pair-"+str(pair)+"-info.txt"
##            self.info_paths[pair] = path
##            o = open(path, 'w')
##            o.write(self.pairs[pair]["name"]+"\n")
##            o.write(self.pairs[pair]["labels"]+"\n")
##            o.write(self.pairs[pair]["fwd"]+"\n")
##            o.write(self.pairs[pair]["rev"]+"\n")
##            o.write(self.pairs[pair]["length"]+"\n")
##            o.write("\n")
##            o.close()
##
##    def write_primer_files(self):
##        for pair in self.get_pair_nums():
##            path_fwd = self.outdir+"primer-1-"+str(pair)+".fa"
##            path_rev = self.outdir+"primer-2-"+str(pair)+".fa"
##            self.primer_paths[pair] = [path_fwd,path_rev]
##            o = open(path_fwd, 'w')
##            o.write(">p1\n")
##            o.write(self.get_fwd(pair)+"\n")
##            o.close()
##            o = open(path_rev, 'w')
##            o.write(">p2\n")
##            o.write(self.get_rev(pair)+"\n")
##            o.close()
##
##    def primerSpecificity(self, verbose=False):
##        for pair in pairs.get_pair_nums():
##            if verbose:
##                sys.stderr.write("\t Pair " + str(pair) + "\n")
##            cmd = self._ps_cmd(pair)
##            result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
##            self.primerSpecResults[pair] = result
##        self.primer_test = True
##
##
##    def productSpecificity(self, verbose=False):
##        if self.blastdb is not None:
##            for pair in pairs.get_pair_nums():
##                if verbose:
##                    sys.stderr.write("\t Pair " + str(pair) + "\n")
##                cmd = self._target_cmd(pair)
##                result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
##                self.productSpecResults[pair] = result
##            self.product_test = True
##
##    def write_primerSpecResults(self):
##        if self.primer_test:
##            for pair in self.get_pair_nums():
##                path = self.outdir+"primer-pair-"+str(pair)+".specificity.txt"
##                self.primerSpecOut_paths[pair] = path
##                o = open(path, 'w')
##                o.write(self.primerSpecResults[pair]+"\n")
##                o.close()
##
##    def write_productSpecResults(self):
##        if self.product_test:
##            for pair in self.get_pair_nums():
##                path = self.outdir+"primer-pair-"+str(pair)+".product.txt"
##                self.productSpecOut_paths[pair] = path
##                o = open(path, 'w')
##                o.write(self.productSpecResults[pair]+"\n")
##                o.close()
##
##    def parse_primerSpecResult(self, pair):
##        if self.primer_test:
##            result = self.primerSpecResults[pair].split("\n")
##            i = 0
##            while result[i] != "FINAL REPORT:":
##                i+=1
##            parsed_result = []
##            for j in range(i,len(result)):
##                parsed_result.append(result[j])
##            parsed_result = Result(parsed_result)
##            return parsed_result
##
##    def parse_productSpecResult(self, pair):
##        if self.product_test:
##            result = self.productSpecResults[pair].split("\n")
##            i = 0
##            while result[i] != "PCR TARGET BLAST RESULTS":
##                i+=1
##            parsed_result = []
##            for j in range(i,len(result)):
##                parsed_result.append(result[j])
##            parsed_result = Result(parsed_result)
##            return parsed_result
##
##
##    def analyze_primerSpecResult(self,pair):
##        if self.primer_test:
##            f = self.parse_primerSpecResult(pair)
##            counts = []
##            f.next() # burn off "FINAL REPORT:"
##            f.next() #
##            counts.append( f.next().split()[4] )
##            counts.append( f.next().split()[4] )
##            counts.append( f.next().split()[4] )
##            counts.append( f.next().split()[4] )
##            f.next()
##            f.next()
##            counts.append( f.next().split()[2] )
##            counts.append( f.next().split()[2] )
##            f.next()
##            f.next()
##            counts.append( f.next().split()[2] )
##            counts.append( f.next().split()[2] )
##            return counts
##
##    def analyze_productSpecResult(self,pair):
##        if self.product_test:
##            f = self.parse_productSpecResult(pair)
##            counts = []
##            f.next() # burn off "PCR TARGET BLAST RESULTS"
##            counts.append( f.next().split()[5] )
##            return counts
##
##    def analyze_primerSpecResults(self):
##        if self.results is None:
##            self._initialize_results()
##        if self.primer_test:
##            for pair in self.get_pair_nums():
##                counts = self.analyze_primerSpecResult(pair)
##                self.results[pair]['primers'] = counts
##
##    def analyze_productSpecResults(self):
##        if self.results is None:
##            self._initialize_results()
##        if self.product_test:
##            for pair in self.get_pair_nums():
##                counts = self.analyze_productSpecResult(pair)
##                self.results[pair]['target'] = counts
##
##    def analyze_results(self):
##        if not self.primer_test:
##            self.primerSpecificity()
##        if not self.product_test:
##            self.productSpecificity
##        self.analyze_primerSpecResults()
##        self.analyze_productSpecResults()
##
##    def print_results(self):
##        for pair in self.get_pair_nums():
##            results = [str(pair), self.get_fwd(pair), self.get_rev(pair)] + [e for e in self.results[pair]['primers'] + self.results[pair]['target']]
##            print ("\t").join(results)
##
##    def _initialize_results(self):
##        self.results = {}
##        for pair in self.get_pair_nums():
##            self.results[pair] = {'primers':[], 'target':[]}
##            
##        
##
##    def _ps_cmd(self, pair):
##        return ["primerSpecificity.sh",  self.bt2, pairs.get_primer_path(pair)[0], pairs.get_primer_path(pair)[1]]
##
##    def _target_cmd(self, pair):
##        return ["productSpecificity.sh",  self.bt2, pairs.get_primer_path(pair)[0], pairs.get_primer_path(pair)[1], self.blastdb, self.reference]
##
##
## 
