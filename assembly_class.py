from Bio import SeqIO
import numpy as np
from helper_functions import *

class Assembly(object):
    def __init__(self, fastxpath, fastx="fastq"):
        self.assembly = {}
        self.fastxpath = fastxpath
        self.fastx = fastx
    def load_assembly(self):
        for line in SeqIO.parse(self.fastxpath, self.fastx):
            self.assembly[line.name] = line
    def add_to_assembly(self, contig):
        ''' contig is SeqIO fastx object
            --> can use this if it makes more sense to add assembly on the fly while computing other calculations'''
        self.assembly[line.name] = line
    def extract_info(self, contig, start, end):
        return self.assembly[contig][start:end]
    def is_correct_seq(self, contig, start, end, seq, makeupper=False):
        info = self.extract_info(contig, start, end)
        if makeupper:
            return str(info.seq).upper() == seq.upper()
        return str(info.seq) == seq
    def get_qv_stats(self, contig, start, end):
        if self.fastx == "fastq":
            info = self.extract_info(contig, start, end)
            qv_mean = np.mean(info.letter_annotations['phred_quality'])        
            qv_median = np.median(info.letter_annotations['phred_quality'])
            qv_min = np.min(info.letter_annotations['phred_quality'])
            qv_max = np.max(info.letter_annotations['phred_quality'])
            return qv_mean, qv_median, qv_min, qv_max
        return "-","-","-","-"
