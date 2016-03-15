import sys
import numpy as np
from Bio import SeqIO
from helper_functions import *

class Change(object):
    def __init__(self, line):
        self.change = line.strip().split()
        self.change_type = None
        self.old_len = None
        self.new_len = None
    def get_change(self):
        return self.change
    def get_coord(self):
        contig = self.change[0].split(":")[0]
        start = int(self.change[0].split(":")[-1].split("-")[0])
        end = int(self.change[0].split(":")[-1].split("-")[-1])
        ##make pythonese
        start -= 1
        return contig, start, end
    def get_old_seq(self):
        return self.change[2]
    def get_new_seq(self):
        return self.change[3]
    def get_change_type(self):
        if self.change_type == None:
            self._get_change_type()
        return self.change_type
    def get_old_len(self):
        if self.old_len == None:
            self._get_old_len()
        return self.old_len
    def get_new_len(self):
        if self.new_len == None:
            self._get_new_len()
        return self.new_len
    ################################################
    ###### '''internal functions''' ################
    ################################################
    def _get_change_type(self):
        if self.change[2] == ".":
            self.change_type = "ins"
        elif self.change[3] == ".":
            self.change_type = "del"
        else:
            self.change_type = "sub"
    def _get_old_len(self):
        if self.change[2] == ".":
            self.old_len = 0
        else:
            self.old_len = len(self.change[2])
    def _get_new_len(self):
        if self.change[3] == ".":
            self.new_len = 0
        else:
            self.new_len = len(self.change[3])



##class Changes(object):
##    def __init__(self):
##        self.changes = {}
##        self.numchanges = 0
##        self.change_types = {'sub':0, 'ins':0, 'del':0}
##    def add_change(self, change):
##        ''' change is of class Change'''
##        self.numchanges += 1
##        self.changes[self.numchanges] = change
##        self.change_types[change.get_change_type()] += 1
##
##class Change_Stats(object):
##    def __init__(self):
##        self.numchanges = 0
##        self.change_types = {'sub':{}, 'ins':{}, 'del':{}}
##        for key in self.change_types.keys():
##            self.change_types[key]['count'] = 0
##            self.change_types[key]['oldlen'] = []
##            self.change_types[key]['newlen'] = []
##        self.subs = {'len':[]}
##        self.ins = {'len':[]}
##        self.dels = {'len':[]}
##    def add_change(self, change, assembly):
##        ''' change is of class Change'''
##        ''' assembly is of class Assembly - with assembly pre-loaded: assembly.load_assembly()'''
##        self.numchanges += 1
##        ch_type = change.get_change_type()
##        self.change_types[ch_type]['count'] += 1
##        self.change_types[ch_type]['oldlen'].append( change.get_old_len() )
##        self.change_types[ch_type]['newlen'].append( change.get_new_len() )
        
        
        
def characterize_changes(changes, assembly):
    ## line by line -- not loaded into memory like Changes() class approach
    ''' changes = path to pilon.changes file'''
    ''' assembly = Assembly Class object - with assembly pre-loaded: assembly.load_assembly()'''
    for line in changes:
        change = Change(line)
        seq = change.get_old_seq()
        oldlen = change.get_old_len()
        newlen = change.get_new_len()
        change_type = change.get_change_type()
        contig, start, end = change.get_coord()
        qv_mean, qv_median, qv_min, qv_max = assembly.get_qv_stats(contig, start, end) ##note the values for insertions are not of what was inserted.. but the area its inserted in
        rc = "-"
        if change_type == "sub":
            if seq == revcomp(change.get_new_seq()) and oldlen > 1:
                rc = "True"
            else:
                rc = "False"
        out = ("\t").join(change.get_change() + [str(e) for e in [change_type, oldlen, newlen, qv_mean, qv_median, qv_min, qv_max, rc]])
        sys.stdout.write(out+"\n")
