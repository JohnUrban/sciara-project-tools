#!/usr/bin/env python2.7
import sys
import argparse
import numpy as np
import re

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in bnx files, returns single line RMAP format (e.g. maligner input).
    
    RMAP format (columns):
    1. MAP NAME (STRING): A unique identifier for the map.
    2. MAP LENGTH (INT): The total length of the map, in base pairs.
    3. NUMBER OF FRAGMENTS (INT): The number of fragments in the map.
    4. FRAGMENT 1 (INT): The length of the first restriction fragment in the map, in base pairs.
    5. FRAGMENT 2 (INT): The length of the second restriction fragment in the map. ... and so on.
    .. .....
    .. .....
    N. Fragment X (INT): where X = N-4+1; length of Xth restriction fragment in the map. ... and so on.

    Fields 4 onwards provide the ordered listing of restriction fragment lengths in bp units.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-b', "--bnx",
                   type=str, required=True,
                   help='''BNX file from BioNano Irys''')

parser.add_argument("--snr",
                   type=int, default=0,
                   help='''minimum SNR to convert read. Default = 0.''')

parser.add_argument("--intensity",
                   type=int, default=0,
                   help='''Minimum intensity to convert read. Default = 0.''')


parser.add_argument("--bpp",
                   type=int, default=None,
                   help='''Minimum intensity to convert read. Default = 0.''')

args = parser.parse_args()



class RunInfo(object):
    def __init__(self):
        self.runinfo = {}
        self.run_IDs = []
        self.date_re_str = "[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}\s[0-9]{1,2}:[0-9]{1,2}:[0-9]{1,2}\s[A|P]M"
        self.date_re = re.compile(self.date_re_str)
        self.rundata_type = [int, float, float, int, str, int, str, float, float, int]
    def add_run_data(self, data):
        '''assumes line string already parsed and typed'''
        runID = data[-1]
        self.run_IDs.append(runID)
        self.runinfo[runID] = data
    def add_run_data_line(self, line):
        d = self.date_re.findall(line)[0]
        data = line.strip().split(d)[-1].strip().split()
        for i in range(len(data)):
            data[i] = self.rundata_type[i](data[i])
        self.add_run_data(data)
    def get_bpp(self, runID):
        return self.runinfo[runID][2]
    def get_stretch_factor(self, runID):
        return self.runinfo[runID][1]
    def get_number_of_scans(self, runID):
        return self.runinfo[runID][3]
    def get_run_IDs(self):
        return self.run_IDs

class BNX_Read(object):
    def __init__(self, channel_0, channel_1, QX11, QX12):
        ''' channel_0, channel_1, QX11, QX12 are just 4 unprocessed lines from bnx file describing the read'''
        self.ch0types = [int, int, float, float, float, int, int, int, int, str, int, int, int]
        self.__add_channel_0(channel_0)
        self.__add_channel_1(channel_1)
        self.__add_QX11(QX11)
        self.__add_QX12(QX12)
        
    def __add_channel_0(self, line):
        line = line.strip().split()
        assert len(line) == len(self.ch0types)
        self.moleculeID = int(line[1])
        self.length = float(line[2])
        self.avgintensity = float(line[3])
        self.snr = float(line[4])
        self.number_of_labels = int(line[5])
        self.original_molecule_ID = int(line[6])
        self.scan_number = int(line[7])
        self.scan_direction = int(line[8])
        self.chip_ID = str(line[9])
        self.flowcell = int(line[10])
        self.runID = int(line[11])
        self.global_scan_number = int(line[11])

    def __add_channel_1(self, line):
        line = [float(e) for e in line.strip().split()]
        self.label_positions = line[1:]
        assert len(self.label_positions) == self.number_of_labels + 1 ## last "label pos" is end of backbone and equals molecule length

    def __add_QX11(self,line):
        self.label_snr = [float(e) for e in line.strip().split()[1:]]
        assert len(self.label_snr) == self.number_of_labels

    def __add_QX12(self, line):
        self.label_intensity = [float(e) for e in line.strip().split()[1:]]
        assert len(self.label_intensity) == self.number_of_labels
        
    def get_number_of_labels(self):
        return self.number_of_labels
    def get_number_of_frags(self):
        return self.number_of_labels + 1
    def get_molecule_ID(self):
        return self.moleculeID
    def get_length(self):
        return self.length
    def get_avg_backbone_intensity(self):
        return self.avgintensity
    def get_backbone_snr(self):
        return self.snr
    def get_label_positions(self):
        return self.label_positions
    def get_label_snr(self):
        return self.label_snr
    def get_label_intensity(self):
        return self.label_intensity
    def get_run_ID(self):
        return self.runID
        
def bnx2rmap(read, runinfo):
    runid = read.get_run_ID()
##    bpp = runinfo.get_bpp(runid)
##    print bpp
##    print read.get_length()
    map_name = read.get_molecule_ID()
##    map_len_bp = int(round(read.get_length())) #int(bpp * read.get_length())
    num_frags = read.get_number_of_frags()
    frags = np.array(read.get_label_positions())
    lengths = np.zeros(num_frags)
##    print read.get_label_positions()
##    print frags
    assert num_frags == len(frags)
    lengths[0] = frags[0]
    for i in range(1,num_frags):
        lengths[i] = frags[i]-frags[i-1]
##    print lengths
    lengths = list(lengths)
##    print lengths
##    print read.get_length(), sum(lengths)
    for i in range(num_frags):
        lengths[i] = int(round(lengths[i])) #int( bpp * frags[i] )
##    print lengths
##    print map_len_bp, sum(lengths)
        map_len_bp = sum(lengths)
    return [map_name, map_len_bp, num_frags] + lengths
    
        
runinfo = RunInfo()

with open(args.bnx) as f:
    for line in f:
        if line[0] == "#": #header
            if "# BNX File Version" in line:
                file_version = float(line.strip().split()[-1])
            elif "# Label Channels:" in line:
                label_channel = int(line.strip().split()[-1])
            elif "# Bases per Pixel:" in line:
                global_bpp = float(line.strip().split()[-1])
            elif "#rh" in line:
                rh = line.strip().split()[4:]
            elif "# Run Data" in line:
                runinfo.add_run_data_line(line)
        else:
            ch0 = line
            ch1 = f.next()
            qx11 = f.next()
            qx12 = f.next()
            read = BNX_Read(ch0, ch1, qx11, qx12)
            rmap = bnx2rmap(read, runinfo)
            print ("\t").join([str(e) for e in rmap])
##            print sum(rmap[3:]), len(rmap[3:])

        
                
