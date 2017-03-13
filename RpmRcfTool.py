#!/usr/bin/env python
import sys
import argparse


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Compute RCF from RPM or RPM from RCF.
    Output rounded to nearest int.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--rpm', '-rpm',
                   type=str, default=False,
                   help='''Provide single RPM integer value OR comma-separated range of integer values as min,max,stepsize (e.g. 1600,13000,100)''')

parser_input.add_argument('--rcf', '-rcf',
                   type=str, default=False,
                   help='''Provide single RCF integer value OR comma-separated range of integer values as min,max,stepsize ''')

parser.add_argument('--radius', '-r',
                   type=float, default=75, required=True,
                   help='''Provide radius of rotor in mm. Default is 75 mm.''')

parser.add_argument('--range', '-R',
                   type=list, default=False, required=False,
                   help='''Provide range numbers: ''')


args = parser.parse_args()


### functions
def rpm2rcf(rpm, radius):
    rcf = 1.12*radius*(rpm/1000.0)**2
    return int(round(rcf,0))

def rcf2rpm(rcf, radius):
    rpm = 1000.0*(rcf/(1.12*radius))**0.5
    return int(round(rpm,0))

def answer(arg, function, radius):
    ''' Provide args.rpm or args.rcf as arg, rpm2rcf or rcf2rpm as function, and args.radius as radius'''
    inputlist = [int(e) for e in arg.split(",")]
    if len(inputlist) == 1:
        print function(inputlist[0], radius)
    elif len(inputlist) == 3:
        for inputval in range(inputlist[0], inputlist[1], inputlist[2]):
            print inputval, function(inputval, radius)

### EXECUTE   
if args.rpm:
    answer(args.rpm, rpm2rcf, args.radius)
elif args.rcf:
    answer(args.rcf, rcf2rpm, args.radius)
    
