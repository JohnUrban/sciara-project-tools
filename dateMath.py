#!/usr/bin/env python2.7
import argparse
import datetime as dt

parser = argparse.ArgumentParser(description="""

    Get days, weeks, months -- etc -- info from dates and questions.

    
    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('-d','--date',
                   type=str, default=None, 
                   help='''Provide date-of-interest as YYYY-MM-DD.''')

parser.add_argument('-b', '--basedate', 
                   type=str, default=None, 
                   help='''Second date to add or subtract from date given. Provide date-of-interest as YYYY-MM-DD. Defaults to current date.''')

parser.add_argument('-a','--add', 
                   type=int, default=0, 
                   help='''Tell me the date it will be in this many days. ''')

parser.add_argument('-s','--subtract', 
                   type=int, default=0, 
                   help='''Tell me the date it was this many days ago. ''')

parser.add_argument('-M','--minimaloutput', 
                   action='store_true', default=False, 
                   help='''. ''')

parser.add_argument('-S','--sciaraInfo', 
                   action='store_true', default=False, 
                   help='''Given a mating date, calculates estimated dates for
hatch/larval stage start, pupal stage start, adult stage start, approximate death dates, and cold room dates according to Jacob's method.
Only uses one date to calculate these with the following priority --date > --basedate > nothing.
In other words, --date or --basedate can be used, but using both is the same as using --date, and using neither specifies today's date.
''')


args = parser.parse_args()


def datemath(x, y):

    ans = str(x-y).split()
    ans2 = ans[2].split(":")

    weeks = float(ans[0])/7.0

    months =  float(ans[0])/30.0

    years =  float(ans[0]) / 365.0

    print "\tDifference between dates is", ans[0], ans[1], ans2[0], 'hours,', ans2[1], 'minutes,', ans2[2], 'seconds.'
    print "\tThat's roughly %f weeks or %f months or %f years." % (weeks, months, years)

def add(d, args, verbose=True):
    ans = d+dt.timedelta(days=args.add)
    if verbose:
        print "\tIn", args.add, "days from", d.strftime('%Y-%m-%d'), "it will be", ans.strftime('%Y-%m-%d')
    return ans

def sub(d, args, verbose=True):
    ans = d-dt.timedelta(days=args.subtract)
    if verbose:
        print "\tOn", args.subtract, "days before", d.strftime('%Y-%m-%d'), "it was", ans.strftime('%Y-%m-%d')
    return ans

def strToDate(x):
    l = [int(e) for e in x.strip().split('-')]
    return dt.datetime(l[0], l[1], l[2], 0, 0, 0)
        
def default(args):
    print 
    now = dt.datetime.now()
    if args.basedate is None:
        args.basedate = now.strftime('%Y-%m-%d')
        base = now
    else:
        #basel = [int(e) for e in args.basedate.strip().split('-')]
        #base = dt.datetime(basel[0], basel[1], basel[2], 0, 0, 0)
        base = strToDate(args.basedate)
    
    if args.date is not None:
        print "Given provided date and base date (today if not specified):"
        #datel = [int(e) for e in args.date.strip().split('-')]
        #date = dt.datetime(datel[0], datel[1], datel[2], 0, 0, 0)
        date = strToDate(args.date)
        datemath(date, base)
        if args.add or args.subtract:
            print "Given provided date and days to add/sub:"
        if args.add > 0:
            add(date, args)
        if args.subtract > 0:
            sub(date, args)


    if args.add or args.subtract:
            print "Given today's date and days to add/sub:"
    if args.add > 0:
        add(base, args)

    if args.subtract > 0:
        sub(base, args)
    print


def dateDecision(args):
    if args.date is not None:
        DATE = strToDate(args.date)
    elif args.basedate is not None:
        DATE = strToDate(args.basedate)
    else:
        DATE = dt.datetime.now()

    return DATE


def datemsg(msg, date0, date1, days, duration, sfx=''):
    print '\t'.join([str(e) for e in [msg, date0.strftime('%Y-%m-%d'), date1.strftime('%Y-%m-%d'), str(days)+sfx, str(duration)+sfx]])

def dadd(d, x):
    return d+dt.timedelta(days=x)

#pre-24
#ee-25
#8-26
#10-27
#12-28
#14-29
#EE30
#DJ31
def sciara(args, hatch0=7, hatch1=10, larval0=7, larval1=31,
           earlyLarval0=7, earlyLarval1=20, lateLarval0=21, lateLarval1=31,
           eyespot0=22, eyespot1=25, pupal0=28, pupal1=36, adult0=35,
           adult1=43, cold0=16, cold1=22, post0=8, post1=20):
    NOW=dt.datetime.now()
    DATE = dateDecision(args)
    print '\t'.join(['Event/Stage', 'EarliestSeen', 'LatestSeen', 'DPM', 'Duration'])
    datemsg('Today_Date',
            NOW,
            NOW,
            '-',
            '-')
    datemsg('Mating_Date',
            DATE,
            DATE,
            '-',
            '-')
    datemsg('Hatch_Date',
            dadd(DATE, hatch0),
            dadd(DATE, hatch1),
            '-'.join(str(e) for e in [hatch0, hatch1]),
            '-')
    datemsg('Larval_stages',
            dadd(DATE, larval0),
            dadd(DATE, larval1),
            '-'.join(str(e) for e in [larval0, larval1]),
            '21')

    datemsg('Instar_1-3',
            dadd(DATE, earlyLarval0),
            dadd(DATE, earlyLarval1),
            '-'.join(str(e) for e in [earlyLarval0, earlyLarval1]),
            '14')
    datemsg('Instar_4',
            dadd(DATE, lateLarval0),
            dadd(DATE, lateLarval1),
            '-'.join(str(e) for e in [lateLarval0, lateLarval1]),
            '7-10')
    datemsg('Early_eyepot',
            dadd(DATE, eyespot0),
            dadd(DATE, eyespot1),
            '-'.join(str(e) for e in [eyespot0, eyespot1]),
            1)
    datemsg('Eyespot_8x4',
            dadd(DATE, eyespot0+1),
            dadd(DATE, eyespot1+1),
            '-'.join(str(e+1) for e in [eyespot0, eyespot1]),
            1)
    datemsg('Eyespot_10x5',
            dadd(DATE, eyespot0+2),
            dadd(DATE, eyespot1+2),
            '-'.join(str(e+2) for e in [eyespot0, eyespot1]),
            1)
    datemsg('Eyespot_12x6',
            dadd(DATE, eyespot0+3),
            dadd(DATE, eyespot1+3),
            '-'.join(str(e+3) for e in [eyespot0, eyespot1]),
            1)
    datemsg('Eyespot_14x7',
            dadd(DATE, eyespot0+4),
            dadd(DATE, eyespot1+4),
            '-'.join(str(e+4) for e in [eyespot0, eyespot1]),
            1)
    datemsg('Edge_Eye',
            dadd(DATE, eyespot0+5),
            dadd(DATE, eyespot1+5),
            '-'.join(str(e+5) for e in [eyespot0, eyespot1]),
            1)
    datemsg('Drop_Jaw',
            dadd(DATE, eyespot0+6),
            dadd(DATE, eyespot1+6),
            '-'.join(str(e+6) for e in [eyespot0, eyespot1]),
            1)
    
    datemsg('Pupal_stages',
            dadd(DATE, pupal0),
            dadd(DATE, pupal1),
            '-'.join(str(e) for e in [pupal0, pupal1]),
            5)
    datemsg('Adult_stages',
            dadd(DATE, adult0),
            dadd(DATE, adult1),
            '-'.join(str(e) for e in [adult0, adult1]),
            adult1-adult0+1)
    
    datemsg('Cold_room',
            dadd(DATE, cold0),
            dadd(DATE, cold1),
            '-'.join(str(e) for e in [cold0, cold1]),
            cold1-cold0+1)
    ## Adult stage if date given was date taken from cold room
    datemsg('Adult_cold',
            dadd(DATE, post0),
            dadd(DATE, post1),
            '-'.join(str(e) for e in [post0, post1]),
            '8-20')




def run(args):
    if args.sciaraInfo:
        sciara(args)
    else:
        default(args)

    





#### EXECUTE
run(args)
