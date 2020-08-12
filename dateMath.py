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


args = parser.parse_args()


def datemath(x, y):

    ans = str(x-y).split()
    ans2 = ans[2].split(":")

    weeks = float(ans[0])/7.0

    months =  float(ans[0])/30.0

    years =  float(ans[0]) / 365.0

    print "\tDifference between dates is", ans[0], ans[1], ans2[0], 'hours,', ans2[1], 'minutes,', ans2[2], 'seconds.'
    print "\tThat's roughly %f weeks or %f months or %f years." % (weeks, months, years)

def add(d, args):
    ans = d+dt.timedelta(days=args.add)
    print "\tIn", args.add, "days from", d.strftime('%Y-%m-%d'), "it will be", ans.strftime('%Y-%m-%d')
    
def sub(d, args):
    ans = d-dt.timedelta(days=args.subtract)
    print "\tOn", args.subtract, "days before", d.strftime('%Y-%m-%d'), "it was", ans.strftime('%Y-%m-%d')
    

def run(args):
    print 
    now = dt.datetime.now()
    if args.basedate is None:
        args.basedate = now.strftime('%Y-%m-%d')
        base = now
    else:
        basel = [int(e) for e in args.basedate.strip().split('-')]
        base = dt.datetime(basel[0], basel[1], basel[2], 0, 0, 0)
    
    if args.date is not None:
        print "Given provided date and base date (today if not specified):"
        datel = [int(e) for e in args.date.strip().split('-')]
        date = dt.datetime(datel[0], datel[1], datel[2], 0, 0, 0)
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

    





#### EXECUTE
run(args)
