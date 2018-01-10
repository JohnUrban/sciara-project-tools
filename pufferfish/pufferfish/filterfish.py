import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg
from normalize import protocol1, protocol2, protocol3, protocol4, protocol5, protocol6, normalize




def run(parser, args):
    if args.protocol1:
        protocol=1
    elif args.protocol2:
        protocol=2
    elif args.protocol3:
        protocol=3
    elif args.protocol4:
        protocol=4
    elif args.protocol5:
        protocol=5
    elif args.protocol6:
        protocol=6
    
    if not args.skipnorm:
        late = normalize(latestage=args.latestage, protocol=protocol, earlystage=args.earlystage, pseudo=args.pseudo, bandwidth=args.bandwidth, quiet=args.quiet)
    else:
        late = CovBed(args.latestage)
    
    if args.counts:
        if not args.quiet:
            bdgmsg("final normalized late stage counts", args.collapsed)
        o = open(args.counts + ".bedGraph", 'w')
        o.write(late.get_bdg(late.count, args.collapsed))
        o.close()

        

    if args.stdev_above: ## want relation in terms of SD units
        sd = late.get_sd()
        value = late.get_mean() + args.value * sd
    elif args.stdev_below: ## want relation in terms of SD units
        sd = late.get_sd()
        value = late.get_mean() - args.value * sd
    elif args.mean:  ## want relation in terms MU units
        mu = late.get_mean()
        value = args.value * mu
    else: ## want relation wrt given value 
        value = args.value

    if not args.quiet:
            bdgmsg("filtered bdg values s.t. col4 is " + args.relation + " " + str(value) +"...", collapsed=False)
    if args.stdev_above or args.stdev_below or args.mean:
        newmsg("Mean = " + str(late.get_mean()))
        newmsg("SD = " + str(late.get_sd()))
        
    sys.stdout.write(late.filtered_bdg(relation = args.relation, value = value, bdg=None))
        
