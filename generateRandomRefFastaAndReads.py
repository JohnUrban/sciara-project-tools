#!/usr/bin/env python2.7
import sys, random, argparse
import numpy as np
from makeRandomSeqTools import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION

    Don't use this.

    So many better tools!

    I just needed this for a moment.

        

    """,
    formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-O', '--refpre',
                    type=str, default='reference',
                    help='''Default reference.''')

parser.add_argument('-o', '--readpre',
                    type=str, default='reads',
                    help='''Default reads.''')

parser.add_argument('-L', '--reflen',
                    type=int, default=100000,
                    help='''Default 100000.''')
parser.add_argument('-N', '--nrefs',
                    type=int, default=1,
                    help='''Default 1.''')

parser.add_argument('-n', '--nreads',
                    type=int, default=10,
                    help='''Default 10.''')

parser.add_argument('-m', '--meanreadlen',
                    type=int, default=5000,
                    help='''Default 5000.''')

parser.add_argument('-s', '--sdreadlen',
                    type=int, default=500,
                    help='''Default 500.''')

parser.add_argument('-l', '--minreadlen',
                    type=int, default=500,
                    help='''Default 500.''')

parser.add_argument('-X', '--mismatchrate', type=float, default=0, help=''' 0-1''')
parser.add_argument('-D', '--delrate', type=float, default=0, help=''' 0-1''')
parser.add_argument('-I', '--insrate', type=float, default=0, help=''' 0-1''')

parser.add_argument('-DEL', '--lgdelrate', type=float, default=0, help=''' 0-1''')
parser.add_argument('-INS', '--lginsrate', type=float, default=0, help=''' 0-1''')
parser.add_argument('-INV', '--lginvrate', type=float, default=0, help=''' 0-1''')
parser.add_argument('-DUP', '--lgduprate', type=float, default=0, help=''' 0-1''')
parser.add_argument('-TRA', '--trarate', type=float, default=0, help=''' 0-1''')

parser.add_argument("-A", "--alphabet",
                   type=str, default='A,C,G,T',
                   help='''Comma-separated letters to make barcode with. Default: A,C,G,T''')
parser.add_argument("-w", "--nonself",
                   type=int, default=3,
                   help='''Nonself weight -- given last letter was b, higher weights give lower probablities that next letter will be b. Default: 2. w=1 gives uniform. w=0 gives homopolymers of length L.''')

parser.add_argument("-W", "--weights",
                    type=str, default=False,
                    help=''' Provide comma-separated list of integer weights for each letter in alphabet. If not used, they are weighted uniformly.''')

args = parser.parse_args()




alphabet = args.alphabet.strip().split(",")

weights = {k:1 for k in alphabet}

if args.weights:
    args.weights = [int(e) for e in args.weights.strip().split(",")]
    for i in range(len(alphabet)):
        weights[alphabet[i]] = args.weights[i]



init = alphabet
with open(args.refpre+'.fasta','w') as f:
    d={}
    reflen = {}
    for i in range(args.nrefs):
        f.write( '>'+args.refpre+'_'+str(i)+'\n')
        seq = init[random.randint(0,len(init)-1)]
        model = getmodel(seq[-1],args.nonself,alphabet,weights)
        for j in range(1,args.reflen):
            seq += model[random.randint(0,len(model)-1)]
            model = getmodel(seq[-1],args.nonself,alphabet,weights)
        f.write(seq+'\n')
        d[i] = seq
        reflen[i] = len(seq)



tra_possible = len(d.keys()) > 1


with open(args.readpre+'.fasta','w') as f:
    for i in range(args.nreads):
        rlen = int( np.random.normal(args.meanreadlen,args.sdreadlen) )
        while rlen < args.minreadlen:
            rlen = int( np.random.normal(args.meanreadlen,args.sdreadlen) )

        # Add large events
        refsel = np.random.choice(d.keys())
        minsvsize = 200
        SV = ''
        #sys.stderr.write(str(reflen[refsel])+'\n')
        if np.random.binomial(1,args.lgdelrate):
            largest = int(reflen[refsel]*0.5)
            size = np.random.randint(minsvsize, largest)
            finalstart = int(reflen[refsel] - size - (0.01*rlen+1)) # should give just a tiny bit of wiggle room for the breakpoint to actually be on the read in edge case
            #finalstart = reflen[refsel] - min(minsvsize*3, reflen[refsel]-minsvsize)
            svstart = int( np.random.uniform(0, finalstart) )
            #largest = max(minsvsize, reflen[refsel]-finalstart-minsvsize)
            #size = np.random.randint(minsvsize, largest)
            #sys.stderr.write(str(finalstart) + ' ' + str(svstart) + ' ' + str(largest)+' ' + str(size) + '\n')
            svend = svstart + size
            newref = d[refsel][:svstart] + d[refsel][svend:]
            newreflen = len(newref)
            newrlen = min(rlen, newreflen)
            fudge = int(max(newrlen-200, newrlen*0.9))
            relearliestlstart = min(max(0, svstart-fudge), newreflen-newrlen)
            rellastend = min(newreflen, svstart+fudge)
            start = int( np.random.uniform( relearliestlstart, max(relearliestlstart,svstart-1) ) )
            end = start + newrlen
            breakpoint_in_read = svstart-start
            remainder = newrlen - breakpoint_in_read
            read = newref[start:end]
            rlen = newrlen
            assert svend-svstart == size
            SV += "DEL:"+str(svstart)+"-"+str(svend)+"_len:"+str(svend-svstart)+"_bpir/LHS:"+str(breakpoint_in_read)+'_RHSlen:'+str(remainder)
        elif np.random.binomial(1,args.lginsrate):
            ## for now I am limiting lg insertions to be 50% RL
            ##  and be 100% internal to read only
            largest = int(0.5*rlen)
            minsvsize = min(minsvsize, largest)
            size = np.random.randint(minsvsize, largest)
            refrlen = rlen-size
            start = int( np.random.uniform(0, len(d[refsel])-rlen) )
            end = start + refrlen
            refread = d[refsel][start:end]
            insertref = np.random.choice(d.keys())
            insert_start = int( np.random.uniform(0, reflen[insertref]-size) )
            insert_end = insert_start+size
            insert = d[insertref][insert_start:insert_end]
            breakpoint_in_read = np.random.randint(0, refrlen)
            refstart = start+breakpoint_in_read
            #sys.stderr.write(str(breakpoint_in_read)+'\n')
            read = refread[:breakpoint_in_read] + insert + refread[breakpoint_in_read:]
            SV += "INS:"+str(refstart-1)+"-"+str(refstart+1)+'from:'+str(insertref)+":"+str(insert_start)+"-"+str(insert_end)+"_len:"+str(size)+"_bpir/LHS:"+str(breakpoint_in_read)+'_RHSlen:'+str(breakpoint_in_read+size)
        elif np.random.binomial(1,args.lginvrate):
            #largest = int(reflen[refsel]*0.5) ## This allows too many reads to be inside the INV completely
            largest = int(rlen*0.95)
            size = np.random.randint(minsvsize, largest)
            finalstart = reflen[refsel] - size
            svstart = int( np.random.uniform(0, finalstart) )
            svend = svstart + size
            newref = d[refsel][:svstart] + revcomp(d[refsel][svstart:svend]) + d[refsel][svend:]
            fudge = int(max(rlen-200, rlen*0.9))
            earlieststart = min(max(0, svstart-fudge), reflen[refsel]-rlen)
            fudge = int(min(200, rlen*0.1))
            laststart = min(reflen[refsel]-rlen, svend-fudge)
            #sys.stderr.write(str(earlieststart)+' '+str(laststart)+'\n')
            start = np.random.randint( earlieststart, laststart ) 
            end = start + rlen
            if start > svstart:
                bp1 = '5pToRead'
            else:
                bp1 = svstart - start
            if end < svend:
                bp2 = '3ptoRead'
            else:
                bp2 = svend - start
            read = newref[start:end]
            assert svend-svstart == size
            SV += "INV:"+str(svstart)+"-"+str(svend)+"_len:"+str(size)+"_bp1-2_wrtr:"+str(bp1)+'-'+str(bp2)
        elif np.random.binomial(1,args.lgduprate):
            ## WORKING ON 20190214
            largest = int(rlen*0.95)
            size = np.random.randint(minsvsize, largest)
            finalstart = reflen[refsel] - size
            svstart = int( np.random.uniform(0, finalstart) )
            svend = svstart + size
            newref = d[refsel][:svstart] + 2*d[refsel][svstart:svend] + d[refsel][svend:]
            newreflen = len(newref)
            fudge = int(max(rlen-200, rlen*0.9))
            earlieststart = min(max(0, svstart-fudge), newreflen-rlen)
            fudge = int(min(200, rlen*0.1))
            laststart = min(newreflen-rlen, svend+size-fudge)
            #sys.stderr.write(str(earlieststart)+' '+str(laststart)+'\n')
            start = np.random.randint( earlieststart, laststart ) 
            end = start + rlen
            if start > svstart:
                bp1 = '5pToRead'
            else:
                bp1 = svstart - start
            if end < svend+size:
                bp2 = '3ptoRead'
            else:
                bp2 = svend + size - start
            read = newref[start:end]
            #assert svend-svstart == size
            SV += "DUP:"+str(svstart)+"-"+str(svend)+"_len:"+str(size)+"_bp1-2_wrtr:"+str(bp1)+'-'+str(bp2)

        elif np.random.binomial(1,args.trarate) and tra_possible:
            pass
        else:
            start = int( np.random.uniform(0, len(d[refsel])-rlen) )
            end = start+rlen
            read = d[refsel][start:end]

        ## Add small "technology" errors
        mutread = ''
        cnt = {'M':0, 'X':0, 'D':0, 'I':0}
        errors = 0
        alnlen = 0
        for b in read:
            newb = b
            error = 0
            aln = 1
            if np.random.binomial(1,args.mismatchrate):
                ntweights = {k:v for k,v in weights.iteritems()}
                ntweights[b] = 0
                model = getmodel(b,1,alphabet,ntweights)
                newb = model[random.randint(0,max(len(model)-1,0))]
                cnt['X']+=1
                error = 1
            if np.random.binomial(1,args.delrate):
                newb = ''
                cnt['D']+=1
                error = 1
            if np.random.binomial(1,args.insrate):
                model = getmodel(b,1,alphabet,weights)
                if len(newb) > 0:
                    aln = 2
                newb = newb + model[random.randint(0,len(model)-1)]
                cnt['I']+=1
                error = 1
            if not error:
                newb = b
                cnt['M']+=1
            mutread += newb
            errors += error
            alnlen += aln
        finalrlen = len(mutread)
        name = args.readpre +'_' + str(i) + '\tref:' + str(refsel)+';start:'+str(start)
        name += ';rlen_init:'+str(rlen)+';rlen_final:'+str(finalrlen)
        name += ';' + SV
        approx_pctsim = cnt['M']/float(sum(cnt.values()))
        alt_approx_pctsim = (alnlen-errors)/float(alnlen)
        name += ';' + 'approx_pctsim:' + str(approx_pctsim) + ';' + 'alt_approx_pctsim:' + str(alt_approx_pctsim)
        for key in cnt.keys():
            name += ';' + key + ':' + str(cnt[key])
        # could also add MXDI counts, orig read len, new readlen, etc

        f.write('>'+name+'\n'+mutread+'\n')
        






## ALT FOR INSERTION
##        elif np.random.binomial(1,args.lginsrate):
##            ## for now I am limiting lg insertions to be 90% RL
##            ##  and be 100% internal to read only
##            largest = int(0.9*rlen)
##            minsvsize = min(minsvsize, largest)
##            size = np.random.randint(minsvsize, largest)
##            #rlen = rlen-size
##            #start = int( np.random.uniform(0, len(d[refsel])-rlen) )
##            #end = start+rlen
##            #read = d[refsel][start:end]
##            insertref = np.random.choice(d.keys())
##            insert_start = int( np.random.uniform(0, reflen[insertref]-size) )
##            insert_end = insert_start+size
##            insert = d[insertref][insert_start:insert_end]
##            #breakpoint_in_read =
##            svstart = int( np.random.uniform(0, len(d[refsel])-minsvsize) )
##            newref = d[refsel][:svstart] + insert + d[refsel][svstart:]
##            newreflen =len(newref)
##            earliest_start = max(0, svstart-readlen+10)
##            latest_start = max(earliest_start, min(svstart-size-10, newreflen-rlen))
##            start = int( np.random.uniform(earliest_start, latest_start) )
##            end = start+rlen
##            read = newref[start:end]
##            breakpoint_in_read = svstart-start
##            remainder = newrlen - breakpoint_in_read
##            SV += "INS:"+str(svstart-1)+"-"+str(svstart+1)+"_len:"+str(size)+"_bpir/LHS:"+str(breakpoint_in_read)+'_RHSlen:'+str(remainder)

