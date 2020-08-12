#!/usr/local/bin/python
import sys
##import numpy as np
## adapted from MACS2


def binomial_cdf_inv ( cdf, a, b ):
    """BINOMIAL_CDF_INV inverts the binomial CDF. For lower tail only!

    """
    if cdf < 0 or cdf >1:
        raise Exception("CDF must >= 0 or <= 1")

    cdf2 = 0

    for x in xrange(0,a+1):
        pdf = binomial_pdf (x,a,b)
        cdf2 += pdf
        if cdf < cdf2:
            return x
    return a


def binomial_pdf( x, a, b ):
    """binomial PDF by H. Gene Shin
    
    """
    
    if a<1:
        return 0
    elif x<0 or a<x:
        return 0
    elif b==0:
        if x==0:
            return 1
        else:
            return 0
    elif b==1:
        if x==a:
            return 1
        else:
            return 0

    if x>a-x:
        p=1-b
        mn=a-x
        mx=x
    else:
        p=b
        mn=x
        mx=a-x
    pdf=1
    t = 0
    for q in xrange(1, mn+1):
        pdf *= (a-q+1)*p/(mn-q+1)
        if pdf < 1e-100:
            while pdf < 1e-3:
                pdf /= 1-p
                t -= 1
        if pdf > 1e+100:
            while pdf > 1e+3 and t < mx:
                pdf *= 1-p
                t += 1

    for i in xrange(mx-t):
        pdf *= 1-p
        
    pdf=float("%.10e" % pdf)
    return pdf


def cal_max_dup_tags ( genome_size, tags_number, p=1e-5 ):
    """Calculate the maximum duplicated tag number based on genome
    size, total tag number and a p-value based on binomial
    distribution. Brute force algorithm to calculate reverse CDF no
    more than MAX_LAMBDA(100000).
    
    """
    return binomial_cdf_inv(1-p,tags_number,1.0/genome_size)





genome_size = float(sys.argv[1])
tags_number = int(sys.argv[2])
p = float(sys.argv[3])
maxtags=cal_max_dup_tags (genome_size, tags_number, p)



## idea...
## for each bin of some width...
## obtain 5pe info until reach end of bin
## 
