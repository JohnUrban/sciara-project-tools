#!/bin/bash

if [ $# -eq 0 ]; then echo "
Usage:
primerStockCalculator.sh nmol
-- where nmol is how many nmols come in dry primer tube
-- this script just tells you how much ultra pure water to add to make 100 uM
-- this just calculates for 100 uM
-- Note: that a simple trick for a human to do is just moving the decimal point in nmol over to the right 1 decimal place
-- i.e. multiply nmol by 10
-- That gives how much UPW to add in ul
-- e.g.: 17.6 nmol --> 176 ul
"
exit
fi

A=`echo $1*0.000000001 | bc -l`
B=`echo 100*0.000001 | bc -l`
echo 1000000*$A/$B | bc -l
echo 10*$1 | bc -l
