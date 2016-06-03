#!/bin/bash

pycalc () { python -c "print $1"; }

G=300e6; 
for f in g*/; 
 do echo $f; 
 N=`awkSum $f/correction/*.gkpStore/readlengths.txt`; 
 pycalc $N/$G; 
 N=`awkSum $f/trimming/*.gkpStore/readlengths.txt`; 
 pycalc $N/$G; 
 N=`awkSum $f/unitigging/*.gkpStore/readlengths.txt`; 
 pycalc $N/$G; 
 grep Guess $f/trimming/0-*/*err
 grep Guess $f/unitigging/0-*/*err; 
 echo; 
done
