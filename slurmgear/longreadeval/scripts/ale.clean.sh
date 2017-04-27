#!/bin/bash

rm *.param
A=`echo *.ALE.txt`
head -n50 $A | grep "^#" > $A.1
mv $A.1 $A
