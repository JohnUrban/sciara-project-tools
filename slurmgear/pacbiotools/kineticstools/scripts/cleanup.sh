#!/bin/bash


date
echo DIR SPACE
diskUsageThisDirAndAllSubDir
echo NUM FILES
find . -name "*" | wc -l

echo MAIN DIR
ls
rm -r temp/
rm out_all.cmp.h5
echo
echo
echo CMPFILES
ls cmpfiles/
rm cmpfiles/*

echo DIR SPACE
diskUsageThisDirAndAllSubDir
echo NUM FILES
find . -name "*" | wc -l
echo MAINDIR
ls
echo cmpfiles
ls cmpfiles/
date
