#!/bin/bash

while read line; do echo $line; grep slurm $line; echo; done < <(find . -name "*.out")
while read line; do a=`grep -c slurm $line`; if [ $a -ge 1 ]; then echo $line; grep -A 2 slurm $line; echo; fi; done < <(find . -name "*.out")
