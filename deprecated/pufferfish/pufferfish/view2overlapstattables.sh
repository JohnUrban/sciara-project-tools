#!/bin/bash

F1=$1
F2=$2

paste <(awk 'OFS="\t" {print $1,$2,$3}' $F1) <(awk '{print $3,$4}' $F2) | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' 

