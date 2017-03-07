#!/bin/bash


if [ $# -eq 0 ] || [ $1 == "-h" ] || [ $1 == "-help" ] || [ $1 == "--help" ]; then
    echo "
    Usage: bash $0 names longtabl1 ... longtablN

    ...where names is the names.txt file in slurmgear/evalwrapper/utils dir.
    ...where longtables were produced by get-long-table.sh or one of its derivatives/relatives.
    "
    exit
fi

paste $@ | awk '{gsub(/\t/," \\& "); print $0 " \\\\"}'
