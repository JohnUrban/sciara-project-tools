#!/bin/bash

PATH=~/software/sciaratools/sciara-project-tools/alenano/src/:$PATH
PATH=~/software/alenano/src/:$PATH

ALE ${BAM} ${REF} ${BASE}.ALE.txt >> ${BASE}.ale.err
