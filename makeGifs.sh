#!/bin/bash

if [ $# -eq 0 ] ; then echo -e "\nmakeGIF.sh delay out.gif input1 .... inputN\n\tDelay = integer - e.g. 10\n\tout.gif = output filename with gif extension\n\tInputs = collection of PNG of JPGs.\n" ; exit ; fi

DELAY=${1}
OUTGIF=${2}
INPUTS=${@:3}

for VAR in DELAY OUTGIF INPUTS; do echo ${VAR} ${!VAR} ; done

convert -delay ${DELAY} -loop 0 ${INPUTS} ${OUTGIF}
