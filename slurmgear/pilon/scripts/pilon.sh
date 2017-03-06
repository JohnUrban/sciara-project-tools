#!/bin/bash


##--export=NOSTRAYS=${NOSTRAYS},JX2=${JX2},PILONJAR=${PILONJAR},ASM=${ASM},READS=${}

echo NOSTRAYS, $NOSTRAYS
if $NOSTRAYS; then NO="--nostrays"; else NO=""; fi

echo "java -Xmx${JX} -jar $PILONJAR --genome $ASM --output $PRE --changes --frags ${READS} --diploid --fix ${FIX} ${NO}"


if $NOSTRAYS; then
    echo "Path1: nostrays=true, fix=$FIX"
    java -Xmx${JX} -jar $PILONJAR --genome $ASM --output $PRE --changes --frags ${READS} --diploid --fix $FIX --nostrays
else
    echo "Path2: nostrays=false, fix=$FIX"
    java -Xmx${JX} -jar $PILONJAR --genome $ASM --output $PRE --changes --frags ${READS} --diploid --fix $FIX
fi
