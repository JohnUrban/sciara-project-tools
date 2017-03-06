#!/bin/bash


##--export=NOSTRAYS=${NOSTRAYS},JX2=${JX2},PILONJAR=${PILONJAR},ASM=${ASM},READS=${}

echo NOSTRAYS, $NOSTRAYS
echo "java -Xmx${JX} -jar $PILONJAR --genome $ASM --output $PRE --changes --frags ${READS} --diploid --fix bases --nostrays"

if $NOSTRAYS; then
    java -Xmx${JX} -jar $PILONJAR --genome $ASM --output $PRE --changes --frags ${READS} --diploid --fix bases --nostrays
else
    java -Xmx${JX} -jar $PILONJAR --genome $ASM --output $PRE --changes --frags ${READS} --diploid --fix bases
fi
