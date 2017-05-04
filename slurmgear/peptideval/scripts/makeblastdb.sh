#!/bin/bash

module load blast/2.2.30+ 

echo ASM $ASM
echo DBTYPE $DBTYPE
echo OUT $OUT
echo

makeblastdb -in $ASM -dbtype $DBTYPE -out $OUT
