#!/bin/bash

echo ASM $ASM
echo DBTYPE $DBTYPE
echo OUT $OUT
echo

makeblastdb -in $ASM -dbtype $DBTYPE -out $OUT
