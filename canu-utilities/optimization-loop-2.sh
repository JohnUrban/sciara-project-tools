#!/bin/bash

bin=~/software/canu/canu/Linux-amd64/bin
GKP=`readlink -f *.gkpStore`
OVL=`readlink -f *.ovlStore`
G1=210000000
G2=350000000
A=`grep repeatdetect 4-unitigger/unitigger.sh | awk '{print $2}'`
B1=`grep repeatdetect 4-unitigger/unitigger.sh | awk '{print $3}'`
X=`grep Guess 0-mercounts/*err | awk '{print $5}'`
B2=`echo $X*2 | bc` 
B3=`echo $X*3 | bc`
C=`grep repeatdetect 4-unitigger/unitigger.sh | awk '{print $4}'`

echo Initial Messages:
echo bin, $bin
echo GKP, $GKP
echo OVL, $OVL
echo G1, $G1
echo G2, $G2
echo X, $X
echo A, $A
echo B1, $B1
echo B2, $B2
echo B3, $B3
echo C, $C
echo
echo Starting Loop....
echo


for ovl in 500 1000 1500 2000 2500 3000 3500 4000 4500 5000; do
   for ee in 0.025 0.030 0.035 0.040 0.045 0.050 0.055 0.060 0.065 0.070 0.080 0.090 0.100 ; do
      for B in $B1 $B2 $B3; do 
        echo "ASM, OVL $ovl, EE $ee, B $B"
        mkdir 4-unitigger-RS-NS-CS-${ovl}-${ee}-${B} && cd 4-unitigger-RS-NS-CS-${ovl}-${ee}-${B} && echo `pwd` && $bin/bogart -G $GKP -O $OVL -T asm.tigStore -o asm -B 19824 -el ${ovl}  -eg ${ee} -eb 0.0625 -em 0.0375 -er ${ee} -threads 4 -M 16 -RS -NS -CS -unassembled 2 1000 0.75 0.75 2 -repeatdetect $A $B $C && cd ..
        echo $G1
        $bin/tgStoreDump -contigs -bubbles -sizes -s $G1 -G $GKP -T 4-unitigger-RS-NS-CS-${ovl}-${ee}-${B}/asm.tigStore 1 | awk -v "var=parameters::: G: $G1 ovl: $ovl ee: $ee B: $B" '{print $0 "\t" var}' #> 4-unitigger-RS-NS-CS-${ovl}-${ee}-${B}/${G1}.g1.out
        echo
        echo $G2
        $bin/tgStoreDump -contigs -bubbles -sizes -s $G2 -G $GKP -T 4-unitigger-RS-NS-CS-${ovl}-${ee}-${B}/asm.tigStore 1 | awk -v "var=parameters::: G: $G2 ovl: $ovl ee: $ee B: $B" '{print $0 "\t" var}' #> 4-unitigger-RS-NS-CS-${ovl}-${ee}-${B}/${G2}.g2.out
        echo; echo
      done
   done
done
