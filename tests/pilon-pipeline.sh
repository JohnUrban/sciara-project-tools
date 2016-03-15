#!/bin/bash

if [ $# -eq 0 ]; then echo "
Usage:
sh pilon-pipeline.sh pathto_pilon-pipeline.py
"
exit
fi


echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --bam mappedreads.bam"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --bam mappedreads.bam
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --bam mappedreads.bam"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --bam mappedreads.bam
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --bam mappedreads.bam --mkdup_bam mappedreads.mkdup.bam"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --bam mappedreads.bam --mkdup_bam mappedreads.mkdup.bam
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --mkdup_bam mappedreads.mkdup.bam"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --R1 embryoR1.fq --R2 embryoR2.fq --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --mkdup_bam mappedreads.mkdup.bam
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --mkdup_bam mappedreads.mkdup.bam"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --flagstat --mkdupjar mkdup.jar --clean --bt2 bt2index --mkdup_bam mappedreads.mkdup.bam
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --flagstat --mkdupjar mkdup.jar --clean --mkdup_bam mappedreads.mkdup.bam"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --flagstat --mkdupjar mkdup.jar --clean --mkdup_bam mappedreads.mkdup.bam
echo
echo
echo
echo "python $1 --dry --pilonjar pilon.jar --reference fake.fasta --mkdupjar mkdup.jar --clean --mkdup_bam mappedreads.mkdup.bam"
python $1 --dry --pilonjar pilon.jar --reference fake.fasta --mkdupjar mkdup.jar --clean --mkdup_bam mappedreads.mkdup.bam
