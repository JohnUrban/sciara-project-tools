#!/bin/bash

source /users/jurban/data/software/hisat/source.sh

echo G $G
echo PRE $PRE
echo

hisat2-build $G $PRE
