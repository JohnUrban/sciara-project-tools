#!/bin/bash

echo BAM $BAM
echo BEDPE $BEDPE
echo MINSUPPORT ${MINSUPPORT}
echo

source ~/software/sniffles/source.sh
sniffles -m $BAM -b $BEDPE.bedpe --min_support $MINSUPPORT
