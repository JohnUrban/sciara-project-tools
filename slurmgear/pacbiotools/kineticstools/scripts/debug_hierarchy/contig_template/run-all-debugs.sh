#!/bin/bash

SUBDIRRUN=run-kinetics-contig-nonKpipe-debug-SUBDIRs.sh

bash create-config.sh

echo debug1
bash run-kinetics-contig-nonKpipe-debug.sh

for SUBDIR in debug*/; do
  echo ${SUBDIR}
  if [ -d ${SUBDIR} ]; then
    cp ${SUBDIRRUN} ${SUBDIR}/
    cd ${SUBDIR}
    bash ${SUBDIRRUN}
    cd ../
  fi
done

