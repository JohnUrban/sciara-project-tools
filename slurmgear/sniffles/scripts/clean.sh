#!/bin/bash


if [ -f $GATEFILE ]; then
  echo Gate file, $GATEFILE, exists....
  LC=`cat $GATEFILE | wc -l`
  if [ $LC -gt 0 ]; then
    echo Gate file is non-empty so deleting delfile, $DELFILE
    rm $DELFILE
  else
    echo Gate-file is empty so keeping delfile, $DELFILE
  fi
else
  echo Gate file, $GATEFILE, does not exist so not deleting $DELFILE
fi
