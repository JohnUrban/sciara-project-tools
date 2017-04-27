#!/bin/bash


if [ -f $GATEFILE1 ] && [ -f $GATEFILE2 ]; then
  echo Both $GATEFILE1 and $GATEFILE2 exist....
  LC=`cat $GATEFILE1 | wc -l`
  LC2=`cat $GATEFILE2 | wc -l`
  if [ $LC -gt 0 ] && [ $LC2 -gt 0 ]; then
    echo Gate files are non-empty so deleting delfile, $DELFILE
    echo $GATEFILE1 $LC
    echo $GATEFILE2 $LC2
    rm $DELFILE
  else
    echo Gate-files are empty so keeping delfile, $DELFILE
    echo $GATEFILE1 $LC
    echo $GATEFILE2 $LC2
  fi
else
  echo Gate files (${GATEFILE1} and ${GATEFILE2}) do not exist so not deleting $DELFILE
fi
