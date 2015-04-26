#!/bin/bash

# This script checks prints out the extension of a given filename
# (it is vaguely reminiscent of the basename command)

if [ "$#" != "1" ]
then
  echo "$0: expecting 1 arg (filename to split into base and extn)"
  exit 1
fi

file=`basename "$1"`

suffix=`echo "$file" | cut -f 2- -d "."`
if [ "$file" = "$suffix" ]
then
  exit 0
fi

while true
do
  prev="$suffix"
  suffix=`echo "$suffix" | cut -f 2- -d "."`
  if [ "$prev" = "$suffix" ]
  then
    break
  fi
done
echo "$suffix"
