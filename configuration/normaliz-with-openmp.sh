#! /bin/sh

if [ "$#" != "1" ]
then
  echo "ERROR: $0 requires 1 arg (full path of libnormaliz.a)" > /dev/stderr
  exit 1
fi

nm "$1" | fgrep -l " omp" > /dev/null
if [ $? -ne 1 ]
then
  echo "OpenMP"
fi
