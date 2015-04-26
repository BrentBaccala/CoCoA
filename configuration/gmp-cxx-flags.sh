#! /bin/bash

# This script prints out the -m*** flags used for compiling GMP
# (or exits with code 1 if an error occurred).

if [ $# != 0 ]
then
  echo "ERROR: $0 expects no arguments"
  exit 1
fi

EXTLIBS=configuration/ExternalLibs
if [ \! -d configuration -o \! -d $EXTLIBS -o \! -d $EXTLIBS/include ]
then
  echo "ERROR: $0 expects the $EXTLIBS/ subtree to exist"
  exit 1
fi

# if [ "$#" != "1" ]
# then
#   echo "ERROR: $0 requires 1 arg (full path of gmp.h)" > /dev/stderr
#   exit 1
# fi

# # Assume the file exists and is readable -- already checked in configure script.
# GMP_H="$1"
# GMP_INC_DIR=`dirname "$GMP_H"`

# Below we create a small C++ program for printing out the GMP compilation flags.
# We use two temporary files.
TMP_DIR=configuration/get-gmp-cxxflags-$UID@$HOSTNAME-$$
/bin/rm -rf $TMP_DIR
mkdir $TMP_DIR
pushd $TMP_DIR  >/dev/null
cat > prog.C <<EOF
#include "gmp.h"
#include <iostream>

int main()
{
  std::cout << __GMP_CFLAGS << std::endl;
}
EOF

$CXX -I ../ExternalLibs/include prog.C -o prog 2> /dev/null
#$CXX -I "$GMP_INC_DIR" prog.C -o prog 2> /dev/null
# if [ $? -ne 0 ]
# then
#   cd ..
#   /bin/rm -rf $TMPDIR
#   exit 1
# fi
GMP_CXXFLAGS=`./prog`
#if [ $? -ne 0 ]
#then
#  # Deliberately leave $TMPDIR to assist debugging.
#  echo "ERROR: $0 test program crashed! "
#  exit 1
#fi
popd  > /dev/null
/bin/rm -rf $TMP_DIR

# GMP_CXXFLAGS contains all the compilation flags used when building GMP.
# We pick out just the compilation options which begin with -m
CoCoALib_CXXFLAGS=
for opt in $GMP_CXXFLAGS
do
  case $opt in
  ( -m* )
    CoCoALib_CXXFLAGS="$CoCoALib_CXXFLAGS $opt";;
  esac
done

echo $CoCoALib_CXXFLAGS
