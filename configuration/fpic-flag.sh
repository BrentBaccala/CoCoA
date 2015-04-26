#!/bin/bash

# Auxiliary script for CoCoALib configuration process.

# Script to see whether the -fPIC flag produces annoying compiler warnings.
# If no warning is produced, the script prints -fPIC; otherwise it prints nothing.

if [ $# -ne 1 ]
then
  echo "***ERROR*** $0 needs 1 arg (name of C++ compiler)"
  exit 1
fi

CXX="$1"

FPIC_FLAG=-fPIC

umask 22
TMP_DIR=fpic-check-$HOSTNAME-$UID-$$
/bin/rm -rf "$TMP_DIR" && mkdir "$TMP_DIR"
if [ $? -ne 0 ]
then
  echo "***ERROR*** $0: unable to create temp directory $TMP_DIR"  
  exit 2
fi
cd "$TMP_DIR"

cat > test.C <<EOF
int f(int x)
{
  return (x+1)*x+41;
}
EOF


COMPILER_MESG=`"$CXX" $FPIC_FLAG -c -o test.o test.C 2>& 1`
if [ $? -ne 0 ]
then
  echo "***ERROR*** $0: test compilation failed"
  exit 3
fi

cd ..
/bin/rm -rf "$TMP_DIR"
if [ -z "$COMPILER_MESG" ]
then
  echo $FPIC_FLAG
fi
