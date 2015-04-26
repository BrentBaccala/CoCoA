#!/bin/bash

# This script checks compatibility of user supplied CXXFLAGS with
# the GMP library.
# Exit code is 0 iff CXXFLAGS and GMP are compatible, o/w non-zero.

# This script expects 1 arg: the CXXFLAGS used to compile GMP
# Also expects that CXX and CXXFLAGS have been exported from caller.
# Also expects that the ExternalLibs/ subtree has been created.

if [ $# != 1 ]
then
  echo "$0: ERROR: requires 1 arg (CXXFLAGS_FOR_GMP)"
  exit 1
fi

CXXFLAGS_FOR_GMP="$1"
GMP_LDLIB=-lgmp
if [ -f ExternalLibs/lib/libgmp-symlink.a ]
then
  GMP_LDLIB=-lgmp-symlink
fi

umask 22

# Create tmp directory, put C prog in it, compile and run.
# TMP_DIR depends on hostname, userid, and process number to avoid unfortunate
# name clashes if several people try to install CoCoALib simultaneously.
TMP_DIR=gmp-check-cxxflags-$HOSTNAME-$$
/bin/rm -rf $TMP_DIR && mkdir $TMP_DIR
if [ $? -ne 0 ]; then
  echo "ERROR: $0 failed to create temporary directory \"$TMP_DIR\""
  exit 1
fi

cd $TMP_DIR
cat > TestProg.c <<EOF
#include "gmp.h"

int main()
{
  mpz_t A; mpz_init_set_ui(A,1);
  mpz_t B; mpz_init_set_ui(B,1);
  mpz_add(A, A, B);
  if (mpz_cmp_ui(A, 2) != 0)
    return 1;
  return 0;
}
EOF


# Compile test prog using given CXX and CXXFLAGS, and GMP header and library
$CXX $CXXFLAGS $CXXFLAGS_FOR_GMP TestProg.c -o TestProg -I ../ExternalLibs/include -L../ExternalLibs/lib $GMP_LDLIB 2> /dev/null

# Check whether compilation failed; if so, complain.
if [ $? -ne 0 ]
then
  echo "ERROR: $0 failed to compile/link TestProg; command used was"
  echo "ERROR: $CXX $CXXFLAGS $CXXFLAGS_FOR_GMP TestProg.c -o TestProg -I ../ExternalLibs/include -L../ExternalLibs/lib $GMP_LDLIB"
#  cd ..
#  /bin/rm -rf $TMP_DIR
  exit 2
fi

# Compilation succeeded, so try running $PROG.
./TestProg  2>/dev/null

# Check whether execution failed; if so, complain (probably linker problems).
if [ $? -ne 0 ]
then
  echo "ERROR: $0 TestProg crashed"
#  cd ..
#  /bin/rm -rf $TMP_DIR
  exit 3
fi

cd ..
/bin/rm -rf $TMP_DIR
exit 0
