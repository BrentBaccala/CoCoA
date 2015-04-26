#! /bin/sh

# This script tests whether CXX can compile with gmp (and gmpxx)
# by just adding -lgmp (and -lgmpxx) as flags.
# The check entails compiling very simple source files.

# Exit code is non-zero if "-lgmp" does not work (or an error occurs).
# Exit code is 0 if "-lgmp" works; and first part of output is GMP version,
# if "-lgmpxx" works as well, the the output contains just 1 part, otherwise
# if "-lgmpxx" does not work then there is a second part.

if [ -z "$CXX" ]
then
  echo "$0: ERROR: environment variable CXX not set."
  exit 1
fi

# Now create a temp dir in which to work:
# I use a temp dir so all junk created by the compiler can be removed
# simply by removing the directory -- also I do not need to worry about
# stomping on existing files.
if [ \! -w . ]
then
  echo "$0: ERROR: current directory must be writable."
  exit 1
fi

umask 22
TMP_DIR=gmp-try-default-$HOSTNAME-$UID-$$
/bin/rm -rf $TMP_DIR  &&  mkdir $TMP_DIR  2>/dev/null
if [ $? -ne 0 ]
then
  echo "ERROR: Unable to create working directory ($TMP_DIR)"
  exit 1
fi

cd $TMP_DIR

# Here is the simple source code we shall use to test for gmp:
cat > TestGMP.C <<EOF
#include "gmp.h"
#include <iostream>

int main()
{
  mpz_t i;
  mpz_init(i);
  mpz_set_si(i, 12345);
  std::cout << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << std::endl;
}
EOF

# Try plain compiler (without CXXFLAGS):
$CXX $CXXFLAGS TestGMP.C -o TestGMP -lgmp > /dev/null 2>&1
if [ $? -ne 0 -o \! -f TestGMP -o \! -x TestGMP ]
then
  cd ..; /bin/rm -rf $TMP_DIR
  echo "No default gmp installation"
  exit 2
fi

GMP_VER=`./TestGMP`
echo $GMP_VER

# Here is the simple source code we shall use to test for gmpxx:
cat > TestGMPXX.C <<EOF
#include "gmpxx.h"
#include <iostream>

int main()
{
  mpz_class i(12345);
  std::cout << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << std::endl;
}
EOF

# Try plain compiler (without CXXFLAGS):
$CXX $CXXFLAGS TestGMPXX.C -o TestGMPXX -lgmpxx -lgmp > /dev/null 2>&1
if [ $? -ne 0 -o \! -f TestGMPXX -o \! -x TestGMPXX ]
then
  cd ..; /bin/rm -rf $TMP_DIR
  echo "No default gmpxx installation"
  exit 0
fi

cd ..; /bin/rm -rf $TMP_DIR
